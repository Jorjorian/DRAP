version development

task bwaMem {
    input {
        Array[File] reference
        File fq1
        File fq2
        String tag
        String platform = "illumina"
        Int cpu = 8
        String docker = "docker.io/oblivious1/drap:bwa"
    }
    command <<<
        source activate base
        conda activate bwa_env
        bwa mem -t ~{cpu} \
            -R "~{'@RG\\tID:' + tag + '\\tSM:' + tag + '\\tLB:' + tag + '\\tPL:' + platform}" \
            ~{reference[0]}  \
            -k 12 -r 1 -B 2 -T 20 -O 3 \
            ~{fq1} ~{fq2} |
            samtools sort --reference ~{reference[0]} \
                --threads ~{cpu} \
                -o ~{tag + '.bam'}
        samtools index ~{tag + '.bam'}
    >>>

    output {
        Array[File] bam = [
          tag + ".bam",
          tag + ".bam.bai"
        ]
    }

    runtime {
        memory: "6 GB"
        cpu: cpu
        docker: "~{docker}"
    }
}

task bamStat {
    input {
        Array[File] bam
        String tag = ""
        Boolean empty = false
        Int read_counts = 0
        String docker = "docker.io/oblivious1/drap:bwa"
    }
    String infix = if tag == "" then tag else "." + tag
    command <<<
        source activate base
        conda activate bwa_env
        # Initialize variables for JSON data
        echo "{" > individual_data.json
        samtools stats ~{bam[0]} > ~{tag + ".stats"}

        n_mapped_reads=$(grep 'reads mapped:' ~{tag + ".stats"} | cut -f 3)
        total_reads=$(grep 'raw total sequences:' ~{tag + ".stats"} | cut -f 3)
        error_rate=$(grep 'error rate:' ~{tag + ".stats"} | cut -f 3)
        average_length=$(grep 'average length:' ~{tag + ".stats"} | cut -f 3)
        echo $total_reads > total_reads
        echo $n_mapped_reads > mapped_reads
        echo $total_reads > total_reads
        echo $error_rate > error_rate
        echo $average_length > average_length
        if [ ~{read_counts} -gt 0 ]; then
        mapping_rate=$(bc -l <<< "(100 * $n_mapped_reads / ( 2 * ~{read_counts})) ")
        echo $mapping_rate > mapping_rate
        else
        awk -F"\t" '/^SN\treads mapped:/ {mapped=$3} /^SN\traw total sequences:/ {total=$3} END {print (mapped/total)*100}' ~{tag + ".stats"} > mapping_rate
        fi

        samtools depth -a ~{bam[0]} > ~{tag + ".all.depth"}

        sed -n '/^SN\tinsert size average:/{p;q}' ~{tag + ".stats"} |cut -f3 > insert_size_mean
        sed -n '/^SN\tinsert size standard deviation:/{p;q}' ~{tag + ".stats"} |cut -f3 > insert_size_sd

        # Python script to calculate median insert size
        python <<CODE
        import pandas as pd
        import numpy as np

        def weighted_median(values, weights, quantiles=0.5, interpolate=False):
        i = values.argsort()
        sorted_weights = weights[i]
        sorted_values = values[i]
        Sn = sorted_weights.cumsum()
        if interpolate:
        Pn = (Sn - sorted_weights/2 ) / Sn[-1]
        return np.interp(quantiles, Pn, sorted_values)
        else:
        return sorted_values[np.searchsorted(Sn, quantiles * Sn[-1])]

        # Read the file into a DataFrame
        df = pd.read_csv('insert_sizes_distribution.txt', sep='\t', header=None, names=['size', 'count'])
        # Calculate the median insert size
        if not df.empty:
        median_insert_size = weighted_median(df['size'].values ,df['count'].values)
        else:
        median_insert_size = ''
        # Print the median
        with open('insert_size_median', 'w') as o:
        o.write(str(median_insert_size))
        CODE

        cat ~{tag + ".all.depth"} |datamash median 3 > median_depth
        cat ~{tag + ".all.depth"} |datamash mean 3 > mean_depth
    >>>
    output {
        String sample = tag
        File stats = tag + ".stats"
        File all_depth = tag + ".all.depth"
        Float insert_size_mean = if size("insert_size_mean") > 1 then read_float("insert_size_mean") else -1.0
        Float insert_size_sd = if size("insert_size_sd") > 1 then read_float("insert_size_sd") else -1.0
        Float mapping_rate = if size("mapping_rate") > 1 then read_float("mapping_rate") else -1.0
        Float median_depth = if size("median_depth") > 1 then read_float("median_depth") else -1.0
        Float mean_depth = if size("mean_depth") > 1 then read_float("mean_depth") else -1.0
        Float average_length = if size("average_length") > 1 then read_float("average_length") else -1.0
        Float error_rate = if size("error_rate") > 1 then read_float("error_rate") else -1.0
        Float n_mapped_reads = if size("mapped_reads") > 1 then read_float("mapped_reads") else -1.0
        Float total_reads = if size("total_reads") > 1 then read_float("total_reads") else -1.0
        File bam_json_full = write_json({
                                            "tag": tag,
                                            "% Mapped Reads": mapping_rate,
                                             "average_length": average_length,
                                            "error_rate": error_rate,
                                            "n_mapped_reads": n_mapped_reads,
                                            "total_reads": total_reads,
                                            "insert_size_mean": insert_size_mean,
                                            "insert_size_sd": insert_size_sd,
                                            "median_depth": median_depth,
                                            "mean_depth": mean_depth
                                        })
        File bam_json_empty = write_json({"tag": tag,
                                             "% Mapped Reads": -1,
                                             "average_length": -1,
                                             "error_rate": -1,
                                             "n_mapped_reads": -1,
                                             "total_reads": -1,
                                             "insert_size_mean": -1,
                                             "insert_size_sd": -1,
                                             "median_depth": -1,
                                             "mean_depth": -1
                                         })
        File bam_json = if empty then bam_json_empty else bam_json_full
        File overall_json = 'individual_data.json'
    }
    runtime {
        memory: "8 GB"
        docker: "~{docker}"
    }
}

workflow loose_mapping_metrics {
    input {
        File fq1
        File fq2
        String tag
        Array[File] ref
    }
    call bwaMem {
        input:
            fq1 = fq1,
            fq2 = fq2,
            reference = ref,
            tag = tag
    }
    call bamStat {
        input:
            bam = bwaMem.bam,
            tag = tag
    }
    output {
        File bam_json = bamStat.bam_json
        File stats = bamStat.stats
        File depth = bamStat.all_depth
    }
}

