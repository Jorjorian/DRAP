version development

task bwaMem {
    input {
        Array[File] reference
        File fq1
        File fq2
        String tag
        String platform = "illumina"
        Int cpu = 3
        String docker = "gcr.io/clearlabs-science/janus_test"
    }
    command <<<
        source /opt/conda/etc/profile.d/conda.sh
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
    }
    String infix = if tag == "" then tag else "." + tag
    String prefix = basename(bam[0], ".bam") + infix
    command <<<
        conda activate sorbet_env

        # Initialize variables for JSON data
        echo "{" > individual_data.json
        samtools stats ~{bam[0]} > ~{prefix + ".stats"}

        n_mapped_reads=$(grep 'reads mapped:' ~{prefix + ".stats"} | cut -f 3)
        echo $n_mapped_reads > mapped_reads
        if [ ~{read_counts} -gt 0 ]; then
        mapping_rate=$(bc -l <<< "(100 * $n_mapped_reads / ( 2 * ~{read_counts})) ")
        echo $mapping_rate > mapping_rate
        else
        awk -F"\t" '/^SN\treads mapped:/ {mapped=$3} /^SN\traw total sequences:/ {total=$3} END {print (mapped/total)*100}' ~{prefix + ".stats"} > mapping_rate
        fi

        samtools depth -a ~{bam[0]} > ~{prefix + ".all.depth"}

        sed -n '/^SN\tinsert size average:/{p;q}' ~{prefix + ".stats"} |cut -f3 > insert_size_mean
        sed -n '/^SN\tinsert size standard deviation:/{p;q}' ~{prefix + ".stats"} |cut -f3 > insert_size_sd

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

        cat ~{prefix + ".all.depth"} |datamash median 3 > median_depth
        cat ~{prefix + ".all.depth"} |datamash mean 3 > mean_depth

        # Finalize JSON file
        echo "\"average_bamstats\": {\"depth\": $(<mean_depth), \"median_depth\": $(<median_depth), \"insert_size_mean\": $(<insert_size_mean), \"insert_size_sd\": $(<insert_size_sd)}" >> individual_data.json
        echo "}" >> individual_data.json
    >>>
    output {
        String sample = tag
        File stats = prefix + ".stats"
        File all_depth = prefix + ".all.depth"
        Float insert_size_mean = if size("insert_size_mean") > 1 then read_float("insert_size_mean") else -1.0
        Float insert_size_sd = if size("insert_size_sd") > 1 then read_float("insert_size_sd") else -1.0
        Float mapping_rate = if size("mapping_rate") > 1 then read_float("mapping_rate") else -1.0
        Float median_depth = if size("median_depth") > 1 then read_float("median_depth") else -1.0
        Float mean_depth = if size("mean_depth") > 1 then read_float("mean_depth") else -1.0
        File bam_json_full = write_json({
                                            "tag": tag,
                                            "% Mapped Reads": mapping_rate,
                                            "insert_size_mean": insert_size_mean,
                                            "insert_size_sd": insert_size_sd,
                                            "median_depth": median_depth,
                                            "mean_depth": mean_depth
                                        })
        File bam_json_empty = write_json({"tag": tag,
                                             "% Mapped Reads": -1,
                                             "insert_size_mean": -1,
                                             "insert_size_sd": -1,
                                             "median_depth": -1,
                                             "mean_depth": -1
                                         })
        File bam_json = if empty then bam_json_empty else bam_json_full
        File overall_json = 'individual_data.json'
    }
    runtime {
        memory: "2 GB"
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
            tag = tag,
            reference = ref
    }
    output {
        File bam_json = bamStat.bam_json
    }
}

