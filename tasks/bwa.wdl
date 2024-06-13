version development

task bwaMem {
    input {
        Array[File] reference
        File fq1
        File fq2
        String sample = basename(fq1, "_1.fq.gz")
        String platform = "illumina"
        Int cpu = 3
        String docker = "gcr.io/clearlabs-science/janus_test"
    }
    command <<<
        conda activate sorbet_env
        # -k Minimum seed length 20 -> 12
        # -r Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. This is a key heuristic parameter for tuning the performance. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy
        # -B Mismatch penalty. 4 -> 2
        # -T Don’t output alignment with score lower than INT. 30 -> 20
        # -O Gap open penalty 6 -> 3
        conda activate sorbet_env
        bwa mem -t ~{cpu} \
            -R "~{'@RG\\tID:' + sample + '\\tSM:' + sample + '\\tLB:' + sample + '\\tPL:' + platform}" \
            ~{reference[0]}  \
            -k 12 -r 1 -B 2 -T 20 -O 3 \
            ~{fq1} ~{fq2} |
            samtools sort --reference ~{reference[0]} \
                --threads ~{cpu} \
                -o ~{sample + '.bam'}
        samtools index ~{sample + '.bam'}
    >>>

    output {
        Array[File] bam = [
            sample + ".bam",
            sample + ".bam.bai"
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
        Array[File]? reference
        String tag = ""
        String sample_id = ""
        String docker = "gcr.io/clearlabs-science/wastewater"
        Boolean empty = false
        Int read_counts = 0
    }
    String infix = if tag == "" then tag else "." + tag
    String prefix = basename(bam[0], ".bam") + infix
    command <<<
        conda activate sorbet_env
        # Extract reference sequences from BAM file

        # Initialize variables for JSON data
        echo "{" > individual_data.json
        total_depth=0
        total_coverage=0
        total_coverage_10=0
        total_coverage_100=0
        count=0
        samtools stats -c 10,10000,10 -x -p ~{bam[0]} > ~{prefix + ".stats"}
        samtools depth -as ~{bam[0]} > ~{prefix + ".effective.depth"}
        samtools idxstats ~{bam[0]} | grep -v '\*' | cut -f1 > reference_sequences.txt
        echo "{" > individual_data.json
        echo "\"sample_index\": \"~{sample_id}\"," >> individual_data.json
        reference_seq_count=`cat reference_sequences.txt | wc -l`
        if [ $reference_seq_count -lt 1000 ]; then
            for seq in $(cat reference_sequences.txt); do
                echo $seq
                grep $seq ~{prefix + ".effective.depth"} > ${seq}.depth
                # Calculate depth for the sequence
                echo `cat ${seq}.depth |wc -l` \
                `cat ${seq}.depth |awk '$3 >= 1' |wc -l` \
                `cat ${seq}.depth |awk '$3 >= 10' |wc -l` \
                `cat ${seq}.depth |awk '$3 >= 100' |wc -l` \
                > _tmp_depth
                cat _tmp_depth
                ln_count=$(cat _tmp_depth |awk '{print $1}')
                if [ $ln_count -gt 0 ]; then
                    coverage_value=$(awk '{print ($2 / $1) * 100}' _tmp_depth | awk '{printf "%f", $0}')
                    coverage_value_10=$(awk '{print ($3 / $1) * 100}' _tmp_depth | awk '{printf "%f", $0}')
                    coverage_value_100=$(awk '{print ($4 / $1) * 100}' _tmp_depth | awk '{printf "%f", $0}')
                    depth_value=$(cat tmp.depth |datamash mean 3 | awk '{printf "%f", $0}')
                else
                    coverage_value=-1
                    coverage_value_10=-1
                    coverage_value_100=-1
                    depth_value=-1
                fi
                # Calculate coverage for the sequence
                # Add to JSON
                echo "\"${seq}_bamstats\": {\"depth\": $depth_value, \"coverage\": $coverage_value, \"coverage_value_10\": $coverage_value_10,  \"coverage_value_100\": $coverage_value_100}," >> individual_data.json
                if [[ ! $seq =~ INTERNAL_CONTROL ]]; then
                        total_depth=$(echo $total_depth + $depth_value | bc)
                        total_coverage=$(echo $total_coverage + $coverage_value | bc)
                        total_coverage_10=$(echo $total_coverage_10 + $coverage_value_10 | bc)
                        total_coverage_100=$(echo $total_coverage_100 + $coverage_value_100 | bc)
                        count=$((count+1))
                    fi
            done
        fi
        # Compute average depth and coverage
        if [ $count -gt 0 ]; then
            average_depth=$(echo $total_depth / $count | bc -l | awk '{printf "%f", $0}')
            average_coverage=$(echo $total_coverage / $count | bc -l | awk '{printf "%f", $0}')
            average_coverage_10=$(echo $total_coverage_10 / $count | bc -l | awk '{printf "%f", $0}')
            average_coverage_100=$(echo $total_coverage_100 / $count | bc -l | awk '{printf "%f", $0}')
        else
            average_depth=-1
            average_coverage=-1
            average_coverage_10=-1
            average_coverage_100=-1
        fi
        echo $average_depth > mean_depth
        echo $average_coverage > depth_ge_1
        echo $average_coverage_10 > depth_ge_10
        echo $average_coverage_100 > depth_ge_100



    # Finalize JSON file
        echo "\"average_bamstats\": {\"depth\": $average_depth, \"coverage\": $average_coverage, \"coverage_10x\": $average_coverage_10, \"coverage_100x\": $average_coverage_100}" >> individual_data.json
        echo "}" >> individual_data.json
        grep ^IS  ~{prefix + ".stats"} | cut -f 2,3   > insert_sizes_distribution.txt
        python <<CODE
import pandas as pd
import numpy as np
#def weighted_median(values, weights):
#    print(values, weights)
#    i = np.argsort(values)
#    c = np.cumsum(weights)
#    return values[i[np.searchsorted(c, 0.5 * c[-1])]]
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
print(df)
# Calculate the median insert size
if not df.empty:
    median_insert_size = weighted_median(df['size'].values ,df['count'].values)
else:
    median_insert_size = ''
# Print the median
with open('insert_size_median', 'w') as o:
    o.write(str(median_insert_size))
CODE
        n_mapped_reads=$(grep 'reads mapped:' ~{prefix + ".stats"} | cut -f 3)
        echo $n_mapped_reads > mapped_reads
        if [ ~{read_counts} -gt 0 ]; then
          mapping_rate=$(bc -l <<< "(100 * $n_mapped_reads / ( 2 * ~{read_counts})) ")
          echo $mapping_rate > mapping_rate
        else
           awk -F"\t" '/^SN\treads mapped:/ {mapped=$3} /^SN\traw total sequences:/ {total=$3} END {print (mapped/total)*100}' ~{prefix + ".stats"} > mapping_rate
        fi
        samtools coverage  ~{bam[0]} > ~{prefix + ".coverage"}

        samtools depth -a ~{bam[0]} > ~{prefix + ".all.depth"}
        cat ~{prefix + ".effective.depth"} | wc -l > reference_size
        # check if bedfile is empty
        echo $average_depth > depth
        echo $average_coverage > coverage
        sed -n '/^SN\traw total sequences:/{p;q}' ~{prefix + ".stats"} |cut -f3 > total_reads
        sed -n '/^SN\treads properly paired:/{p;q}' ~{prefix + ".stats"} |cut -f3 > proper_reads
        sed -n '/^SN\tinsert size average:/{p;q}' ~{prefix + ".stats"} |cut -f3 > insert_size_mean
        sed -n '/^SN\tinsert size standard deviation:/{p;q}' ~{prefix + ".stats"} |cut -f3 > insert_size_sd
        echo `cat ~{prefix + ".effective.depth"} |wc -l` \
            `cat ~{prefix + ".effective.depth"} |awk '$3 >= 1' |wc -l` \
            `cat ~{prefix + ".effective.depth"} |awk '$3 >= 10' |wc -l` \
            `cat ~{prefix + ".effective.depth"} |awk '$3 >= 100' |wc -l` \
            > _tmp_depth
        awk '$1 > 0{print ($2 / $1) * 100} $1 <= 0{print -1}' _tmp_depth > all_depth_ge_1
        awk '$1 > 0{print ($3 / $1) * 100} $1 <= 0{print -1}' _tmp_depth > all_depth_ge_10
        awk '$1 > 0{print ($4 / $1) * 100} $1 <= 0{print -1}' _tmp_depth > all_depth_ge_100
        cat ~{prefix + ".effective.depth"} |datamash median 3 > median_depth
        cat ~{prefix + ".effective.depth"} |datamash mean 3 > mean_depth
        cat ~{prefix + ".effective.depth"} |datamash sum 3 > effective_bases
        cat ~{prefix + ".all.depth"} |datamash sum 3 > all_bases
        paste effective_bases reference_size |awk 'FNR == 1{print $1/$2}' > sequencing_depth
        # more similar to wgs_pipeline
        #   https://github.com/clearlabs/qaptr/blob/master/wgs_pipeline/artic_wgs_pipeline_paralel.py#L240
        echo $average_depth > coverage_depth
        echo $average_depth
        # Remove unmapped reads and create a filtered BAM file
        samtools view -F 4 -b ~{bam[0]} > ~{prefix + ".filtered.bam"}

        # Create .bai index file for the filtered BAM
        samtools index ~{prefix + ".filtered.bam"}
    >>>
    output {
        String sample = sample_id
        File stats = prefix + ".stats"
        File all_depth = prefix + ".all.depth"
        File effective_depth = prefix + ".effective.depth"
#        File amplicon_cov = prefix + ".amplicon.cov"
#        File amplicon_depth = prefix + ".amplicon.depth"
        Array[File] filtered_bam = [prefix + ".filtered.bam", prefix + ".filtered.bam.bai"]
        Int total_reads = if size("total_reads") > 1 then read_int("total_reads") else -1
        Int proper_reads = if size("proper_reads") > 1 then read_int("proper_reads") else -1
        Float insert_size_mean = if size("insert_size_mean") > 1 then read_float("insert_size_mean") else -1.0
        Float insert_size_median = if size("insert_size_median") > 1 then read_float("insert_size_median") else -1.0
        Float insert_size_sd = if size("insert_size_sd") > 1 then read_float("insert_size_sd") else -1.0
        Float depth_ge_1 = if size("depth_ge_1") > 1 then read_float("depth_ge_1") else -1.0
        Float depth_ge_10 = if size("depth_ge_10") > 1 then read_float("depth_ge_10") else -1.0
        Float mapping_rate = if size("mapping_rate") > 1 then read_float("mapping_rate") else -1.0
        Float depth_ge_100 = if size("depth_ge_100") > 1 then read_float("depth_ge_100") else -1.0
        Float median_depth = if size("median_depth") > 1 then read_float("median_depth") else -1.0
        Float mean_depth = if size("mean_depth") > 1 then read_float("mean_depth") else -1.0
        Float coverage = if size("coverage") > 1 then read_float("coverage") else -1.0
        Float sequencing_depth = if size("sequencing_depth") > 1 then (if read_string("sequencing_depth") == "-nan" then -1.0 else read_float("sequencing_depth")) else -1.0
        Float coverage_depth = if size("mean_depth") > 1 then read_float("mean_depth") else -1.0
        Int effective_bases = if size("effective_bases") > 1 then read_int("effective_bases") else -1
        Int mapped_reads = if size("mapped_reads") > 1 then read_int("mapped_reads") else -1
        Int all_bases = if size("all_bases") > 1 then read_int("all_bases") else -1
        File bam_json_full = write_json({
        "sample_index": sample_id,
        "% Mapped Reads": mapping_rate,
         "mapped_reads": mapped_reads,
        "total_reads": total_reads,
        "proper_reads": proper_reads,
        "insert_size_median": insert_size_median,
        "insert_size_mean": insert_size_mean,
        "insert_size_sd": insert_size_sd,
        "Genome Coverage (>=1x)": depth_ge_1,
        "Genome Coverage (>=10x)": depth_ge_10,
        "Genome Coverage (>=100x)": depth_ge_100,
        "median_depth": median_depth,
        "mean_depth": mean_depth,
        "Sequencing Depth": sequencing_depth,
        "effective_bases": effective_bases,
        "all_bases": all_bases
    })
        File bam_json_empty = write_json({"sample_index": sample_id,
        "% Mapped Reads": -1,
        "mapped_reads": -1,
        "total_reads": -1,
        "proper_reads": -1,
        "insert_size_median": -1,
        "insert_size_mean": -1,
        "insert_size_sd": -1,
        "Genome Coverage (>=1x)": -1,
        "Genome Coverage (>=10x)": -1,
        "Genome Coverage (>=100x)": -1,
        "median_depth": -1,
        "mean_depth": -1,
        "Sequencing Depth": -1,
        "effective_bases": -1,
        "all_bases": -1
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
            sample = tag
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

