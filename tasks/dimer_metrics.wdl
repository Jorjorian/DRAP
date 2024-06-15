version development

task dimer_metrics {
    input {
        File fq_1
        File fq_2
        File adapter_list
        File primer_list
        String tag
    }

    command <<<
        # Create a directory for output
        mkdir -p ${tag}_output

        # Run BBDuk to trim adapters and primers and get counts
        bbduk.sh \
        in1=${fq_1} in2=${fq_2} \
        ref=${adapter_list},${primer_list} \
        ktrim=r k=23 mink=11 hdist=1 \
        stats=${tag}_output/bbduk_stats.txt \
        out1=${tag}_output/trimmed_1.fastq out2=${tag}_output/trimmed_2.fastq

        # Parse the BBDuk stats to extract dimer counts
        awk 'BEGIN {FS="\t"}; NR>1 {print $1, $2}' ${tag}_output/bbduk_stats.txt > ${tag}_output/dimer_counts.txt

        # Convert dimer counts to JSON format
        python <<CODE
import json

dimer_counts = {}
with open('${tag}_output/dimer_counts.txt') as f:
for line in f:
    adapter, count = line.strip().split()
    dimer_counts[adapter] = int(count)

with open('${tag}_dimer_metrics.json', 'w') as f:
    json.dump(dimer_counts, f, indent=4)
CODE
    >>>

    output {
        File dimer_metrics_json = "${tag}_dimer_metrics.json"
    }
    runtime {
        memory: "10 GB"
        cpu: 4
        docker: "quay.io/biocontainers/bbmap:39.06--h92535d8_1"
    }
}