version development

task dimer_metrics {
    input {
        File fq_1
        File fq_2
        File adapter_list
        File? primer_list
        String tag
    }

    command <<<

        # Create a directory for output
        mkdir -p ~{tag}_output

        # Run BBDuk to trim adapters and primers and get counts
        bbduk.sh \
        in1=~{fq_1} in2=~{fq_2} \
        ref=~{adapter_list} ~{', ' + primer_list} \
        stats=~{tag}_output/~{tag}.bbmap.stats \
        out1=~{tag}_output/filtered_1.fastq out2=${tag}_output/filtered_2.fastq


        # Convert dimer counts to JSON format

        python <<CODE
import json
import os
dimer_counts = {}
if os.path.exists("~{tag}_output/~{tag}.bbmap.stats"):
    with open("~{tag}_output/~{tag}.bbmap.stats") as f:
        for line in f:
            if line.startswith("#"):
                continue
            adapter, count, perc = line.strip().split()
            dimer_counts[adapter] = {"count": int(count), "percentage": float(perc.strip("%"))}
else:
    dimer_counts = {}

with open('~{tag}_dimer_metrics.json', 'w') as f:
    json.dump(dimer_counts, f, indent=4)
# strip first 3 lines from stats file
if os.path.exists("~{tag}_output/~{tag}.bbmap.stats"):
    with open("~{tag}_output/~{tag}.bbmap.stats") as f:
        lines = f.readlines()
        lines = lines[3:]
else:
    lines = []
with open("~{tag}_output/~{tag}.bbmap.stats", "w") as f:
    f.writelines(lines)
CODE
    >>>

    output {
        File dimer_metrics_json = "~{tag}_dimer_metrics.json"
        File bbduk_stats = "~{tag}_output/~{tag}.bbmap.stats"
    }
    runtime {
        memory: "10 GB"
        cpu: 1
        docker: "docker.io/oblivious1/drap:bbmap"
    }
}

task combine_dimer_metrics {
    input {
        Array[File] dimer_metrics_files
    }

    command <<<

    # Combine dimer metric JSONs into a single MultiQC-compatible JSON format

    python <<CODE
import json
import os
combined_data = {}

for dimer_file in "~{sep=' ' dimer_metrics_files}".split(' '):
  with open(dimer_file) as f:
     dimer_data = json.load(f)
  sample_name = os.path.basename(dimer_file).split('.')[0].split('_')[0]
  combined_data[sample_name] = {}
  for adapter, metrics in dimer_data.items():
     combined_data[sample_name][adapter] = metrics["count"]
     multiqc_data = {
          "id": "custom_data_bargraph",
          "section_name": "Dimer Metrics",
          "description": "This plot shows dimer metrics from multiple samples.",
          "plot_type": "bargraph",
          "pconfig": {
                "id": "dimer_metrics_bargraph",
                "title": "Dimer Metrics",
                "ylab": "Count",
                "xlab": "Adapters",
                "xDecimals": False },
          "data": combined_data
         }

with open("dimermetrics_mqc.json", 'w') as f:
  json.dump(multiqc_data, f, indent=4)
CODE
>>>

output {
  File multiqc_dimer_metrics_json = "dimermetrics_mqc.json"
 }

runtime {
    memory: "4 GB"
    cpu: 1
    docker: "python:3.8-slim"
    }
}
