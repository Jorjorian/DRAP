version 1.0

task fastqStage {
    input {
        File  source_fq1
        File? source_fq2
        String prefix
    }
    command <<<
        cp ~{source_fq1} ~{prefix + "_1.fq.gz"}
        if [ -r "~{source_fq2}" ]; then
            cp ~{source_fq2} ~{prefix + "_2.fq.gz"}
        fi

    >>>
    output {
        File  fq1 = "${prefix}_1.fq.gz"
        File? fq2 = "${prefix}_2.fq.gz"
    }
    runtime {
        docker: "gcr.io/clearlabs-science/janus_test"
    }
}
task passing_filter {
    input {
        Array[File] bcl2fastq_demux_stats
        Array[File] bcl2fastq_fastq_stats
        Array[String] sequencers_used
        String sample_id
        String docker = "gcr.io/clearlabs-science/janus_test"
        Boolean empty = false
    }

    command <<<
    python <<CODE
    true = True
    false = False
    import csv

    # Input values and files
    bcl2fastq_demux_stats_files = ["~{sep='","' bcl2fastq_demux_stats}"]
    bcl2fastq_fastq_stats_files = ["~{sep='","' bcl2fastq_fastq_stats}"]
    sequencers_used = [~{sep=',' sequencers_used}]
    sample_id = "~{sample_id}"

    # Map sample_id to SampleNumber
    sample_number = None
    mapping_list = []
    if ~{empty} == True:
        with open("pf_reads", 'w') as f:
            f.write(str(-1))
        with open("percent_reads", 'w') as f:
            f.write(str(-1.0))
        exit(0)
    for demux_file in bcl2fastq_demux_stats_files:
        with open(demux_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            i =0
            row_list = []
            for row in reader:
                row_list.append(row)
                i+=1
                if i > 1:
                    break
        mapping_list.append({row_list[1][i]: row_list[0][i] for i in range(len(row_list[0]))})
    for fastq_file in bcl2fastq_fastq_stats_files:
        with open(fastq_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader)  # skip header
            for row in reader:
                if row[0] == sample_id:
                    sample_number = int(row[0])
                    break
    # Extract and sum the NumberOfReadsPF values for the given sample_id
    pf_reads_for_sample = 0
    total_pf_reads = 0
    for i, fastq_file in enumerate(bcl2fastq_fastq_stats_files):
        if not sequencers_used[i]:
            continue
        with open(fastq_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader)  # skip header
            for row in reader:
                print(row[3])
                total_pf_reads += int(row[3])
                if (int(row[0]) == int(mapping_list[i][sample_id])) and sequencers_used[i]:
                    pf_reads_for_sample += int(row[3])

    # Calculate % Reads Identified (PF)
    percent_reads = (pf_reads_for_sample / total_pf_reads) * 100 if total_pf_reads > 0 else 0

    # Write outputs to files
    with open("pf_reads", 'w') as f:
        f.write(str(pf_reads_for_sample))

    with open("percent_reads", 'w') as f:
        f.write(str(percent_reads))

    CODE

    >>>
    output {
        Int pf_reads = if size("pf_reads") > 0 then read_int("pf_reads") else -1
        Float percent_reads = if size("percent_reads") > 0 then read_float("percent_reads") else -1.0
        File pf_json = write_json({
        "sample_index": sample_id,
        "PF Reads for Sample": pf_reads,
        "% Reads Identified (PF)": percent_reads
        })

    }
    runtime {
        docker: "~{docker}"
    }
}

task fastaIndex {
    input {
        File fasta
        String docker = "gcr.io/clearlabs-science/janus_test"
    }
    String bsname = basename(fasta)
    String prefix = sub(bsname, "\.[^.]+$", "")

    command <<<
        conda activate sorbet_env
        cp ~{fasta} ~{bsname}
        bwa index ~{bsname}
        samtools faidx ~{bsname}
        samtools dict -o ~{prefix + ".dict"} ~{bsname}
    >>>
    output {
        Array[File] index = [
            bsname, bsname + ".fai", prefix + ".dict"
        ]
        Array[File] bwa_index = flatten([
            [bsname],
            prefix(bsname, [".amb", ".ann", ".bwt", ".pac", ".sa"])
        ])
        Array[File] all_index = flatten([
            [bsname],
            prefix(bsname, [".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"]),
            [prefix + ".dict"]
        ])
    }
    runtime {
        memory: "2 GB"
        docker: "~{docker}"
    }
}
task fastqDump {
    input {
        String sra
    }

    command <<<

        # Prefetch the SRA data
        prefetch ~{sra}

        # Split into paired-end fastq files
        fastq-dump --split-files ~{sra} --gzip
    >>>

    output {
        File fq1 = "~{sra}_1.fastq.gz"
        File fq2 = "~{sra}_2.fastq.gz"
    }

    runtime {
        docker: "gcr.io/clearlabs-science/sra_bash"
        memory: "5 GB"
        cpu: 1
        maxRetries: 4
    }
}
# task startServer {
#     input {}
#     command <<<
#         curl -v localhost:8000/engine/v1/status 2>&1 |
#             grep '200 OK' > status
#         if [ ! -s "status" ]; then
#             nohup java -Dbackend.providers.Local.config.concurrent-job-limit=2 \
#                 -jar /opt/cromwell.jar server &
#         fi
#     >>>
#     output {
#         File? log = "nohup.out"
#     }
#     runtime {
#         memory: "1 GB"
#     }
# }

# task submitWorkflow {
#     input {
#         File wdl
#         File param
#     }
#     command <<<
#         java -jar /opt/cromwell.jar submit -t WDL -v 1.0 \
#             ~{wdl} -i ~{param}
#     >>>
#     runtime {
#         memory: "1 GB"
#     }
# }
