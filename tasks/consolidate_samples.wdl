version development

task consolidate_samples {
    input {
        File zip_file
        Array[String] group_on
        File? default_reference
    }

    command <<<
        set -e  # Exit on error

        OUTPUT_DIR=$(pwd)/output_fastqs
        mkdir -p ${OUTPUT_DIR}
        unzip ~{zip_file} -d extracted_contents
        cd $(find extracted_contents -type d | sort | tail -n 1)
        for group in ~{sep=' ' group_on}; do
        R1_files=$(find . \( -name "*${group}*R1*.fastq.gz" -o -name "*${group}*1.fq.gz" \) | tr '\n' ' ')
        R2_files=$(find . \( -name "*${group}*R2*.fastq.gz" -o -name "*${group}*2.fq.gz" \) | tr '\n' ' ')

            # Concatenate the fastq files
            cat ${R1_files} > ${OUTPUT_DIR}/${group}_R1.fastq.gz
            cat ${R2_files} > ${OUTPUT_DIR}/${group}_R2.fastq.gz

            # Save the group name to the names.txt file
        done
        ls ${OUTPUT_DIR}/*_R1.fastq.gz | sed 's/.*\///' | sed 's/_R1.fastq.gz//' > ${OUTPUT_DIR}/names.txt

        # Concatenate any fasta files
        if [ -n "~{default_reference}" ]; then
            cp ~{default_reference} ${OUTPUT_DIR}/combined.fasta
        else
            touch ${OUTPUT_DIR}/combined.fasta
        fi

        for fasta_file in $(find . -name "*.fasta"); do
            cat ${fasta_file} >> ${OUTPUT_DIR}/combined.fasta
        done
    >>>

    output {
        Array[File] R1_files = glob("output_fastqs/*_R1.fastq.gz")
        Array[File] R2_files = glob("output_fastqs/*_R2.fastq.gz")
        Array[String] names = read_lines("output_fastqs/names.txt")
        File reference = "output_fastqs/combined.fasta"

        Array[Pair[File, File]] fastq_pairs = zip(R1_files, R2_files)
        Array[Pair[String, Pair[File, File]]] output_fastq_pairs = zip(names, fastq_pairs)
    }

    runtime {
        docker: "docker.io/oblivious1/drap:bash"
        memory: "4G"
        cpu: 2
    }
}



