version development

task consolidate_samples {
    input {
        Array[Pair[String, Array[Pair[File, File]]]] dropped_read_groups
    }

    # Convert the dropped_read_groups to a JSON string
    File dropped_read_groups_json = write_json(dropped_read_groups)

    command <<<
        set -e  # Exit on error

        mkdir -p output_fastqs

        # Write the dropped_read_groups JSON string to a file safely
        cp ~{dropped_read_groups_json} dropped_read_groups.json

        # Process the JSON file to consolidate samples
        jq -c '.[]' dropped_read_groups.json | while read pair; do
        name=$(echo ${pair} | jq -r '.left')
        files=$(echo ${pair} | jq -c '.right')

        # Extract R1 and R2 files
        R1_files=$(echo ${files} | jq -r '.[] | .left' | tr '\n' ' ')
        echo ${R1_files}

        R2_files=$(echo ${files} | jq -r '.[] | .right' | tr '\n' ' ')
        echo ${R2_files}

        # Concatenate the files into the output directory
        cat ${R1_files} > output_fastqs/${name}_R1.fastq.gz
        cat ${R2_files} > output_fastqs/${name}_R2.fastq.gz

        # Save the name to the names.txt file
        echo "${name}" >> output_fastqs/names.txt
        done
    >>>

    output {
        Array[File] R1_files = glob("output_fastqs/*_R1.fastq.gz")
        Array[File] R2_files = glob("output_fastqs/*_R2.fastq.gz")
        Array[String] names = read_lines("output_fastqs/names.txt")

        Array[Pair[File, File]] fastq_pairs = zip(R1_files, R2_files)
        Array[Pair[String, Pair[File, File]]] output_fastq_pairs = zip(names, fastq_pairs)
    }

    runtime {
        docker: "docker.io/oblivious1/drap:bash"
        memory: "4G"
        cpu: 2
    }
}



