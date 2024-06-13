version development

task consolidate_samples {
    input {
    Array[Pair[String, Array[Pair[File, File]]]] dropped_read_groups
}

command <<<
    mkdir -p output_fastqs
    for pair in ${write_lines(select_first([dropped_read_groups]))}; do
    name=$(echo ${pair} | cut -d' ' -f1)
    files=$(echo ${pair} | cut -d' ' -f2)
    R1_files=$(echo ${files} | jq -r '.[] | .left')
    R2_files=$(echo ${files} | jq -r '.[] | .right')

    cat ${R1_files} > output_fastqs/${name}_R1.fastq
    cat ${R2_files} > output_fastqs/${name}_R2.fastq
    done
>>>

output {
    Array[Pair[File, File]] output_fastq_pairs =

    }
runtime {
    docker: "ubuntu:latest"
                    memory: "4G"
    cpu: 2
}
}