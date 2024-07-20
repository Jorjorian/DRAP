version development
task multiQC {
    input {
    Array[Array[File?]] misc_file_arrays
    Array[File?] misc_files
    String name
    File config_file
}

Array[File?] misc_file_list = flatten(misc_file_arrays)

command <<<
    mkdir -p data
    for file in ~{sep(' ', misc_file_list)}; do
        cp $file data
        echo $file >> data_files.txt
    done
    for file in ~{sep(' ', misc_files)}; do
        cp $file data
        echo $file >> data_files.txt
    done
    echo "Contents of data directory:"
    ls -R data
    multiqc  --file-list data_files.txt -o . -n ~{name} --config ~{config_file} --interactive --verbose
>>>

output {
    File report = "~{name}.html"
                         }

runtime {
    docker: "multiqc/multiqc:dev"
}
}
