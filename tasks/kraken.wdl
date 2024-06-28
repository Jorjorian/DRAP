version development

task runKraken {
  input {
    File fq_1
    File fq_2
    String tag
    File Kraken2_tar
    String kraken2_db = "kraken2_db"
    String docker
    Int memory = 8
    Int cpu = 4
  }
  command <<<
    ls ~{kraken2_db}
    mkdir Kraken_output
    touch Kraken_output/krona_input.txt
    touch Kraken_output/~{tag}.kraken_output
    echo ~{Kraken2_tar}
    if [ ! -d "~{kraken2_db}" ]; then
        echo "Creating kraken2_db directory"
        mkdir -p ~{kraken2_db}
    fi
    tar -xvf ~{Kraken2_tar} -C ~{kraken2_db}
    kraken2 --db ~{kraken2_db} --threads ~{cpu} --output Kraken_output/~{tag}.kraken_output --report Kraken_output/~{tag}.kreport2 ~{fq_1} ~{fq_2}
    cat Kraken_output/~{tag}.kraken_output | cut -f 2,3 > Kraken_output/krona_input.txt
    python3 <<CODE
    import json
    with open("Kraken_output/~{tag}.kreport2") as infile:
        H = {"sample_index": "~{tag}"}
        for line in infile:
            line = line.strip().split()
            perc, count, clade, taxon = line[0], line[1], line[3], line[-1]
            if clade in ["C", "O", "F", "G", "S"]:
                try:
                    assert H.get(clade, False)
                except AssertionError:
                    # Add a new entry to dict
                    H[clade] = [taxon, perc, count]
                else:
                    if float(H[clade][1]) < float(perc) and int(H[clade][2]) < int(
                        count
                    ):
                        # If new line in kreport file indicates higher abundance percent
                        # and count, then overwrite clade information with this new line
                        # information
                        H[clade] = [taxon, perc, count]
    with open("~{tag}_kraken.json", "w") as OUT:
        json.dump(H, OUT)
    CODE

  >>>
  output {
    File kraken_report = "Kraken_output/~{tag}.kreport2"
    File kraken_file = "Kraken_output/~{tag}.kraken_output"
    File krona_in = "Kraken_output/krona_input.txt"
    File kraken_json = "~{tag}_kraken.json"
  }
  runtime {
    docker: "staphb/kraken2:latest"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk 10 SSD"
    preemptible: 0
  }
}

task runKrona {
  input {
    File krona_file
    String tag
#    String strain_id
    Directory taxonomy_db
    String docker = "quay.io/biocontainers/krona:2.8--pl526_1"
    Int memory = 8
    Int cpu = 4
    }
  command <<<
    touch ~{tag}.krona.html
    ktImportTaxonomy ~{krona_file} -tax ~{taxonomy_db} -o ~{tag}.krona.html

  >>>
  output {
    File html_report = "~{tag}.krona.html"
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk 10 SSD"
    preemptible: 0
    }
}
