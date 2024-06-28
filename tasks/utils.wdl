version 1.0


task fastaIndex {
    input {
        File fasta
        String docker = ""
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

