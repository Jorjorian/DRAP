version 1.0

task fastqc {
    input {
        File fq
        String tag
        String read
        Int cpu = 3
    }
    command <<<
        touch "~{tag}_~{read}_fastqc.html"
        touch "~{tag}_~{read}_fastqc.zip"
        cp ~{fq} ~{tag}_~{read}.fastq.gz
        fastqc --threads ~{cpu} --nogroup --noextract --outdir . ~{tag}_~{read}.fastq.gz

    >>>
    output {
        File html = "~{tag}_~{read}_fastqc.html"
        File zip = "~{tag}_~{read}_fastqc.zip"
    }
    runtime {
        docker: "staphb/fastqc"
        memory: "2 GB"
        cpu: cpu
        returnCodes: "*"



    }
    parameter_meta {
        fq: {description: "input fastq file", category: "required"}
    }
}