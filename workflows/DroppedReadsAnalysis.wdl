version development
import "../tasks/fastqc.wdl" as fastqc
import "../tasks/kraken.wdl" as kraken
import "../tasks/utils.wdl" as utils
import "../tasks/consolidate_samples.wdl" as consolidate_samples
import "../tasks/flowcell_metrics.wdl" as flowcell_metrics
import "../tasks/dimer_metrics.wdl" as dimer_metrics
import "../tasks/multiQC.wdl" as multiQC
import "../tasks/bwa.wdl" as bwa


workflow DRAP {

  input {
      Array[Pair[String, Array[Pair[File, File]]]] dropped_read_groups
      File kraken_tar
      File reference_fasta
      File multiqc_config
      File? primer_list
      File adapter_list
  }

  call consolidate_samples.consolidate_samples as consolidate_samples {
      input:
          dropped_read_groups = dropped_read_groups
  }
  call utils.fastaIndex as ref {
        input:
            fasta = reference_fasta,
    }

  ## scatter bcl2fastq processing
  scatter(group in consolidate_samples.output_fastq_pairs){
#      call contamination_blast.blast as blast {
#          input:
#              read_group = group
#          }
      call kraken.runKraken as contamination_kraken {
          input:
              tag = group.left,
              fq_1 = group.right.left,
              fq_2 = group.right.right,
              Kraken2_tar = kraken_tar
          }
    call fastqc.fastqc as r1_fastqc {
        input:
            fq = group.right.left,
            tag = group.left,
            read = '1'
    }
    call fastqc.fastqc as r2_fastqc {
        input:
            fq = group.right.right,
            tag = group.left,
            read = '2'
    }
    call flowcell_metrics.flowcell_metrics as flowcell_metrics {
        input:
            fq_1 = group.right.left,
            fq_2 = group.right.right,
            tag = group.left
    }
    call dimer_metrics.dimer_metrics as dimer_metrics {
        input:
            fq_1 = group.right.left,
            fq_2 = group.right.right,
            tag = group.left,
            adapter_list = adapter_list
    }

    call bwa.loose_mapping_metrics as loose_mapping {
        input:
            fq1 = group.right.left,
            fq2 = group.right.right,
            tag = group.left,
            ref = ref.bwa_index
    }
}
    call dimer_metrics.combine_dimer_metrics as combine_dimer_metrics {
        input:
            dimer_metrics_files = dimer_metrics.dimer_metrics_json
    }
    call flowcell_metrics.consolidate_flowcell_metrics as consolidate_flowcell_metrics {
        input:
            flowcell_metrics_files = flowcell_metrics.metrics
    }
    call multiQC.multiQC as multiQC {
        input:
            misc_file_arrays = [loose_mapping.stats, loose_mapping.bam_json, loose_mapping.depth,
                                dimer_metrics.dimer_metrics_json, flowcell_metrics.metrics, flowcell_metrics.overall_png,
                                flowcell_metrics.clusters, contamination_kraken.kraken_report,
                                contamination_kraken.kraken_file, r1_fastqc.html, r1_fastqc.zip, r2_fastqc.html,
                                r2_fastqc.zip],
            misc_files = [consolidate_flowcell_metrics.multiqc_flowcell_metrics_json,
                          combine_dimer_metrics.multiqc_dimer_metrics_json],
            name = 'multiqc',
            config_file = multiqc_config
    }
    output {
        Array[File] stats = loose_mapping.stats
        File multiqc_html = multiQC.report
      }
}