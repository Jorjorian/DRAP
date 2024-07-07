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
      Directory kraken_db
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
  scatter(group in consolidate_samples.groups){
      call contamination_blast.blast as blast {
          input:
              read_group = group
          }
      call kraken.kraken as contamination_kraken {
          input:
              sample_id = group[0],
              fq_1 = group[1][0],
              fq_2 = group[1][1],
              kraken_db = kraken_db
          }
    call fastqc.fastqc as r1_fastqc {
        input:
            fq_1 = group[1][0],
            tag = group[0]

    }
    call fastqc.fastqc as r2_fastqc {
        input:
            fq_1 = group[1][1],
            tag = group[0]
    }
    call flowcell_metrics.flowcell_metrics as flowcell_metrics {
        input:
            fq_1 = group[1][0],
            fq_2 = group[1][1],
            tag = group[0]
    }
    call dimer_metrics.dimer_metrics as dimer_metrics {
        input:
            fq_1 = group[1][0],
            fq_2 = group[1][1],
            tag = group[0],
            primer_list = primer_list,
            adapter_list = adapter_list

    }
    call bwa.loose_mapping_metrics as loose_mapping {
        input:
            fq_1 = group[1][0],
            fq_2 = group[1][1],
            tag = group[0],
            ref = ref.bwa_index
    }


}
    call multiQC.multiQC as multiQC {
        input:
            misc_files = [],
            misc_file_arrays = [],
            name = 'multiqc',
            config_file = multiqc_config
    }
    output {
#        File? metrics = process_samples.metrics
        File final_report = multiQC.report
      }
}