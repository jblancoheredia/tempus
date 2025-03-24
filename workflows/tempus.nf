/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT PLUGINS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                                                              } from 'plugin/nf-schema'
include { samplesheetToList                                                             } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                                                                       } from '../modules/nf-core/multiqc/main'
include { BWAMEM2_MEM                       as BWAMEM2_MEM_NORMAL                       } from '../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_MEM                       as BWAMEM2_MEM_TUMOUR                       } from '../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_SORT                     as SAMTOOLS_SORT_NORMAL                     } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SAMTOOLS_SORT_TUMOUR                     } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX                    as SAMTOOLS_INDEX_NORMAL                    } from '../modules/nf-core/samtools/index/main' 
include { SAMTOOLS_INDEX                    as SAMTOOLS_INDEX_TUMOUR                    } from '../modules/nf-core/samtools/index/main' 
include { BEDTOOLS_BAMTOFASTQ               as BEDTOOLS_BAMTOFASTQ_NORMAL               } from '../modules/local/bedtools/bamtofastq/main'
include { BEDTOOLS_BAMTOFASTQ               as BEDTOOLS_BAMTOFASTQ_TUMOUR               } from '../modules/local/bedtools/bamtofastq/main'
include { TELOMEREHUNTER_FULL                                                           } from '../modules/local/telomerehunter/full/main'
include { TELFUSDETECTOR_FULL                                                           } from '../modules/local/telfusdetector/full/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                 IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMultiqc                                                          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                                        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                                        } from '../subworkflows/local/utils_nfcore_tempus_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                CREATE CUSTOM CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_fai                                          = Channel.fromPath(params.fai).map                              { it -> [[id:it.Name], it] }.collect()
ch_bwa2                                         = Channel.fromPath(params.bwa2).map                             { it -> [[id:it.Name], it] }.collect()
ch_fasta                                        = Channel.fromPath(params.fasta).map                            { it -> [[id:it.Name], it] }.collect()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                 RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TEMPUS {

    take:
    ch_samplesheet
    
    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samplesheet
        .map { normal_meta, bam_n, bai_n, tmuour_meta, bam_t, bai_t -> tuple(normal_meta, bam_n) }
        .set { ch_normal_raw_bam }

    ch_samplesheet
        .map { normal_meta, bam_n, bai_n, tmuour_meta, bam_t, bai_t -> tuple(tmuour_meta, bam_t) }
        .set { ch_tumour_raw_bam }

    //
    // MODULE: Run SAMtools Sorting by Name
    //
    SAMTOOLS_SORT_NORMAL(ch_normal_raw_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_NORMAL.out.versions)
    ch_normal_name_sorted_bam = SAMTOOLS_SORT_NORMAL.out.bam
    
    //
    // MODULE: Run SAMtools Sorting by Name
    //
    SAMTOOLS_SORT_TUMOUR(ch_tumour_raw_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_TUMOUR.out.versions)
    ch_tumour_name_sorted_bam = SAMTOOLS_SORT_TUMOUR.out.bam

    //
    // MODULE: Run BEDtools BAMtoFASTQ
    //
    BEDTOOLS_BAMTOFASTQ_NORMAL(ch_normal_name_sorted_bam)
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOFASTQ_NORMAL.out.versions)
    ch_normal_fastq = BEDTOOLS_BAMTOFASTQ_NORMAL.out.fastq

    //
    // MODULE: Run BEDtools BAMtoFASTQ
    //
    BEDTOOLS_BAMTOFASTQ_TUMOUR(ch_tumour_name_sorted_bam)
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOFASTQ_TUMOUR.out.versions)
    ch_tumour_fastq = BEDTOOLS_BAMTOFASTQ_TUMOUR.out.fastq

    //
    // MODULE: Run BWAMEM2
    //
    BWAMEM2_MEM_NORMAL(ch_normal_fastq, ch_bwa2, ch_fasta, true)
    ch_versions = ch_versions.mix(BWAMEM2_MEM_NORMAL.out.versions)
    ch_normal_bam = BWAMEM2_MEM_NORMAL.out.bam

    //
    // MODULE: Run BWAMEM2
    //
    BWAMEM2_MEM_TUMOUR(ch_tumour_fastq, ch_bwa2, ch_fasta, true)
    ch_versions = ch_versions.mix(BWAMEM2_MEM_TUMOUR.out.versions)
    ch_tumour_bam = BWAMEM2_MEM_TUMOUR.out.bam

    //
    // MODULE: Run TelomereHunter
    //
    TELOMEREHUNTER_FULL(ch_tumour_bam, ch_normal_bam, params.banding)
    ch_versions = ch_versions.mix(TELOMEREHUNTER_FULL.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'tempus_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
