#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

date=new Date().format( 'yyMMdd' )
user="$USER"

// Preset parameters:

params.server                           ="lnx01"
params.outdir                           ='RNA_results'
params.rundir                           ="${launchDir.baseName}" 
params.genome                           ="hg38" //default assembly is hg38, unless --genome is set!
params.rsem_strand                      ="reverse" // TO DO: INCLUDE STRANDEDNESS SELECTION!!!
outdir_full_path                        = "${launchDir}/${params.outdir}/"
params.gatk                             ="new"

// Unset parameters:
params.help                             =false
params.data                             =null
params.fastq                            =null
params.primary                          =null
// skip or include individual analysis - unset parameters:
params.dupradar                         =null   // not run by default
params.qualimap                         =null   // not run by default
params.skipVariants                     =null
params.skipTrim                         =null
params.skipQC                           =null
// sub-workflow deselect:
params.skipAltSplicing                  =null
params.skipFusion                       =null
params.skipExpression                   =null
// Assembly-independent variables:

//inhouse_genelist="/data/shared/genomes/databases/genelists/tumortarget/tumortarget.inhouse.v2.127genes.txt"

inhouse_genelist="/data/shared/genomes/databases/genelists/tumortarget/240123.inhouse.MOMA.Fusion.241genes.for.grep.txt"

multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
switch (params.gatk) {

    case 'danak':
    gatk_image="gatk419.sif";
    break;
    case 'new':
    gatk_image="gatk4400.sif";
    break;
    default:
    gatk_image="gatk4500.sif";
    break;
}


switch (params.server) {
    case 'lnx02':
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/,/fast/:/fast/,/lnx01_data3/:/lnx01_data3/";
        simgpath="/data/shared/programmer/simg";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/wgs_splitinterval_BWI_subdivision3/*.interval_list";
        tmpDIR="/fast/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
        tank_storage="/home/mmaj/tank.kga/data/data.storage.archive/";
        data_archive="/lnx01_data2/shared/dataArchive/";
        genomes_dir="/fast/shared/genomes"
        //modules_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules/";
    break;  
    case 'lnx01':
       // syspath="/data/shared";
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/";
        simgpath="/data/shared/programmer/simg";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/wgs_splitinterval_BWI_subdivision3_GATK_callable/*.interval_list";
        tmpDIR="/data/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/gatk4261.sif gatk";
        multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
        tank_storage="/home/mmaj/tank.kga2/data/data.storage.archive";
        data_archive="/lnx01_data2/shared/dataArchive/";
        modules_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules";
    break;
    case 'kga01':
        simgpath="/data/shared/programmer/simg";
      //  syspath="/data/shared";
        s_bind="/data/:/data/";
        tmpDIR="/data/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/gatk4261.sif gatk";
        tank_storage="/home/mmaj/tank.kga/data/data.storage.archive";
        data_archive="/data/shared/dataArchive/";
        modules_dir="/home/mmaj/LNX01_mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules";
    break;
}
/*
*/

def helpMessage() {
    log.info"""

    Generel info:
    Requires a samplesheet containing 5 columns in specific order (tab separated), without headerline:
    1) caseID, 2) NPN normal WES, 3) NPN tumor WES, 4) NPN tumor RNA, 5) PCGR tumor value

    Example samplesheet:

    johnDoe 112217976652	111184925465	111184925473    23

    The script will automatically look for fastq in subfolders at KG Vejle data archive (fastq storage incl novaRuns).

    The user can point to a specific folder containing raw data (FastQ) using the --fastq option 
    This is only needed if raw data only exists outside the data archive (e.g. if data are in personal folders or at other KG Vejle analysis servers).


    Main options:
      --help                print this help message

      --genome              hg19 or hg38
                                Default: hg38

      --outdir              Select which folder to write output to.
                                Default: RNA_results

      --samplesheet         path to case samplesheet. Can contain multiple patients/cases (each case in a separate line). 

      --server              Select which server the analysis is performed on (lnx01 or lnx02)
                                Default: lnx01

      --fastq               Path to dir with fastq files
                                Default: data storage dirs at lnx01 server or kga01 server

    QC options:
      --skipQC              Skip all QC steps except MultiQC
      --qualimap            Run qualimap BAMQC
                                Default: Do not run BAMQC (timeconsuming!)
      --dupradar            Run Dupradar (not available this version of the script)
                                Default: Do not run Dupradar (timeconsuming!)

    Analysis selection options:
      --skipVariants        Do not call variants (HaplotypeCaller) on RNA seq data
      --skipFusion          Do not call Fusions
      --skipAltSplicing     Do not call alternative splicevariants (e.g. exon skipping)
      --skipExpression      Do not analyse geneexpression (i.e. skip RSEM, HTseq and FeatureCounts)
      
    """.stripIndent()
}
if (params.help) exit 0, helpMessage()



// config files:


///////////////////////// SAMPLESHEET CHANNELS /////////////////////

// Samplesheet cols (fixed order)
// 0: CaseID, 1: WES.blood, 2: WES.tumor, 3: RNA tumor, 4: pcgr_tumor_code

////////////////////////////////////////////////////////////////////////////////

channel.fromPath(params.samplesheet)
    .splitCsv(sep:'\t')
    .map { row -> tuple(row[3], row[0])}
    .set {case_rna_ss}                                     //sampleID(NPN), caseID (from samplesheet)

if (params.fastq) {
    params.reads="${params.fastq}/**{.,-}{RV1}{.,-}*R{1,2}*{fq.gz,fastq,fastq.gz}"
    
    channel.fromFilePairs(params.reads, checkIfExists: true)
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .map { it -> [it[0], file(it[1][0]),file(it[1][1])] }
        .set { readpairs_ch }                           //sampleID(NPN), r1, r2
}

if (!params.fastq) {
    params.reads="${data_archive}/{lnx01,kga01_novaRuns,tank_kga_external_archive}/**/*{.,-}{RV1}{.,-}*R{1,2}*{fq,fastq}.gz"

    channel.fromFilePairs(params.reads, checkIfExists: true)
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .map { it -> [it[0], file(it[1][0]),file(it[1][1])] }
        .set { readpairs_ch }                           //sampleID(NPN), r1, r2
}

// TO DO FROM HERE TO DSL2
case_rna_ss
    .join(readpairs_ch)
    .map {tuple(it[1],it[0]+"_RNA",it[2],it[3])}
    .set {case_sample_reads_ch}                    //caseID, sampleID, r1, r2
  
case_sample_reads_ch
    .map {tuple(it[0],it[2],it[3])}
    .set { trinity_splicing_readinput_ch }          // caseID, r1, r2

/*    
process jaffa_kga01 {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${params.outdir}/genefusions/JAFFA/", mode: 'copy'

    conda '/data/shared/programmer/miniconda3/envs/bpipe'
    input:
    tuple val(caseID), val(sampleID), file(r1), file(r2)

    output:
    path("${sampleID}* /")
    script:
    """
    bpipe run \
    /data/shared/programmer/JAFFA-version-2.3/JAFFA_direct.groovy \
    ${r1} ${r2}
    """
}
*/

include { // SUB workflows:
         SUB_RNA_PREPROCESS;                    // trim, align with STAR 2-pass basic (twice: for Arriba and for star-fusion)
         SUB_RNA_QC;                            // various QC tools
         SUB_RNA_EXPRESSION;                    // expresssion tools (RSEM, FeatureCounts, HTSeq)    
         SUB_RNA_ALT_SPLICING;
         SUB_RNA_FUSION } from "./modules/tumorBoard.rna.modules.v2.nf" 


workflow {
    SUB_RNA_PREPROCESS(case_sample_reads_ch)
   // jaffa_kga01(case_sample_reads_ch)

    if (!params.skipQC) {
        SUB_RNA_QC(SUB_RNA_PREPROCESS.out.star_bam)
    }
    if (!params.skipExpression) {
        SUB_RNA_EXPRESSION(SUB_RNA_PREPROCESS.out.star_bam,
                       SUB_RNA_PREPROCESS.out.star_rsem_bam)
    }
    
    if (!params.skipFusion) {
    SUB_RNA_FUSION(case_sample_reads_ch,
                   SUB_RNA_PREPROCESS.out.star_arriba_bam,
                   SUB_RNA_PREPROCESS.out.star_chimeric_junctions)
    }
    // set trinity splicing input channel: tuple val(caseID), val(sampleID), path(bam),path(bai), path(junction), path(r1),path(r2),path(sj_tab):
    
    if (!params.skipAltSplicing) {
    SUB_RNA_PREPROCESS.out.star_bam
    .join(SUB_RNA_PREPROCESS.out.star_junctions)
    .join(trinity_splicing_readinput_ch)
    .join(SUB_RNA_PREPROCESS.out.star_sjtab)
    .set{trinity_splicing_input}
 
    SUB_RNA_ALT_SPLICING(trinity_splicing_input)
    }    
        //SUB_RNA_PREPROCESS.out.star_bam.join(SUB_RNA_PREPROCESS.out.star_junctions).join(trinity_splicing_readinput_ch).join(SUB_RNA_PREPROCESS.out.star_sjtab))
}







