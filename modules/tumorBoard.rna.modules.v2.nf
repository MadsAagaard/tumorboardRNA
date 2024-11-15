#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


date=new Date().format( 'yyMMdd' )
user="$USER"
runID="${date}.${user}"


//inhouse_genelist="/data/shared/genomes/databases/genelists/tumortarget/tumortarget.inhouse.v2.127genes_for_grep.txt"

inhouse_genelist="/data/shared/genomes/databases/genelists/tumortarget/240123.inhouse.MOMA.Fusion.241genes.for.grep.txt"

inhouse_splicing_genelist="/data/shared/genomes/databases/genelists/tumortarget/tumortarget.inhouse.v2.127genes_for_grep.txt"

//////////////////////////// SWITCHES ///////////////////////////////// 

switch (params.gatk) {

    case 'danak':
    gatk_image="gatk419.sif";
    break;
    case 'new':
    gatk_image="gatk4400.sif";
    break;
    default:
    gatk_image="gatk4400.sif";
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
      //  syspath="/data/shared";
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/";
        simgpath="/data/shared/programmer/simg";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/wgs_splitinterval_BWI_subdivision3_GATK_callable/*.interval_list";
        tmpDIR="/data/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/gatk4261.sif gatk";
        multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
        tank_storage="/home/mmaj/tank.kga2/data/data.storage.archive";
        data_archive="/lnx01_data2/shared/dataArchive/";
        genomes_dir="/data/shared/genomes"
    break;
    case 'kga01':
        simgpath="/data/shared/programmer/simg";
  //      syspath="/data/shared";
        s_bind="/data/:/data/";
        tmpDIR="/data/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/gatk4261.sif gatk";
        tank_storage="/home/mmaj/tank.kga/data/data.storage.archive";
        data_archive="/data/shared/dataArchive/";
        genomes_dir="/data/shared/genomes"
    break;
}

switch (params.genome) {
    case 'hg19':
        // Genome assembly files:
        genome_fasta = "${genomes_dir}/hg19/human_g1k_v37.fasta"
        genome_fasta_fai = "${genomes_dir}/hg19/human_g1k_v37.fasta.fai"
        genome_fasta_dict = "${genomes_dir}/hg19/human_g1k_v37.dict"

        // Gene and transcript annotation files:
        gene_bed12 = "${genomes_dir}/hg19/gene.annotations/hg19.refGene.201228.BED12.bed"
        transcript_fasta="${genomes_dir}/GRCh37/gene.annotations/gencode.v19.pc_transcripts.fa"
        //gencode_gtf = "${genomes_dir}/GRCh37/gene.annotations/gencode.v19.annotation.gtf"
        gencode_gtf = "${genomes_dir}/hg19/gene.annotations/hg19.refGene.201228.gtf" // make sure to change this after testing!
         //Program  files:
        msisensor_list="${genomes_dir}/hg19/human_g1k_v37.microsatellites.list"
        genome_lib_starfusion= "${genomes_dir}/hg19/CTAT/GRCh37.Apr032020.PNP/ctat_genome_lib_build_dir/"
        arriba_blacklist = "/data/shared/programmer/arriba_v2.0.0/database/blacklist_hg19_hs37d5_GRCh37_v2.0.0.tsv.gz"
        arriba_cytoband= "/data/shared/genomes/hg19/arriba/cytobands_hg19_hs37d5_GRCh37_2018-02-23.tsv"
        arriba_protein_gff = "/data/shared/genomes/hg19/arriba/protein_domains_hg19_hs37d5_GRCh37_2019-07-05.gff3"

        // Program indexes
        index_rsem = "/data/shared/genomes/GRCh37/rsem/rsem_gencode"
        //index_star = "/data/shared/genomes/GRCh37/star_align_grch37"
        index_star = "/data/shared/genomes/hg19/star_align_hg19"

        //regions:
        qualimap_ROI="/data/shared/genomes/hg19/interval.files/200108.NCBIrefseq.codingexons.nocontig.20bp.merged.sorted.6col.bed"
        ROI="/data/shared/genomes/hg19/interval.files/WES/IDT.exomes.EV7/EV7.ROI.bed"

        break;
    case 'hg38':
        // Genome assembly files:
        genome_fasta = "${genomes_dir}/hg38/GRCh38.primary.fa"
        genome_fasta_fai = "${genomes_dir}/hg38/GRCh38.primary.fa.fai"
        genome_fasta_dict = "${genomes_dir}/hg38/GRCh38.primary.dict"

        // Gene and transcript annotation files:
        gene_bed12="${genomes_dir}/hg38/gene.annotations/gencode.v36.BED12.bed"
        transcript_fasta="${genomes_dir}/hg38/gene.annotations/gencode.v36.transcripts.fa"
        gencode_gtf = "${genomes_dir}/hg38/gene.annotations/gencode.v36.annotation.gtf"
        gencode_gff3 = "${genomes_dir}/hg38/gene.annotations/gencode.v36.annotation.gff3"
        gencode_gtf_collapsed="${genomes_dir}/hg38/gene.annotations/gencode.v36.annotation.collapsed.gtf"
        ensemble102_gtf="${genomes_dir}/hg38/gene.annotations/ensembl102/Homo_sapiens.GRCh38.102.gtf"
        ensemble102_transcript_fasta="${genomes_dir}/hg38/gene.annotations/ensembl102/Homo_sapiens.GRCh38.102.cdna.all.fa.gz"
        //Program  files:
        msisensor_list="/data/shared/genomes/hg38/program_DBs/msisensor/hg38_msisensor_scan.txt"
        genome_lib_starfusion="/data/shared/genomes/hg38/program_DBs/CTAT/GRCh38_v33_Apr062020.PNP/ctat_genome_lib_build_dir/"
        arriba_cytoband="/data/shared/genomes/hg38/program_DBs/arriba/cytobands_hg38_GRCh38_v2.1.0.tsv"
        arriba_blacklist = "/data/shared/genomes/hg38/program_DBs/arriba/blacklist_hg38_GRCh38_v2.1.0.tsv.gz"
        arriba_protein_gff="/data/shared/genomes/hg38/program_DBs/arriba/protein_domains_hg38_GRCh38_v2.1.0.gff3"
        arriba_known_fusions="/data/shared/genomes/hg38/program_DBs/arriba/known_fusions_hg38_GRCh38_v2.1.0.tsv.gz"

        //fusioncallers
        fusioncatcher_db="/data/shared/genomes/hg38/program_DBs/fusioncatcher/human_v102/"
        fusionreport_db="/data/shared/genomes/hg38/program_DBs/fusion_report_db/"
        // Program indexes:
        index_rsem = "/data/shared/genomes/hg38/rsem/rsem_hg38_gencode36"
        index_star = "${genomes_dir}/hg38/STAR/"
        kallisto_index="/data/shared/genomes/hg38/program_DBs/kallisto/Homo_sapiens.GRCh38.102.cdna.all_kallisto_K31.idx"
        //regions:
        qualimap_ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.SM.6col.bed"
        ROI="/data/shared/genomes/genomes/hg38/interval.files/exome.ROIs/211130.hg38.refseq.gencode.fullexons.50bp.SM.bed"
        gencode36_coding_exons="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.SM.bed"
        break;
}








///////////////////////// PREPROCESS MODULES ////////////////////////////////////

process inputFiles_symlinks_fq{
    errorStrategy 'ignore'
    publishDir "${caseID}/${params.outdir}/input_symlinks/", mode: 'link', pattern:'*.{fastq,fq}.gz'
    input:
    tuple val(caseID), val(sampleID), path(r1),path(r2)// from read_input2
    
    output:
    tuple path(r1),path(r2)
    script:
    """
    """
}


process fastp_TRIM {
    publishDir "${caseID}/${params.outdir}/QC/", mode: 'copy', pattern: '*.{html,json}'
    cpus 10
    tag "$sampleID"

    input:
    tuple val(caseID), val(sampleID), file(r1), file(r2)// from NN1

    output:
    path("*.{html,json}"),                                      emit: fastp_results
    tuple val(caseID), val(sampleID), path("*.fastp.fq.gz"),    emit: fastp_out_ch

    script:
    """
    singularity run -B ${s_bind} \
    ${simgpath}/fastp.sif \
    -i ${r1} -I ${r2} \
    -o ${r1.baseName}.fastp.fq.gz -O ${r2.baseName}.fastp.fq.gz \
    --json ${sampleID}.fastp.json \
    --html ${sampleID}.fastp.html \
    -w ${task.cpus}
    """
}

process align_STAR {
    tag "$sampleID"
    cpus 40
    maxForks 3
    publishDir "${caseID}/${params.outdir}/BAM/", mode: 'copy', pattern: '*.STAR.Aligned.*'
    publishDir "${caseID}/${params.outdir}/QC/", mode: 'copy', pattern: '*.Log.*'
    publishDir "${caseID}/${params.outdir}/BAM/", mode: 'copy', pattern: '*.forArriba.Aligned.*'
 
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/star2711b'

    //read_pairs_star_ch
    input:
    tuple val(caseID), val(sampleID), path(reads)// from fastp_out_ch
    
    output:
    tuple val(caseID), val(sampleID), path("${caseID}.${sampleID}.STAR.Aligned.sortedByCoord.*.bam"),path("${caseID}.${sampleID}.STAR.Aligned.sortedByCoord.*.bai"), emit: star_out_bam_ch // into (split_cigar_input_bam, star_bam_ch2,star_bam_ch3,star_bam_ch4, preseq_input_bam,rseqc_input_bam,qualimapRNA_input_bam,qualimapBAMQC_input_bam,dupradar_input_bam, rnaseqc_input_bam, featurecounts_input_bam,htseq_count_input_bam, trinityfusion_input_bam)

    tuple val(caseID), val(sampleID), path("${caseID}.${sampleID}.STAR.Aligned.toTranscriptome.*.bam"), emit:rsem_input_bam    // into rsem_input_bam
    
    tuple val(caseID), path("*.Chimeric.out.junction"),emit: star_junction_out//  into (trinityfusion_input_junction,trinity_splicing_input_junction)

    tuple val(caseID), val(sampleID), path("*.Chimeric.out.junction"),emit: chimeric_junctions_out  //into chimeric_junctions_out

    tuple val(caseID), path("*.STAR.SJ.out.tab"),emit: star_sjtab_out // into star_sjtab_out
    
    tuple val(caseID),val(sampleID), path("${caseID}.${sampleID}.forArriba.*.bam"),path("${caseID}.${sampleID}.forArriba.*.bai"), emit: arriba_input_bam

    path("*.STAR.Log.*")
    path("*.forArriba.Log.out"), emit: trinity_collect // into (trinity_collect_ch1,trinity_collect_ch2)
    script:
    """
    STAR --runThreadN ${task.cpus} \
        --genomeDir ${index_star} \
        --sjdbGTFfile ${gencode_gtf} \
        --readFilesIn ${reads} \
        --outFileNamePrefix ${caseID}.${sampleID}.STAR. \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif --outSAMunmapped Within \
        --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --outSAMattrRGline ID:${date}.${user} LB:library PL:illumina PU:illumina SM:${caseID}.${sampleID} \
        --chimSegmentMin 12 \
        --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 \
        --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 \
        --chimJunctionOverhangMin 12 --chimOutJunctionFormat 1 \
        --peOverlapNbasesMin 12 --peOverlapMMp 0.1 \
        --outReadsUnmapped None \
        --alignInsertionFlush Right \
        --alignSplicedMateMapLminOverLmate 0 \
        --alignSplicedMateMapLmin 30 \
        --quantMode TranscriptomeSAM GeneCounts \
        --readFilesCommand zcat
 
    STAR --runThreadN ${task.cpus} \
        --genomeDir ${index_star} \
        --sjdbGTFfile ${gencode_gtf} \
        --readFilesIn ${reads} \
        --outFileNamePrefix ${caseID}.${sampleID}.forArriba. \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif --outSAMunmapped Within \
        --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --outSAMattrRGline ID:${date}.${user} LB:library PL:illumina PU:illumina SM:${caseID}.${sampleID} \
        --chimSegmentMin 12 \
        --chimOutType WithinBAM \
        --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 \
        --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 \
        --chimJunctionOverhangMin 12 --chimOutJunctionFormat 1 \
        --peOverlapNbasesMin 12 --peOverlapMMp 0.1 \
        --outReadsUnmapped None \
        --readFilesCommand zcat
    
    bamtools index -in ${caseID}.${sampleID}.forArriba.Aligned.sortedByCoord.out.bam
    bamtools index -in ${caseID}.${sampleID}.STAR.Aligned.sortedByCoord.out.bam
    """
}


process GATK_splitNCIGAR {
    tag "$sampleID"

    input:
    tuple val(caseID),val(sampleID), path(bam), path(bai) from split_cigar_input_bam
    output:
    tuple val(caseID), val(sampleID), path("${caseID}.${sampleID}.CIGARsplit.bam"), path("${caseID}.${sampleID}.CIGARsplit.bai") into HC_input_bam
    when:
    !params.skipVariants
    script:
    """
    ${gatk_exec} SplitNCigarReads \
    -R ${genome_fasta} \
    -I ${bam} \
    -O ${caseID}.${sampleID}.CIGARsplit.bam 
    """
}


process GATK_haplotypecaller{
    cpus 10
    tag "$sampleID"
    publishDir "${caseID}/${params.outdir}/variantcalls/", mode: 'copy', pattern: "*.HC.*"
    publishDir "${caseID}/${params.outdir}/gvcf/", mode: 'copy', pattern: "*.g.*"
    input:
    tuple val(caseID), val(sampleID), path(bam), path(bai)// from HC_input_bam
    
    output:
    tuple val(caseID), val(sampleID), path("${caseID}.${sampleID}.g.vcf")// into sample_gvcf_tuple
    tuple val(caseID), val(sampleID),  path("*.${sampleID}.HC.vcf"), path("*.${sampleID}.HC.vcf.idx")// into id_HCvcf_idx

    path("${caseID}.${sampleID}.g.vcf")// into sample_gvcf_list


    when:
    !params.skipVariants
    script:
    """
    ${gatk_exec} HaplotypeCaller \
    -R ${genome_fasta} \
    -I ${bam} \
    -O ${caseID}.${sampleID}.g.vcf \
    -ERC GVCF \
    -L ${ROI} 
            
    ${gatk_exec} GenotypeGVCFs \
    -R ${genome_fasta} \
    -V ${caseID}.${sampleID}.g.vcf \
    -O ${caseID}.${sampleID}.HC.vcf \
    -G StandardAnnotation \
    -G AS_StandardAnnotation
    """
}

process rseqc {
    tag "$sampleID"
    publishDir "${caseID}/${params.outdir}/QC/rseqc/", mode: 'copy'


    conda '/lnx01_data3/shared/programmer/miniconda3/envs/rseqc'

    input:
    tuple val(caseID),val(sampleID), path(bam), path(bai)// from rseqc_input_bam

    output:
    //path("*.{txt,pdf,r,xls}")
    tuple val(caseID),val(sampleID), path("*.{txt,pdf,r,xls}")

    when:
    !params.skipQC

    script:
    """
    infer_experiment.py -i ${bam} -r ${gene_bed12} > ${caseID}.${sampleID}.rseqc.infer_experiment.txt
    bam_stat.py -i ${bam} > ${caseID}.${sampleID}.rseqc.bamstat.txt
    junction_annotation.py -i ${bam} -r ${gene_bed12} -o ${caseID}.${sampleID}.rseqc
    junction_saturation.py -i ${bam} -r ${gene_bed12} -o ${caseID}.${sampleID}.rseqc
    inner_distance.py -i ${bam} -r ${gene_bed12} -o ${caseID}.${sampleID}.rseqc
    read_distribution.py -i ${bam} -r ${gene_bed12} > ${caseID}.${sampleID}.rseqc.read.dist.txt
    read_duplication.py -i ${bam} -o ${caseID}.${sampleID}.rseqc.read.dup
    """
}

/*
process qualimapRNAseq {
    tag "$sampleID"
    publishDir "${caseID}/${params.outdir}/QC/", mode: 'copy'

    input:
    tuple val(caseID),val(sampleID), path(bam), path(bai)// from qualimapRNA_input_bam
    
    output:
    path ("${caseID}.${sampleID}.qualimapRNA/")

    when:
    !params.skipQC

    script:
    """
    qualimap --java-mem-size=10G rnaseq -outdir ${caseID}.${sampleID}.qualimapRNA -bam ${bam} -gtf ${gencode_gtf} -pe 
    """
}
*/

process qualimapRNAseq {
    tag "$sampleID"
    publishDir "${caseID}/${params.outdir}/QC/", mode: 'copy'

    input:
    tuple val(caseID),val(sampleID), path(bam), path(bai)// from qualimapRNA_input_bam
    
    output:
    path ("${caseID}.${sampleID}.qualimapRNA/")

    when:
    !params.skipQC

    script:
    """
    unset DISPLAY
    singularity run -B ${s_bind} ${simgpath}/qualimap.sif qualimap --java-mem-size=10G \
    rnaseq \
    -outdir ${caseID}.${sampleID}.qualimapRNA \
    -bam ${bam} -gtf ${gencode_gtf} -pe 
    """
}


process qualimapBAMQC {
    tag "$sampleID"
    publishDir "${caseID}/${params.outdir}/QC/", mode: 'copy'
    cpus 10

    input:
    tuple val(caseID),val(sampleID), path(bam), path(bai)// from qualimapBAMQC_input_bam
    
    output:
    path("${caseID}.${sampleID}.bamqc/")// into qualimapBAMQC
    
    when:
    !params.skipQC && params.qualimap
    
    script:
    //use_bed = qualimap_ROI ? "-gff ${qualimap_ROI}" : ''
    """
    unset DISPLAY
    singularity run -B ${s_bind} ${simgpath}/qualimap.sif qualimap --java-mem-size=5G \
    bamqc \
    -nt ${task.cpus} \
    -outdir ${caseID}.${sampleID}.bamqc \
    -bam ${bam} -gff ${qualimap_ROI} -sd -sdmode 0
    """
}

process rnaseQC {
    tag "$sampleID"
    publishDir "${caseID}/${params.outdir}/QC/RNASeQC", mode: 'copy'
    errorStrategy 'ignore'
    cpus 15

    input:
    tuple val(caseID), val(sampleID), path(bam), path(bai)

    output:
    path("${sampleID}_rnaseqc/*")
    
    when:
    !params.skipQC

    script:
    """
    rnaseqc242 ${gencode_gtf_collapsed} ${bam} \
    --bed ${gencode36_coding_exons} ${sampleID}_rnaseqc
    """
}


process featureCounts {
    tag "$sampleID"
    errorStrategy 'ignore'

    publishDir "${caseID}/${params.outdir}/genecounts/featureCount", mode: 'copy'
    publishDir "/data/shared/projects/tumortarget/expression/featurecounts/", mode: 'copy'
    publishDir "${caseID}/${params.outdir}/QC/", mode: 'copy', pattern: "*.summary"
    cpus 10

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/subread206'
    input:
    tuple val(caseID), val(sampleID), path(bam), path(bai)
    
    output:
    path("*.{featureCounts, featureCounts.summary}")
    
    script:
    """
    featureCounts -p -T ${task.cpus} -s 2 -t exon -a ${gencode_gtf} -g gene_id ${bam} -o ${caseID}.${sampleID}.featureCounts
    """
}

process htseq_count {
    errorStrategy 'ignore'
    tag "$sampleID"
    publishDir "${caseID}/${params.outdir}/genecounts/htseq_count", mode: 'copy'
    publishDir "/data/shared/projects/tumortarget/expression/htseq/", mode: 'copy'
    //publishDir "/data/shared/genomes/hg19/databases/rna_seq/inhouse/genecounts/htseq_count/", mode: 'copy'
    cpus 10

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/htseq'
    input:
    tuple val(caseID),val(sampleID), path(bam), path(bai) 
    
    output:
    path("${params.rundir}.${caseID}.${sampleID}.htseq_count.txt") 
    
    script:
    """
    htseq-count ${bam} ${gencode_gtf} -s reverse > ${params.rundir}.${caseID}.${sampleID}.htseq_count.txt
    """
}

process rsem {
    tag "$sampleID"
    cpus 20
    errorStrategy 'ignore'
    
    publishDir "${caseID}/${params.outdir}/genecounts/rsem", mode: 'copy', pattern: "*.results"
    publishDir "/data/shared/projects/tumortarget/expression/rsem/genecounts", mode: 'copy', pattern: "*.genes.results"
    publishDir "/data/shared/projects/tumortarget/expression/rsem/transcriptcounts", mode: 'copy', pattern: "*.isoforms.results"
    publishDir "${caseID}/${params.outdir}/QC/", mode: 'copy', pattern: "*.${caseID}.${sampleID}.RSEM.stat/*"

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/rsem133'

    input:
    tuple val(caseID), val(sampleID), path(bam)

    output:
    path("*.results")
    tuple val(caseID), val(sampleID), path("*.genes.results"), emit: rsem_tpm_ch
    path ("*.RSEM.stat/*"),emit: rsem_stats_out

    script:
 
    strandedness=params.rsem_strand ? "--strandedness ${params.rsem_strand}":""
    """
    rsem-calculate-expression \
    --alignments --no-bam-output -p ${task.cpus} \
    --calc-ci --append-names \
    $strandedness \
    --paired-end ${bam} ${index_rsem} ${caseID}.${sampleID}.RSEM
    """
}

process rsem_genecount_TPM {
    publishDir "/data/shared/projects/tumortarget/expression/rsem/genecounts", mode: 'copy', pattern: "*.txt"
    publishDir "${caseID}/${params.outdir}/TumorTarget_files/", mode: 'copy', pattern: "*.txt"
    publishDir "${caseID}/${params.outdir}/genecounts/rsem/", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(caseID), val(sampleID), path(genecounts)// from rsem_tpm_ch

    output:
    path("*.txt")
    script:
    """
    cut -f1,6 ${genecounts} > ${sampleID}.RSEM.tpm.2col.txt
    sed -i "1 s/TPM/${sampleID}/" ${sampleID}.RSEM.tpm.2col.txt
    """
}

process QC_multiQC {
    //publishDir "${caseID}/${params.outdir}/", mode: 'copy'
    publishDir "${launchDir}", mode: 'copy'
    //publishDir "${launchDir}/*/${params.outdir}", mode: 'copy'
    input:
    //path("*_fastqc.*") from fastqc_results.collect()
    path("*.{html}") from fastp_results.collect()
    //path("*.preseq.txt") from preseq_output.collect()
    //path("*.{txt,pdf,r,xls}") from rseqc_out.collect()
    //path("*.{txt,pdf,r,xls}") from rseqc_out.collect()
    //path("*.bamqc/*") from qualimapBAMQC.collect()
    path("*.qualimapRNA/*") from qualimapRNASEQ.collect()
    path("*.featureCounts.summary") from featureCounts_output.collect()
    path("*.RSEM.stat/*") from rsem_stats_out.collect()

    //path("bamQC/*") from bamQCReport.collect()
    output:
    path ("*multiQC.report.html")

    when:
    !params.skipQC

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/multiqc.sif \
    -f -q -c ${multiqc_config} \
    -n ${date}.RNA_results.multiQC.report.html \
    ${launchDir}/*/${params.outdir}/*
    """
}

////////////////// ALTERNATIVE SPLICING ///////////////////////////////

process trinitySplicing { 
    errorStrategy 'ignore'
    tag "$caseID"

    cpus 20
    publishDir "${caseID}/${params.outdir}/trinity_splicing", mode: 'copy'
    publishDir "${caseID}/${params.outdir}/TumorTarget_files", mode: 'copy', pattern: "*.{html,cancer.introns}"

    input:
    tuple val(caseID), val(sampleID), path(bam),path(bai), path(junction), path(r1),path(r2),path(sj_tab) //from trinity_splicing_merged_input_ch

    //path(genome_lib_starfusion)

    output:
    path("${caseID}.${sampleID}_trinity_splicing.*"), emit: results
    tuple val(caseID), val(sampleID), path("${caseID}.${sampleID}_trinity_splicing.cancer.introns"), emit: inhouse_list

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/ctat_splicing_v002.sif /usr/local/src/CTAT-SPLICING/STAR_to_cancer_introns.py \
    --SJ_tab_file ${sj_tab} \
    --chimJ_file ${junction} \
    --bam_file ${bam} \
    --output_prefix ${caseID}.${sampleID}_trinity_splicing \
    --sample_name ${caseID}.${sampleID} \
    --min_total_reads 10 \
    --vis \
    --ctat_genome_lib ${genome_lib_starfusion}
    """
}

process splicing_inhouse_list {
    errorStrategy 'ignore'
    tag "$caseID"

    publishDir "${caseID}/${params.outdir}/TumorTarget_files", mode: 'copy', pattern: "*.INHOUSE.txt"
    publishDir "${caseID}/${params.outdir}/trinity_splicing", mode: 'copy'
    
    input:
    tuple val(caseID), val(sampleID), path(trinity_splicing)// from inhouse_list_splicing_ch

    output:
    path("*.INHOUSE.txt")

    script:
    """
    cat ${trinity_splicing} | grep -w -f ${inhouse_splicing_genelist} > ${caseID}.${sampleID}.trinity_splicing.INHOUSE.txt
    """
}

///////////////////////////////////////////////////////////////////////////

///////////////////// FUSIONS ///////////////////////////////////
/*
process arriba {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${params.outdir}/genefusions/arriba", mode: 'copy'
    publishDir "${caseID}/${params.outdir}/TumorTarget_files", mode: 'copy', pattern: "*.{INHOUSEFUSION_V2.txt,pdf}"
    input:
    tuple val(caseID), val(sampleID), path(bam), path(bai)// from arriba_input_bam_ch //.join(arriba_input_index_ch)
    output:
    path("*.{txt,tsv,pdf}"), emit: arriba_out
    tuple val(caseID), val(sampleID), path("${caseID}.${sampleID}.arriba.1col.finspect.txt"),emit: arriba_fusion_inspect
    tuple val(caseID), path("${caseID}.${sampleID}.arriba.fusions.INHOUSEFUSION_V2.txt"), emit: fusions
    tuple val(caseID), path("${caseID}.${sampleID}.arriba.fusions.txt"), emit: all_fusions

    script:
    """
    arriba -x ${bam} -g ${gencode_gtf} -a ${genome_fasta} -b ${arriba_blacklist} -o ${caseID}.${sampleID}.arriba.fusions.txt -O ${caseID}.${sampleID}.arriba.discardedfusions.tsv
    
    draw_fusions.R \
    --fusions=${caseID}.${sampleID}.arriba.fusions.txt \
    --alignments=${bam} \
    --annotation=${gencode_gtf} \
    --cytobands=${arriba_cytoband} \
    --proteinDomains=${arriba_protein_gff} \
    --output=${caseID}.${sampleID}.arribafusionsRplots.pdf

    cut -f1-2 ${caseID}.${sampleID}.arriba.fusions.txt | sed "s/\t/--/g" > ${caseID}.${sampleID}.arriba.1col.finspect.txt

    cat ${caseID}.${sampleID}.arriba.fusions.txt| grep -w -f ${inhouse_genelist} > ${caseID}.${sampleID}.arriba.fusions.INHOUSEFUSION_V2.txt
    """
}
*/
process arriba {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${params.outdir}/genefusions/arriba", mode: 'copy'
    publishDir "${caseID}/${params.outdir}/TumorTarget_files", mode: 'copy', pattern: "*.{INHOUSEFUSION_V2.txt,pdf}"
    input:
    tuple val(caseID), val(sampleID), path(bam), path(bai)// from arriba_input_bam_ch //.join(arriba_input_index_ch)
    output:
    path("*.{txt,tsv,pdf}"), emit: arriba_out
    tuple val(caseID), val(sampleID), path("${caseID}.${sampleID}.arriba.1col.finspect.txt"),emit: arriba_fusion_inspect
    tuple val(caseID), path("${caseID}.${sampleID}.arriba.fusions.INHOUSEFUSION_V2.txt"), emit: fusions
    tuple val(caseID), path("${caseID}.${sampleID}.arriba.fusions.txt"), emit: all_fusions

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/arriba240.sif /arriba_v2.4.0/arriba \
    -x ${bam} -g ${gencode_gtf} -a ${genome_fasta} -b ${arriba_blacklist} -o ${caseID}.${sampleID}.arriba.fusions.txt -O ${caseID}.${sampleID}.arriba.discardedfusions.tsv
    
    singularity run -B ${s_bind} ${simgpath}/arriba240.sif /arriba_v2.4.0/draw_fusions.R \
    --fusions=${caseID}.${sampleID}.arriba.fusions.txt \
    --alignments=${bam} \
    --annotation=${gencode_gtf} \
    --cytobands=${arriba_cytoband} \
    --proteinDomains=${arriba_protein_gff} \
    --output=${caseID}.${sampleID}.arribafusionsRplots.pdf

    cut -f1-2 ${caseID}.${sampleID}.arriba.fusions.txt | sed "s/\t/--/g" > ${caseID}.${sampleID}.arriba.1col.finspect.txt

    cat ${caseID}.${sampleID}.arriba.fusions.txt| grep -w -f ${inhouse_genelist} > ${caseID}.${sampleID}.arriba.fusions.INHOUSEFUSION_V2.txt
    """
}


process starfusion {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${params.outdir}/genefusions/StarFusion/", mode: 'copy'
    publishDir "${caseID}/${params.outdir}/TumorTarget_files", mode: 'copy', pattern: "*.INHOUSEFUSION_V2.txt"

    input:
    tuple val(caseID), val(sampleID), path(junctions)

    output:
    path("*.{tsv,txt}"), emit: results
    tuple val(caseID), path("${caseID}.${sampleID}.StarFusion.INHOUSEFUSION_V2.txt"), emit: fusions 
    tuple val(caseID), path("${caseID}.${sampleID}_abridged.tsv"), emit: all_fusions
    
    script:
    """
    singularity run -B ${s_bind} ${simgpath}/ctat_starfusion_v1_11_1.sif STAR-Fusion \
    --genome_lib_dir ${genome_lib_starfusion} \
    -J ${junctions} \
    --examine_coding_effect \
    --output_dir .
    mv star-fusion.fusion_predictions.tsv ${caseID}.${sampleID}_star-fusion.tsv
    mv star-fusion.fusion_predictions.abridged.tsv ${caseID}.${sampleID}_abridged.tsv
    mv star-fusion.fusion_predictions.abridged.coding_effect.tsv ${caseID}.${sampleID}_abridged.coding_effect.tsv

    cat ${caseID}.${sampleID}_abridged.tsv | grep -w -f /data/shared/genomes/databases/genelists/tumortarget/tumortarget.inhouse.v2.127genes.txt > ${caseID}.${sampleID}.StarFusion.INHOUSEFUSION_V2.txt
    """
}


///data/shared/genomes/databases/genelists/tumortarget/tumortarget.inhouse.v2.127genes.txt

/*
trinityfusion_input_junction // caseID, junctionfile
    .join(trinityfusion_input_bam) //caseID, SampleID, bam, bai
    .join(trinity_splicing_readinput_ch2) // caseID, R1, R2
    .join(star_sjtab_out) //caseID, STAR SJ tab out
    .into {trinityfusion_merged_input_ch; trinity_splicing_merged_input_ch; ID_junction_bam_reads_ch1}


trinityfusion_input_bam //caseID, SampleID, bam, bai
    .join(trinityfusion_input_junction) // caseID, junctionfile
    .join(trinity_splicing_readinput_ch2) // caseID, R1, R2
    .join(star_sjtab_out) //caseID, STAR SJ tab out
    .into {trinityfusion_merged_input_ch; trinity_splicing_merged_input_ch; ID_junction_bam_reads_ch1}
*/



process kallisto {
    errorStrategy 'ignore'
    tag "$caseID"
    input:
    tuple val(caseID), val(sampleID), file(r1), file(r2)
    
    output:
    // TO DO
    script:
    """
    singularity run -B ${s_bind} ${simgpath}/kallisto-0.46.2.sif kallisto quant \
    -i ${kallisto_index} \
    --fusion \
    -o ${sampleID}_kalisto/ \
    ${r1} ${r2}


    """
}
/*
process kallisto_pizzly {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${params.outdir}/genefusions/pizzly/", mode: 'copy'

    conda '/data/shared/programmer/miniconda3/envs/pizzly'
    
    maxForks 5
    input:
    tuple val(caseID), val(sampleID), file(r1), file(r2)
    
    output:
    tuple val(caseID), path("${caseID}_${sampleID}.pizzlyFusions.INHOUSEFUSION_V2.txt"), emit: fusions
    tuple val(caseID), path("${caseID}_${sampleID}.pizzlyFusions.txt"), emit: all_fusions
    path("${caseID}_kallisto/")
    path("${caseID}_${sampleID}.fusions.fasta")
    path("${caseID}_${sampleID}.json")
   
    script:
    """
    singularity run -B ${s_bind} ${simgpath}/kallisto-0.46.2.sif kallisto quant \
    -i ${kallisto_index} \
    --fusion \
    -o ${caseID}_kallisto/ \
    ${r1} ${r2}

    pizzly -k 31 \
    --gtf ${ensemble102_gtf} \
    --cache index.cache.txt \
    --align-score 2 \
    --insert-size 400 \
    --fasta ${ensemble102_transcript_fasta} \
    --output ${caseID}_${sampleID} ${caseID}_kallisto/fusion.txt
    
    pizzly_flatten_json.py ${caseID}_${sampleID}.json > ${caseID}_${sampleID}_pizzlyFusionsRAW.txt
    
    awk 'BEGIN{FS=OFS="\t"} ($5>4||$6>4)' ${caseID}_${sampleID}_pizzlyFusionsRAW.txt > ${caseID}_${sampleID}_pizzlyFusions.txt

    cat ${caseID}_${sampleID}.pizzlyFusions.txt | grep -w -f ${inhouse_genelist} > ${caseID}_${sampleID}_pizzlyFusions.INHOUSEFUSION_V2.txt
    """
}
*/
process kallisto_pizzly {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${params.outdir}/genefusions/pizzly/", mode: 'copy'

    conda '/data/shared/programmer/miniconda3/envs/pizzly'
    
    maxForks 5
    input:
    tuple val(caseID), val(sampleID), file(r1), file(r2)
    
    output:
    tuple val(caseID), path("${caseID}_${sampleID}.pizzlyFusions.INHOUSEFUSION_V2.txt"), emit: fusions
    tuple val(caseID), path("${caseID}_${sampleID}.pizzlyFusions.txt"), emit: all_fusions
    path("${caseID}_kallisto/")
    path("${caseID}_${sampleID}.fusions.fasta")
    path("${caseID}_${sampleID}.json")
   
    shell:
    '''
    singularity run -B !{s_bind} !{simgpath}/kallisto-0.46.2.sif kallisto quant \
    -i !{kallisto_index} \
    --fusion \
    -o !{caseID}_kallisto/ \
    !{r1} !{r2}

    pizzly -k 31 \
    --gtf !{ensemble102_gtf} \
    --cache index.cache.txt \
    --align-score 2 \
    --insert-size 400 \
    --fasta !{ensemble102_transcript_fasta} \
    --output !{caseID}_!{sampleID} !{caseID}_kallisto/fusion.txt
    
    pizzly_flatten_json.py !{caseID}_!{sampleID}.json > !{caseID}_!{sampleID}_pizzlyFusionsRAW.txt
    
    awk 'BEGIN{FS=OFS="\t"} ($5>4||$6>4)' !{caseID}_!{sampleID}_pizzlyFusionsRAW.txt > !{caseID}_!{sampleID}.pizzlyFusions.txt

    cat !{caseID}_!{sampleID}.pizzlyFusions.txt | grep -w -f !{inhouse_genelist} > !{caseID}_!{sampleID}.pizzlyFusions.INHOUSEFUSION_V2.txt
    '''
}

process jaffa_conda {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${params.outdir}/genefusions/jaffa/", mode: 'copy'

    conda '/data/shared/programmer/miniconda3/envs/bpipe'

    cpus 12
    maxForks 5

    input:
    tuple val(caseID), val(sampleID), file(r1), file(r2)

    output:
    tuple val(caseID), path("${caseID}.jaffaFusions.INHOUSEFUSION_V2.csv"), emit: fusions
    tuple val(caseID), path("${caseID}.jaffaFusions.csv"), emit: all_fusions
    path("${caseID}.jaffaFusions.fasta")
    path("${caseID}.jaffaFusions.csv")
    script:
    """
    bpipe run \
    -n ${task.cpus} \
    /data/shared/programmer/JAFFA-version-2.3/JAFFA_direct.groovy \
    ${r1} ${r2}
    mv jaffa_results.csv ${caseID}.jaffaFusions.csv
    mv jaffa_results.fasta ${caseID}.jaffaFusions.fasta 

    cat ${caseID}.jaffaFusions.csv | grep -w -f ${inhouse_genelist} > ${caseID}.jaffaFusions.INHOUSEFUSION_V2.csv
    """
}

process fusioncatcher {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${params.outdir}/genefusions/", mode: 'copy'
    cpus 30
    maxForks 5
    input:
    tuple val(caseID), val(sampleID), file(r1), file(r2)

    output:
    path("fusioncatcher/")
    tuple val(caseID), path("fusioncatcher/${caseID}.${sampleID}.fusioncatcher.fusions.INHOUSEFUSION_V2.txt"), emit: fusions
    tuple val(caseID), path("fusioncatcher/${caseID}.${sampleID}.fusioncatcher.fusions.txt"), emit: all_fusions
    script:
    """
    singularity run -B ${s_bind} ${simgpath}/fusioncatcher-1.33.sif /opt/fusioncatcher/v1.33/bin/fusioncatcher.py \
    -d ${fusioncatcher_db} \
    -i ${r1},${r2} \
    --skip-blat \
    --skip-star \
    -p ${task.cpus} \
    -o fusioncatcher

    mv fusioncatcher/final-list_candidate-fusion-genes.txt fusioncatcher/${caseID}.${sampleID}.fusioncatcher.fusions.txt

    mv fusioncatcher/final-list_candidate-fusion-genes.vcf fusioncatcher/${caseID}.${sampleID}.fusioncatcher.fusions.vcf

    mv fusioncatcher/summary_candidate_fusions.txt fusioncatcher/${caseID}.${sampleID}.fusioncatcher.candidateFusionSummary.txt

    cat fusioncatcher/${caseID}.${sampleID}.fusioncatcher.fusions.txt | grep -w -f ${inhouse_genelist} > fusioncatcher/${caseID}.${sampleID}.fusioncatcher.fusions.INHOUSEFUSION_V2.txt
    """
}
/*
process fusionreport_inhouse {
    errorStrategy 'ignore'
    tag "$caseID"

    publishDir "${caseID}/${params.outdir}/genefusions/fusionReportINHOUSEFUSION_V2/", mode: 'copy'

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/fusionreport'

    input:
    tuple val(caseID), path(arriba_fusions), path(starfusion_fusions), path(fusioncatcher_fusions), path(pizzly_fusions)//,path(jaffa_fusions)
    output:
    path("fusionreport_output/")
    path("${caseID}.FusionReport.INHOUSEFUSION_V2.html")
    script:
    """
    /data/shared/programmer/fusion-report-2.1.5p7/bin/fusion_report run \
    ${caseID} \
    fusionreport_output \
    ${fusionreport_db} \
    --arriba ${arriba_fusions} \
    --starfusion ${starfusion_fusions} \
    --fusioncatcher ${fusioncatcher_fusions} \
    --pizzly ${pizzly_fusions}

    mv fusionreport_output/index.html fusionreport_output/${caseID}.FusionReport.INHOUSEFUSION_V2.html
    cp fusionreport_output/${caseID}.FusionReport.INHOUSEFUSION_V2.html .
    """
}
*/
process fusionreport_inhouse {
    errorStrategy 'ignore'
    tag "$caseID"

    publishDir "${caseID}/${params.outdir}/genefusions/fusionReportINHOUSEFUSION_V2/", mode: 'copy'

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/fusionreport'

    input:
    tuple val(caseID), path(arriba_fusions), path(starfusion_fusions), path(fusioncatcher_fusions), path(pizzly_fusions)//,path(jaffa_fusions)
    output:
    path("fusionreport_output/")
    path("${caseID}.FusionReport.INHOUSEFUSION_V2.html")
    script:
    """
    fusion_report run \
    ${caseID} \
    fusionreport_output \
    ${fusionreport_db} \
    --arriba ${arriba_fusions} \
    --starfusion ${starfusion_fusions} \
    --fusioncatcher ${fusioncatcher_fusions} \
    --pizzly ${pizzly_fusions}

    mv fusionreport_output/index.html fusionreport_output/${caseID}.FusionReport.INHOUSEFUSION_V2.html
    cp fusionreport_output/${caseID}.FusionReport.INHOUSEFUSION_V2.html .
    """
    //
    //    --jaffa ${jaffa_fusions}
}

process fusionreport_full {
    errorStrategy 'ignore'
    tag "$caseID"

    publishDir "${caseID}/${params.outdir}/genefusions/fusionReportALL/", mode: 'copy'

   conda '/lnx01_data3/shared/programmer/miniconda3/envs/fusionreport'

    input:
    tuple val(caseID), path(arriba_fusions), path(starfusion_fusions), path(fusioncatcher_fusions), path(pizzly_fusions)//,path(jaffa_fusions)
    output:
    path("fusionreportALL/")
    path("${caseID}.FusionReport.ALL.html")
    script:
    """
    fusion_report run \
    ${caseID} \
    fusionreportALL \
    ${fusionreport_db} \
    --arriba ${arriba_fusions} \
    --starfusion ${starfusion_fusions} \
    --fusioncatcher ${fusioncatcher_fusions} \
    --pizzly ${pizzly_fusions}

    mv fusionreportALL/index.html fusionreportALL/${caseID}.FusionReport.ALL.html
    cp fusionreportALL/${caseID}.FusionReport.ALL.html .
    """
}
//    
//--jaffa ${jaffa_fusions}
/////////////////////////////////////////////////////////////////
//////////////////// SUBWORKFLOWS ///////////////////////////////
/////////////////////////////////////////////////////////////////

workflow SUB_RNA_PREPROCESS {

    take:
    case_sample_reads_ch

    main:
    inputFiles_symlinks_fq(case_sample_reads_ch)
    fastp_TRIM(case_sample_reads_ch)
    align_STAR(fastp_TRIM.out.fastp_out_ch) // star_out_bam_ch: caseID, sampleID, bam,bai

    emit:
    star_bam=align_STAR.out.star_out_bam_ch // CaseID, sampleID, bam, bai
    star_rsem_bam=align_STAR.out.rsem_input_bam  // CaseID, sampleID, transcriptomeBAM
    star_arriba_bam=align_STAR.out.arriba_input_bam // Above: caseID, sampleID, bam, bai 
    star_junctions=align_STAR.out.star_junction_out   // CaseID, jounctionfile
    star_chimeric_junctions=align_STAR.out.chimeric_junctions_out // CaseID, sampleID, junc.file       
    star_sjtab=align_STAR.out.star_sjtab_out           // CaseID, SJ.out.tab
    trinity_collect=align_STAR.out.trinity_collect     // Arriba log output
}

workflow SUB_RNA_QC {
    take:
    star_bam                     //caseID, sampleID, bam, bai

    main:
    rseqc(star_bam)
    qualimapRNAseq(star_bam)
    qualimapBAMQC(star_bam)
    //rnaseQC(star_bam)

    emit:
    rseqc_out=rseqc.out
    qualimapRNA_out=qualimapRNAseq.out
    qualimapBAMQC_out=qualimapBAMQC.out
    //rnaseQC_out=rnaseQC.out
}

workflow SUB_RNA_EXPRESSION {

    take:
    star_bam
    star_rsem_bam
    main:
    featureCounts(star_bam)
    htseq_count(star_bam)
    rsem(star_rsem_bam)
    rsem_genecount_TPM(rsem.out.rsem_tpm_ch)

    emit:
    featureCount_exp=featureCounts.out
    htseqCount_exp=htseq_count.out
    rsem_exp=rsem_genecount_TPM.out
    rsem_stats=rsem.out.rsem_stats_out
}

workflow SUB_RNA_ALT_SPLICING {

    take:
    trinity_splicing_input

    main:
    trinitySplicing(trinity_splicing_input)
    splicing_inhouse_list(trinitySplicing.out.inhouse_list)


}


workflow SUB_RNA_FUSION {

    take:
    case_sample_reads_ch // raw input channel: caseID, sampleID, r1, r2
    star_arriba_bam     // caseID, sampleID, bam, bai
    star_chimeric_junctions // caseID, sampleID, junction file (for starfusion)

    main:
    arriba(star_arriba_bam)
    //arriba_240(star_arriba_bam)
    fusioncatcher(case_sample_reads_ch)     // TESTING
   // jaffa_conda(case_sample_reads_ch)       // TESTING
    kallisto_pizzly(case_sample_reads_ch)   // TESTING
    starfusion(star_chimeric_junctions)
   
    arriba.out.fusions
        .join(starfusion.out.fusions)
        .join(fusioncatcher.out.fusions)
        .join(kallisto_pizzly.out.fusions)
    //   .join(jaffa_conda.out.fusions)
    .set{fusionreport_input}
   
    arriba.out.all_fusions
        .join(starfusion.out.all_fusions)
        .join(fusioncatcher.out.all_fusions)
        .join(kallisto_pizzly.out.all_fusions)
      //  .join(jaffa_conda.out.all_fusions)
    .set{fusionreport_allFusions_input}
   
    fusionreport_inhouse(fusionreport_input)
    fusionreport_full(fusionreport_allFusions_input)
}



























/*
arribafusions_finspect //caseid, sampleID, arriba_finspect
    .join(starfusion_finspect) //caseid, star_finspect
    .join(fusioninspector_input_reads)  //caseid, r1, r2
    .into {finspect_input_ch1;finspect_input_ch2}
finspect_input_ch2.view()



process fusion_inspector { 
    cpus 50
    publishDir "${caseID}/${params.outdir}/genefusions/fusionInspector", mode: 'copy'
    publishDir "${caseID}/${params.outdir}/TumorTarget_files", mode: 'copy', pattern: "${caseID}.${sampleID}.Finspector/* /*.INHOUSEFUSION_V2.*"

    errorStrategy 'ignore'
    //container 'trinityctat/fusioninspector:latest'

    input:
    path(genome_lib_starfusion)
    tuple val(caseID), val(sampleID), path(arriba), path(starfusion), file(r1), file(r2)// from finspect_input_ch1
    output:
    path("${caseID}.${sampleID}.Finspector/*")
    //path("${caseID}.${sampleID}.Finspector")

    script:
    // trinity removed for now: ,${trinityfusion}
    """
    singularity run -B ${s_bind} ${simgpath}/ctat_fusioninspector_v280.sif FusionInspector \
    --fusions ${arriba},${starfusion} \
    --genome_lib ${genome_lib_starfusion} \
    --CPU ${task.cpus} \
    --left_fq ${r1} --right_fq ${r2} \
    --annotate --vis \
    --min_junction_reads 5 \
    --out_prefix ${caseID}.${sampleID} -O ${caseID}.${sampleID}.Finspector

    cat ${caseID}.${sampleID}.Finspector/${caseID}.${sampleID}.FusionInspector.fusions.abridged.tsv| grep -w -f ${inhouse_genelist} > ${caseID}.${sampleID}.Finspector/${caseID}.${sampleID}.fusionInspector.INHOUSEFUSION_V2.txt
    """
}

*/

