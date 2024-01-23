# KG Vejle tumorboard RNA script

# Generel info:
Requires a samplesheet containing 5 columns in specific order (tab separated), without headerline:
1) caseID, 2) NPN normal WES, 3) NPN tumor WES, 4) NPN tumor RNA, 5) PCGR tumor value

Example samplesheet:

johnDoe 112217976652	111184925465	111184925473    23

The script will automatically look for fastq in subfolders at KG Vejle data archive (fastq storage incl novaRuns).

The user can point to a specific folder containing raw data (FastQ) using the --fastq option 
This is only needed if raw data only exists outside the data archive (e.g. if data are in personal folders or at other KG Vejle analysis servers).

# Usage examples:

Analyze the cases in samplesheet.txt, starting with Fastq files, let the script find the relevant fastqfiles at the archive:

        nextflow run KGVejle/tumorboardRNA -r main --samplesheet /path/to/samplesheet.txt

Analyze the cases in samplesheet.txt, but do not call variants:

        nextflow run KGVejle/tumorboardRNA -r main --samplesheet /path/to/samplesheet.txt --skipVariants



# Main options:
        --help                print this help message

        --genome              hg19 or hg38
                                Default: hg38

        --outdir              Select which folder to write output to.
                                Default: RNA_results

        --samplesheet         path to case samplesheet. Can contain multiple patients/cases (each case in a separate line). 

        --server              Select which server the analysis is performed on (kga01 or lnx01)
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
        --skipFusions           Do not call Fusions

# 
