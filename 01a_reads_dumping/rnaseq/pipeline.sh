#!/bin/bash

## Setup
##-----------------------------------------------------------------------------
# Define SRA accession ID and output directory

SRA_ID="$1"
OUTPUT_DIR="${SRA_ID}"
GENOME_DIR="GRCz11/STAR"
GTF_DIR="GRCz11/Danio_rerio.11.110.gtf"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
##-----------------------------------------------------------------------------


# Step 1: Download SRA data using fastq-dump and determine read layout
##-----------------------------------------------------------------------------
prefetch "$SRA_ID" \
    -O "$OUTPUT_DIR"

fasterq-dump "$SRA_ID" \
    --outdir "$OUTPUT_DIR" 

# Check if single-end or paired-end reads were downloaded
if [ -e "$OUTPUT_DIR"/"${SRA_ID}_1.fastq" ] && \
    [ -e "$OUTPUT_DIR"/"${SRA_ID}_2.fastq" ]; then
    READ_LAYOUT="paired-end"
elif [ -e "$OUTPUT_DIR"/"${SRA_ID}.fastq" ]; then
    READ_LAYOUT="single-end"
else
    echo "Error: Unable to determine the read layout."
    exit 1
fi

# Remove SRA file for space saving
rm "${OUTPUT_DIR}"/"${SRA_ID}"/"${SRA_ID}.sra"
mkdir -p "${OUTPUT_DIR}/quality_reports"
##-----------------------------------------------------------------------------


# Step 2: Perform quality control with FastQC
##-----------------------------------------------------------------------------
if [ "$READ_LAYOUT" == "single-end" ]; then
    fastqc "$OUTPUT_DIR"/"${SRA_ID}.fastq" \
        -o "$OUTPUT_DIR"/quality_reports
elif [ "$READ_LAYOUT" == "paired-end" ]; then
    fastqc "$OUTPUT_DIR"/"${SRA_ID}_1.fastq" \
        "$OUTPUT_DIR"/"${SRA_ID}_2.fastq" \
        -o "$OUTPUT_DIR"/quality_reports
fi
##-----------------------------------------------------------------------------


# Step 3: Run MultiQC to aggregate FastQC reports
##-----------------------------------------------------------------------------
multiqc "$OUTPUT_DIR"/quality_reports \
    -o "$OUTPUT_DIR"/quality_reports
##-----------------------------------------------------------------------------


# Step 4: Perform quality control and low-quality read removal with fastp
##-----------------------------------------------------------------------------
if [ "$READ_LAYOUT" == "single-end" ]; then
    fastp \
        -i "$OUTPUT_DIR"/"${SRA_ID}.fastq" \
        -o "$OUTPUT_DIR"/"${SRA_ID}_cleaned.fastq" \
        --report_title "${SRA_ID} FastP Report" \
        --html "$OUTPUT_DIR"/"${SRA_ID}_fastp_report.html"
elif [ "$READ_LAYOUT" == "paired-end" ]; then
    fastp \
        -i "$OUTPUT_DIR"/"${SRA_ID}_1.fastq" \
        -I "$OUTPUT_DIR"/"${SRA_ID}_2.fastq" \
        -o "$OUTPUT_DIR"/"${SRA_ID}_cleaned_1.fastq" \
        -O "$OUTPUT_DIR"/"${SRA_ID}_cleaned_2.fastq" \
        --report_title "${SRA_ID} FastP Report" \
        --html "$OUTPUT_DIR"/"${SRA_ID}_fastp_report.html"
fi
##-----------------------------------------------------------------------------


# Step 5: Perform alignment with STAR
##-----------------------------------------------------------------------------
# Modify the reference genome and alignment parameters as needed
if [ "$READ_LAYOUT" == "single-end" ]; then
    STAR \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$OUTPUT_DIR"/"${SRA_ID}_cleaned.fastq" \
        --outFileNamePrefix "$OUTPUT_DIR"/"${SRA_ID}"
elif [ "$READ_LAYOUT" == "paired-end" ]; then
    STAR \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$OUTPUT_DIR"/"${SRA_ID}_cleaned_1.fastq" \
        "$OUTPUT_DIR"/"${SRA_ID}_cleaned_2.fastq" \
        --outFileNamePrefix "$OUTPUT_DIR"/"${SRA_ID}"
fi
##-----------------------------------------------------------------------------


# Step 6: Transform SAM file to BAM
##-----------------------------------------------------------------------------
samtools view -S \
    -bo "$OUTPUT_DIR"/"${SRA_ID}.bam" \
    "$OUTPUT_DIR"/"${SRA_ID}"Aligned.out.sam 
##-----------------------------------------------------------------------------


## Step 7: Sort BAM file
##-----------------------------------------------------------------------------
samtools sort "$OUTPUT_DIR"/"${SRA_ID}.bam" \
    -o "$OUTPUT_DIR"/"${SRA_ID}.sorted.bam"
##-----------------------------------------------------------------------------


# Step 8: Count reads with featureCounts
##-----------------------------------------------------------------------------
# Modify annotation file and parameters as needed
if [ "$READ_LAYOUT" == "single-end" ]; then
    featureCounts -a "$GTF_DIR" \
        -o "$OUTPUT_DIR"/"${SRA_ID}_counts.txt" \
        "$OUTPUT_DIR"/"${SRA_ID}.sorted.bam" \
        -F GTF -t exon -g gene_id
elif [ "$READ_LAYOUT" == "paired-end" ]; then
    featureCounts -a "$GTF_DIR" \
        -o "$OUTPUT_DIR"/"${SRA_ID}_counts.txt" \
        "$OUTPUT_DIR"/"${SRA_ID}.sorted.bam" \
        -F GTF -t exon -g gene_id -p
fi
##-----------------------------------------------------------------------------


# Optional: Clean up intermediate files if needed
rm "$OUTPUT_DIR"/*.fastq
rm "$OUTPUT_DIR"/*.sam
rm "$OUTPUT_DIR"/"${SRA_ID}.bam"

