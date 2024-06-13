#!/bin/bash

# Create a directory for the dummy FASTQ files
mkdir -p dummy_fastqs

# Function to create a dummy FASTQ file
create_dummy_fastq() {
    local filename=$1
    echo "@SEQ_ID" > $filename
    echo "GATTTGGGGTTTAAAGGTTTGAGTTAGTAA" >> $filename
    echo "+" >> $filename
    echo "!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65" >> $filename
}

# Create dummy FASTQ files
create_dummy_fastq "dummy_fastqs/sample1_R1.fastq"
create_dummy_fastq "dummy_fastqs/sample1_R2.fastq"
create_dummy_fastq "dummy_fastqs/sample2_R1.fastq"
create_dummy_fastq "dummy_fastqs/sample2_R2.fastq"
create_dummy_fastq "dummy_fastqs/sample3_R1.fastq"
create_dummy_fastq "dummy_fastqs/sample3_R2.fastq"

echo "Dummy FASTQ files created in the 'dummy_fastqs' directory."
