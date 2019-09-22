#!/bin/bash

LOCATION=verrucomicroba2/data/raw/osd
PROCESSED=verrucomicroba2/data/processed
METAGENOMES=verrucomicroba2/data/raw/osd

for GENOME in `cat "${LOCATION}"/sample_list.txt`;
do
# do the bowtie mapping to get the SAM file:
bowtie2 -x "${PROCESSED}"/03mapping/CC151_full \
		-1 "${METAGENOMES}"/"${GENOME}"_R1.qc.fastq.gz \
		-2 "${METAGENOMES}"/"${GENOME}"_R2.qc.fastq.gz \
		-S "${PROCESSED}"/03mapping/CC151_full_"${GENOME}".sam \
		--threads 20

# covert the resulting SAM file to a BAM file:
samtools view -F 4 -bS "${PROCESSED}"/03mapping/CC151_full_"${GENOME}".sam > "${PROCESSED}"/03mapping/CC151_full_"${GENOME}".sam-RAW.bam

# sort and index the BAM file:
samtools sort "${PROCESSED}"/03mapping/CC151_full_"${GENOME}".sam-RAW.bam -o "${PROCESSED}"/03mapping/CC151_full_"${GENOME}".bam

samtools index "${PROCESSED}"/03mapping/CC151_full_"${GENOME}".bam;

# remove temporary files:
#rm $LOCATION.sam $LOCATION-RAW.bam;
done	
