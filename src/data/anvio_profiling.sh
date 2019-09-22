#!/bin/bash

LOCATION=verrucomicrobia/data/metagenomes/osd
ANVIO=verrucomicrobia/data/processed/anvio
MAPPING=verrucomicrobia/data/processed/mapping

for GENOME in `cat "${LOCATION}"/sample_list.txt`
do
	anvi-profile -c "${ANVIO}"/CONTIGS.db \
		         -i "${MAPPING}"/CC151_full_"${GENOME}".bam \
				 --num-threads 16 \
				 --sample-name "${GENOME}" \
				 -o "${ANVIO}"/"${GENOME}"
done
