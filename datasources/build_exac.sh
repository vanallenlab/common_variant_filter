#!/bin/bash
# Brendan Reardon
# Van Allen lab, Dana-Farber Cancer Institute
# build_exac.sh, Version 1.4.0
# Last updated 5 January 2018

v=1.4
outname='exac.lite-pass'$v'-r1.txt'

if [! -z $1 ]; then
	exac_vcf=$1
else
	exac_vcf=ExAC.r1.sites.vep.vcf
fi

if [! -z $2 ]; then
	reference_genome=$2
else
	reference_genome=~/Desktop/hg19/Homo_sapiens_assembly19.fasta
fi

if [! -z $3 ]; then
	gatk=$3
else
	gatk=~/Desktop/gatk/GATK-3.8/GenomeAnalysisTK.jar
fi

#java -jar $gatk -T VariantsToTable -R $reference_genome -V $exac_vcf -o $outname -F CHROM -F POS -F ID -F REF -F ALT -F AC -F AN -F AC_POPMAX -F AN_POPMAX -F POPMAX

java -jar $gatk -T VariantsToTable -R $reference_genome -V $exac_vcf -o $outname -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AF -F AC -F AC_AFR -F AC_AMR -F AC_EAS -F AC_FIN -F AC_NFE -F AC_OTH -F AC_SAS -F AN -F AN_AFR -F AN_AMR -F AN_EAS -F AN_FIN -F AN_NFE -F AN_OTH -F AN_SAS
