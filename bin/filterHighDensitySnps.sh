#!/bin/bash

#this script filters out high density variants (more than 10 variants found within 20kb of the variant).

#get input
vcf_file=$1

#convert to bcf and index
bcftools view $1 -O b -o temp.bcf
bcftools index temp.bcf

#make a temp file for reading through variants
awk '!/#/' $1 > temp.vcf
#add the vcf header to new file
awk '/#/' $1 > highDensityVariants.vcf

while read -r c p t
do
	chr=$c
	pos=$p
	echo "$pos"
	pos_up=$(( pos - 20000 ))
	pos_down=$(( pos + 20000 ))
	echo "$chr:$pos_up-$pos_down"
	variants_in_range=$(bcftools view temp.bcf -r $chr:$pos_up-$pos_down | awk '!/#/' | wc -l)
	echo "$variants_in_range"
	if [ $variants_in_range -gt 10 ]
	then
		echo "true"
		bcftools view temp.bcf -r $chr:$pos -O v | awk '!/#/' >> highDensityVariants.vcf
	fi
done < temp.vcf

bcftools view highDensityVariants.vcf -O b -o highDensityVariants.bcf
bcftools index highDensityVariants.bcf

#intersect the two files
bcftools isec -C temp.bcf highDensityVariants.bcf -p isec_dir

