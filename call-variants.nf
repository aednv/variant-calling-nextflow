#!/usr/bin/env nextflow

//syntax version
nextflow.enable.dsl=2

/*
 * Script parameters
 */
//sample id
params.id = "brl2"

//fastq read file path
params.reads = "./sequenceData/brl2.fastq"

//reference genome file path
params.ref = "./referenceGenome/Sviridis_500_v2.1/assembly/Sviridis_500_v2.0.fa.gz"

//variant caller to use ('mpileup' for bcftools_mpileup or 'gatk' for gatk_haplotypecaller)
variantCaller = 'mpileup'

//snpEff database name for variant annotation
snpDb = 'Sviridis.v2'

//snpEff path
snpEffPath = "/project/uma_madelaine_bartlett/AmberDeneve/snpEff/snpEff/snpEff.jar"

/*
 * Processes
 */
process fastQCCheck {
	tag {"fastQCCheck $reads"}
	executor 'lsf'
	queue 'short'
	cpus 1
	time '3h'
	
	publishDir "results", mode: 'copy'
	
	input:
		path reads
	
	output:
		path("*")
	
	"""
	module load fastqc/0.11.5
	fastqc ${reads}
	"""
}

process indexRef {
	tag {"indexRef $refGenomeUnzipped"}
	executor 'lsf'
	queue 'short'
	cpus 1
	time '3h'
	
	publishDir "sequenceData", mode: 'copy'
	
	input:
		path refGenomeUnzipped
	
	output:
		path("*.fai")
	
	"""
	module load samtools/1.9
	samtools faidx $refGenomeUnzipped
	"""
}

process bwaIndex {
	tag {"bwaIndex $refGenome"}
	executor 'lsf'
	queue 'long'
	clusterOptions '-R "rusage[mem=15000]" "span[hosts=1]"'
	cpus 1
	time '10h'
	
	publishDir "sequenceData/bwa_index", mode: 'copy'
	
	input:
		path refGenome
	
	output:
		path("*")
	
	"""
	module load bwa/0.7.5a
	bwa index ${refGenome}
	"""
}

process bwaMemAlign {
	tag {"bwaMemAlign ${params.id}"}
	executor 'lsf'
	queue 'long'
	clusterOptions '-R "rusage[mem=2000]" "span[hosts=1]"'
	cpus 16
	time '24h'
	
	input:
		path refGenome
		path reads
	
	output:
		path("*.bam")
	
	script:
	"""
	module load bwa/0.7.5a
	module load samtools/1.9
	INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
	bwa mem -t 16 \$INDEX ${reads} > ${params.id}.aligned.sam
	samtools view -S -b ${params.id}.aligned.sam > ${params.id}.aligned.bam
	"""	
}

process sortBam {
	tag {"sortBam $bamFile"}
	executor 'lsf'
	queue 'long'
	clusterOptions '-R "rusage[mem=5000]" "span[hosts=1]"'
	cpus 4
	time '10h'
	
	input:
		path bamFile
	
	output:
		path("*.sorted.bam")
	
	"""
	module load samtools/1.9
	samtools sort -@ 4 $bamFile > ${params.id}.aligned.sorted.bam
	"""
}

process indexBam {
	tag {"indexBam $sortedBam"}
	executor 'lsf'
	queue 'short'
	clusterOptions '-R "rusage[mem=5000]" "span[hosts=1]"'
	cpus 4
	time '3h'
	
	input:
		path sortedBam
	
	output:
		path("*.bai")
	
	"""
	module load samtools/1.9
	samtools index -@ 4 $sortedBam
	"""
}

process makeRefDict {
	tag {"makeRefDict $refGenome"}
	executor 'lsf'
	queue 'short'
	clusterOptions '-R "rusage[mem=10000]" "span[hosts=1]"'
	cpus 1
	time '3h'
	
	publishDir "sequenceData", mode: 'copy'
	
	input:
		path refGenome
	
	output:
		path("*.dict")
	
	"""
	module load GATK/gatk-4.1.8.1
	source /share/pkg/condas/2018-05-11/bin/activate && conda activate gatk_4.1.8.1
	gatk --java-options "-Xmx4g" CreateSequenceDictionary -R $refGenome
	"""
}

process markDuplicates {
	tag {"markDuplicates $sortedBam"}
	executor 'lsf'
	queue 'long'
	clusterOptions '-R "rusage[mem=5000]" "span[hosts=1]"'
	cpus 4
	time '10h'
	
	input:
		path sortedBam
	
	output:
		path("*.md.bam")
	
	"""
	module load GATK/gatk-4.1.8.1
	source /share/pkg/condas/2018-05-11/bin/activate && conda activate gatk_4.1.8.1
	gatk --java-options "-Xmx4g" MarkDuplicates -I $sortedBam -O ${params.id}.aligned.sorted.md.bam -M md_metrics.txt
	"""
}

process callVariants {
	tag {"callVariants ${params.id}"}
	executor 'lsf'
	queue 'long'
	clusterOptions '-R "rusage[mem=4000]" "span[hosts=1]"'
	cpus 8
	time '24h'
	
	publishDir "sequenceData", mode: 'copy'
	
	input:
		path refGenome
		path refIndex
		path refDict
		path bamIndex
		path rgMdSortedBamFile
	
	output:
		path("*.vcf")
	
	script:
	if( variantCaller == 'gatk' )
		"""
		module load GATK/gatk-4.1.8.1
		source /share/pkg/condas/2018-05-11/bin/activate && conda activate gatk_4.1.8.1
		gatk --java-options "-Xmx4g" HaplotypeCaller -R $refGenome -I $rgMdSortedBamFile -O ${params.id}.gatk.vcf
		"""
	else if( variantCaller == 'mpileup' )
		"""
		module load python3/3.5.0
		module load bcftools/1.9
		bcftools mpileup --threads 8 --min-BQ 30 --per-sample-mF -f $refGenome -Ou $rgMdSortedBamFile | bcftools call -mv -Ov -o ${params.id}.mpileup.vcf 
		"""
	else
		error "Invalid variantCaller parameter: ${variantCaller}"
	
}

process unzipRef {
	tag {"gunzip $refGenome"}
	input:
		path refGenome
	
	output:
		path("*.fa")
	
	script:
	"""
	NAME="$refGenome"
	NO_GZ=\${NAME%.*}
	gunzip -c $refGenome > \$NO_GZ
	"""
}

process addReadGroups {
	tag {"addReadGroups $mdSortedBamFile"}
	executor 'lsf'
	queue 'long'
	clusterOptions '-R "rusage[mem=5000]" "span[hosts=1]"'
	cpus 4
	time '10h'
	
	input:
		path mdSortedBamFile
	
	output:
		path("*.rg.bam")
	
	"""
	module load GATK/gatk-4.1.8.1
	source /share/pkg/condas/2018-05-11/bin/activate && conda activate gatk_4.1.8.1
	gatk --java-options "-Xmx4g" AddOrReplaceReadGroups -I $mdSortedBamFile -O ${params.id}.aligned.sorted.md.rg.bam -LB lib1 -PL Illumina -PU unknown_barcode -SM brl2
	"""
}

process gzipAndSave {
	tag {"gzipAndSave $inputFile"}
	executor 'lsf'
	queue 'short'
	clusterOptions '-R "rusage[mem=5000]" "span[hosts=1]"'
	cpus 1
	time '1h'
	
	publishDir "sequenceData", mode: 'copy'
	
	input:
		path inputFile
	
	output:
		path("*.gz")
	
	"""
	gzip -c $inputFile > ${inputFile}.gz
	"""
}

process splitVariants {
	tag {"splitVariants $vcfFile"}
	executor 'lsf'
	queue 'short'
	clusterOptions '-R "rusage[mem=5000]" "span[hosts=1]"'
	cpus 1
	time '2h'
	
	publishDir "sequenceData", mode: 'copy'
	
	input:
		path vcfFile
	
	output:
		path("*.snps.vcf"), emit: snps
		path("*.indels.vcf"), emit: indels
	
	"""
	module load GATK/gatk-4.1.8.1
	source /share/pkg/condas/2018-05-11/bin/activate && conda activate gatk_4.1.8.1
	gatk --java-options "-Xmx4g" SplitVcfs -I $vcfFile --INDEL_OUTPUT ${params.id}.${variantCaller}.indels.vcf --SNP_OUTPUT ${params.id}.${variantCaller}.snps.vcf --STRICT false
	"""
}

process filterVariants {
	tag {"filterVariants $snps $indels"}
	executor 'lsf'
	queue 'short'
	clusterOptions '-R "rusage[mem=10000]" "span[hosts=1]"'
	cpus 1
	time '2h'
	
	publishDir "sequenceData", mode: 'copy'
	
	input:
		path snps
		path indels
	
	output:
		path("*.snps.filtered.vcf"), emit: snpsFiltered
		path("*.indels.filtered.vcf"), emit: indelsFiltered
	
	script:
	if( variantCaller == 'mpileup' )
		"""
		module load python3/3.5.0
		module load bcftools/1.9
		#snp specific filtering - calculate allele frequency, filter for homozygous mutations.
		bcftools +fill-tags $snps -- -t AF | bcftools view -e 'AF <= 0.5 || QUAL <= 10 || DP < 2 || DP > 100' -Ov -o ${params.id}.${variantCaller}.snps.filtered.vcf
		#indel speciic filtering - calculate allele frequency, filter for homozygous mutations.
		bcftools +fill-tags $indels -- -t AF | bcftools view -e 'AF <= 0.5 || DP < 2 || DP > 100 || IMF < 0.1' -Ov -o ${params.id}.${variantCaller}.indels.filtered.vcf
		"""
		
	else if( variantCaller == 'gatk' )
		"""
		module load python3/3.5.0
		module load bcftools/1.9
		#allele freq already annotated
		bcftools view $snps -e 'AF <= 0.5 || QUAL <= 10 || QD < 2 || INFO/DP > 100' -Ov -o ${params.id}.${variantCaller}.snps.filtered.vcf
		bcftools view $indels -e 'AF <= 0.5 || QD < 2 || INFO/DP > 100' -Ov -o ${params.id}.${variantCaller}.indels.filtered.vcf
		"""
		
	else
		error "Invalid variantCaller parameter."
}

process mergeVariants {
	tag {"mergeVariants $snpsFiltered $indelsFiltered"}
	executor 'lsf'
	queue 'short'
	clusterOptions '-R "rusage[mem=10000]" "span[hosts=1]"'
	cpus 1
	time '1h'
	
	publishDir "sequenceData", mode: 'copy'
	
	input:
		path snpsFiltered
		path indelsFiltered
	
	output:
		path("*.combined.vcf")
	
	"""
	module load GATK/gatk-4.1.8.1
	source /share/pkg/condas/2018-05-11/bin/activate && conda activate gatk_4.1.8.1
	gatk MergeVcfs -I $snpsFiltered -I $indelsFiltered -O ${params.id}.${variantCaller}.combined.vcf
	"""
}

process graphVariantFreq {
	tag {"graphVariantFreq $combinedVcf"}
	executor 'lsf'
	queue 'short'
	clusterOptions '-R "rusage[mem=10000]" "span[hosts=1]"'
	cpus 1
	time '1h'
	
	publishDir "results", mode: 'copy'
	
	input:
		path combinedVcf
	
	output:
		path("*.pdf")
	
	script:
	"""
	module load R/3.6.1
	module load gcc/8.1.0
	module load R/3.6.1_packages/tidyverse/1.3.0
	#extract chromosome lengths from vcf
	head -100 $combinedVcf | grep -oP '(?<=#contig=<ID=).*?(?=,)' > chr_names.txt
	head -100 $combinedVcf | grep -oP '(?<=length=).*?(?=>)' > chr_lengths.txt
	#make a 2 column tsv file recording the chromosome and coordinate location of each variant
	awk '!/##/' $combinedVcf | awk -v OFS='\t' '{print \$1, \$2}' | sed 's/#//' > variantLocations.tsv
	#start R and make graphs for each chromosome
	variant_graphing.R ${combinedVcf}
	
	"""
}

process annotateSnps {
	tag {"annotateSnps $combinedVcf"}
	executor 'lsf'
	queue 'short'
	clusterOptions '-R "rusage[mem=10000]" "span[hosts=1]"'
	cpus 1
	time '2h'
	
	publishDir "results", mode: 'copy'
	
	input:
		path combinedVcf
	
	output:
		path("*")
		
	"""
	module load snpEff_snpSift/4.3T
	java -Xmx8g -jar $snpEffPath $snpDb $combinedVcf > ${params.id}.${variantCaller}.combined.ann.vcf
	cat ${params.id}.${variantCaller}.combined.ann.vcf | grep 'HIGH' > ${params.id}.${variantCaller}.combined.ann.high.vcf
	"""
}

process removeHighDensityVar {
	tag {"removeHighDensity $combinedVcf"}
	executor 'lsf'
	queue 'short'
	cpus 1
	time '2h'
	
	publishDir "results", mode: 'copy'
	
	input:
		path combinedVcf
	
	output:
		path("isec_dir/*.vcf")
	
	"""
	module load bcftools/1.9
	filterHighDensitySnps.sh $combinedVcf
	"""
}

/*
 * Workflow
 */
workflow {
	//create channels
	reads_ch = Channel.fromPath( params.reads, checkIfExists: true )
	ref_ch = Channel.fromPath( params.ref, checkIfExists: true )
	
	//run fastQC and prep file indices
	fastQCCheck( reads_ch )
	bwaIndex( ref_ch )
	unzipRef( ref_ch )
	indexRef( unzipRef.out )
	makeRefDict( ref_ch )
	
	//align reads to genome
	bwaMemAlign( bwaIndex.out, reads_ch )
	
	//sort, index, mark duplicate reads, and set read groups in the bam file
	sortBam( bwaMemAlign.out )
	markDuplicates( sortBam.out )
	addReadGroups( markDuplicates.out )
	indexBam( addReadGroups.out )
	gzipAndSave( addReadGroups.out )
	
	//call variants using HaplotypeCaller
	callVariants( unzipRef.out, indexRef.out, makeRefDict.out, indexBam.out, addReadGroups.out )
	
	//filter variants
	splitVariants( callVariants.out )
	filterVariants( splitVariants.out.snps, splitVariants.out.indels )
	mergeVariants( filterVariants.out.snpsFiltered, filterVariants.out.indelsFiltered )
	removeHighDensityVar( mergeVariants.out )
	//graph homozygous variant frequency by chromosome
	graphVariantFreq( removeHighDensityVar.out )
	
	//snpEff annotation
	annotateSnps( removeHighDensityVar.out )
}
