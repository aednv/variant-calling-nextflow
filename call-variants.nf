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
params.variantCaller = "mpileup"

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
	gzip -c $inputFile > $inputFile.gz
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
		path("*.snps.vcf.gz"), emit: snps
		path("*.indels.vcf.gz"), emit: indels
	
	"""
	module load GATK/gatk-4.1.8.1
	source /share/pkg/condas/2018-05-11/bin/activate && conda activate gatk_4.1.8.1
	gatk --java-options "-Xmx4g" SplitVcfs -I $vcfFile --INDEL_OUTPUT ${params.id}.${variantCaller}.indels.vcf.gz --SNP_OUTPUT ${params.id}.${variantCaller}.snps.vcf.gz
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
		path("*.combined.filtered.vcf.gz")
	
	script:
	if( $variantCaller == 'mpileup' )
		"""
		module load python/3.5.0
		module load bcftools/1.9
		bcftools view -e 'QUAL <= 20 || DP < 2 || DP > 50' $snps -Ov -o snps.filtered.vcf
		bcftools view -e 'DP < 2 || DP > 50 || IMF < 0.1' $indels -Ov -o indels.filtered.vcf
		bcftools concat -n -Oz -o ${params.id}.${variantCaller}.combined.filtered.vcf.gz snps.filtered.vcf indels.filtered.vcf
		"""

	else if( $variantCaller == 'gatk' )
		"""
		module load GATK/gatk-4.1.8.1
		source /share/pkg/condas/2018-05-11/bin/activate && conda activate gatk_4.1.8.1
		gatk VariantFiltration -V $snps -O snps.filtered.gatk.vcf \
			-filter "QD < 2" --filter-name "QD2" \
			-filter "QUAL <= 20" --filter-name "QUAL20"
		gatk VariantFiltration -V $indels -O indels.filtered.gatk.vcf \
			-filter "QD < 2" --filter-name "QD2"
		"""
	else
		error "Invalid variantCaller parameter $variantCaller"
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
	gzipAndSave( callVariants.out )
	
	//filter variants
	splitVariants( callVariants.out )
	filterVariants( splitVariants.out.snps, splitVariants.out.indels )
}
