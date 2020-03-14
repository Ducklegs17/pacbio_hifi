#!/bin/bash

MAX_THREADS = 32
SAMPLES = ["SRR99694"]
PREFIXES = ["SRR9969479.m64011_190430_000122", "SRR9969480.m54119U_190503_125238"]
BBDUK_KMER = ["31"]
BBDUK_COVERAGE = ["0.9"]
CHUNKS = 20
READ_DEPTH = ["10","25","50","75","100"]
LENGTH_CHLOROPLAST = ["134502"]
LENGTH_MITOCHONDRIA = ["415805"]
LENGTH_GENOME = ["5500000"]
ASSEMBLY_TYPE = ["chloroplast"]
ASSEMBLY_TOOLS = ["flye","hifiasm","canu","hicanu"]
READ_SELECTION_METHOD = ["longest"]
HIFIASM_PATH = "../tools/hifiasm/hifiasm"
HICANU_PATH = "/fast/users/a1761942/tools/canu/Linux-amd64/bin/canu"
GFA_TYPES = ["a_ctg","p_ctg","p_utg","r_utg"]

localrules:
	all,
	canu,
	hicanu,

rule all:
	input:
		expand("0_raw/{PREFIX}.subreads.bam.pbi", PREFIX = PREFIXES),
		expand("1_subreads/chunks/{PREFIX}.ccs.{NUM}.bam", NUM = range(1,CHUNKS+1), PREFIX = PREFIXES),
		expand("1_subreads/{PREFIX}.ccs.bam", PREFIX = PREFIXES),
		expand("1_subreads/{SAMPLE}.ccs.bam", SAMPLE = SAMPLES),
		expand("2_{ASS_TYPE}_reads/unsorted/{SAMPLE}_{KMER}_{COV}.fasta", ASS_TYPE = ASSEMBLY_TYPE, SAMPLE = SAMPLES, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE),
		expand("2_{ASS_TYPE}_reads/sorted/{SAMPLE}_{KMER}_{COV}.fasta", ASS_TYPE = ASSEMBLY_TYPE, SAMPLE = SAMPLES, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE),
		expand("3_{ASS_TYPE}_subset/{READ_SELECT}/{SAMPLE}_{KMER}_{COV}_{DEPTH}.fasta", ASS_TYPE = ASSEMBLY_TYPE, READ_SELECT = READ_SELECTION_METHOD, SAMPLE = SAMPLES, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE, DEPTH = READ_DEPTH),
		expand("4_{ASS_TYPE}_assembly/{TOOL}/{SAMPLE}/{READ_SELECT}_{KMER}_{COV}_{DEPTH}/assembly.fasta", ASS_TYPE = ASSEMBLY_TYPE, TOOL = ASSEMBLY_TOOLS, SAMPLE = SAMPLES, READ_SELECT = READ_SELECTION_METHOD, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE, DEPTH = READ_DEPTH),
		expand("4_{ASS_TYPE}_assembly/hifiasm/{SAMPLE}/{READ_SELECT}_{KMER}_{COV}_{DEPTH}/assembly.{GFA_TYPE}.fasta",ASS_TYPE = ASSEMBLY_TYPE, SAMPLE = SAMPLES, READ_SELECT = READ_SELECTION_METHOD, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE, DEPTH = READ_DEPTH, GFA_TYPE = GFA_TYPES),
		expand("quast/assemblytype_{ASS_TYPE}_assemblytool_{TOOL}_readselect_{READ_SELECT}_prefix_{SAMPLE}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}/report.tsv", ASS_TYPE = ASSEMBLY_TYPE, TOOL = ASSEMBLY_TOOLS, READ_SELECT = READ_SELECTION_METHOD, SAMPLE = SAMPLES, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE, DEPTH = READ_DEPTH),
		expand("quast_gfa/assemblytype_{ASS_TYPE}_assemblytool_hifiasmgfa_readselect_{READ_SELECT}_prefix_{SAMPLE}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}_gfatype_{GFA_TYPE}/report.tsv", ASS_TYPE = ASSEMBLY_TYPE, READ_SELECT = READ_SELECTION_METHOD, SAMPLE = SAMPLES, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE, DEPTH = READ_DEPTH, GFA_TYPE = GFA_TYPES),
		expand("mummer/{ASS_TYPE}/prefix_{SAMPLE}_assemblytool_{TOOL}_readselect_{READ_SELECT}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}.mums", ASS_TYPE = ASSEMBLY_TYPE, SAMPLE = SAMPLES, TOOL = ASSEMBLY_TOOLS, READ_SELECT = READ_SELECTION_METHOD, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE, DEPTH = READ_DEPTH),
		expand("mummer/{ASS_TYPE}/prefix_{SAMPLE}_assemblytool_hifiasmgfa_gfatype_{GFA_TYPE}_readselect_{READ_SELECT}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}.mums", ASS_TYPE = ASSEMBLY_TYPE, SAMPLE = SAMPLES, READ_SELECT = READ_SELECTION_METHOD, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE, DEPTH = READ_DEPTH, GFA_TYPE = GFA_TYPES),

#creates an index file of the raw hifi reads
rule pbindex:
	input:
		"0_raw/{prefix, SRR[0-9]*\.\m[0-9]*\.[0-9]*}.subreads.bam",	
	output:
		"0_raw/{prefix}.subreads.bam.pbi",
	log:
		"logs/pbindex/{prefix}.log",
	benchmark:
		"benchmarks/pbindex/prefix_{prefix}.tsv",
	threads:
		1	
	conda:
		"envs/pacbio.yaml"
	shell:
		"""
		pbindex {input}
		"""

# generate high quality consensus reads
rule ccs:
	input:
		sub = "0_raw/{prefix, SRR[0-9]*\.\m[0-9]*\.[0-9]*}.subreads.bam",
		pbi = "0_raw/{prefix}.subreads.bam.pbi",
	output:
		"1_subreads/chunks/{prefix}.ccs.{num}.bam",
	log:
		"logs/ccs/{prefix}_{num}.log"
	benchmark:
		"benchmarks/ccs/prefix_{prefix}_chunknum_{num}.tsv"	
	threads:
		MAX_THREADS
	params:
		chnk = CHUNKS
	conda:
		"envs/pacbio.yaml",
	shell:
		"""
		ccs {input.sub} {output} --chunk {wildcards.num}/{params.chnk} -j {threads}
		"""

#merge chunks
rule pbmerge:
	input:
		expand("1_subreads/chunks/{{prefix}}.ccs.{num}.bam", prefix = PREFIXES, num = range(1,CHUNKS+1)),
	output:
		"1_subreads/{prefix}.ccs.bam",
	log:
		"logs/pbmerge/{prefix}.log",
	benchmark:
		"benchmarks/pbmerge/prefix_{prefix}.log",
	threads:
		MAX_THREADS
	conda:
		"envs/pacbio.yaml",
	shell:
		"""
		pbmerge -o {output} {input}
		"""

#merge bam files
rule pbmerge_sample:
	input:
		expand("1_subreads/{prefix}.ccs.bam", prefix = PREFIXES),
	output:
		"1_subreads/{sample}.ccs.bam",
	log:
		"logs/pbmergesample/{sample}.log",
	benchmark:
		"benchmarks/pbmergesample/prefix_{sample}.tsv",
	threads:
		1
	conda:
		"envs/pacbio.yaml",
	shell:
		"""
		pbmerge -o {output} {input}
		"""	

rule bbduk:
	input:
		bam = "1_subreads/{sample}.ccs.bam",
		ref = "reference/{ass_type}.fasta",
	output:
		out = "2_{ass_type}_reads/unsorted/{sample}_{kmer}_{cov}.fasta",
	log:
		"logs/bbduk/{sample}_{ass_type}_{kmer}_{cov}.log",
	benchmark:
		"benchmarks/bbduk/assemblytype_{ass_type}_prefix_{sample}_kmer_{kmer}_cov_{cov}.tsv",
	threads:
		5
	conda:
		"envs/bbmap.yaml",
	shell:
		"""
		echo "Extracting {wildcards.ass_type} reads with bbduk..."
		bbduk.sh in={input.bam} outm={output.out} ref={input.ref} threads={threads} k={wildcards.kmer} -Xmx1g mincovfraction={wildcards.cov}
		"""

#subset the matching reads down to specified coverage
rule seqtk:
	input:
		"2_{ass_type}_reads/unsorted/{sample}_{kmer}_{cov}.fasta",
	output:
		"3_{ass_type}_subset/random/{sample}_{kmer}_{cov}_{depth}.fasta",
	log:
		"logs/seqtk/random_{ass_type}_{sample}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/seqtk/assemblytype_{ass_type}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		1
	params:
		chlor = LENGTH_CHLOROPLAST,
		mito = LENGTH_MITOCHONDRIA,
	conda:
		"envs/seqtk.yaml",
	shell:
		"""		
		if [ {wildcards.ass_type} == "chloroplast" ]; then
			LENGTH={params.chlor}
		fi
		if [ {wildcards.ass_type} == "mitochondria" ]; then
			LENGTH={params.mito}
		fi

		BP_READS="$(grep -v "^>" {input} | wc | awk "{{print \$3-\$1}}")"
                NO_READS="$(grep '>' {input} | wc -l)"
                echo "Total number of reads: ${{NO_READS}}"
                AVG="$(( ${{BP_READS}} / ${{NO_READS}} ))"
                echo "AVG length of reads: ${{AVG}}"
                NUM="$(( ${{LENGTH}} * {wildcards.depth} / ${{AVG}} ))"
                echo "Number of reads required to achieve depth of {wildcards.depth}: ${{NUM}}"
                MAX_DEPTH="$(( ${{BP_READS}} / ${{LENGTH}} ))"
                if [ ${{NUM}} -gt ${{NO_READS}} ]; then
                        NUM=${{NO_READS}}
                        echo "Number of reads required for requested coverage is greater than the number of reads available."
                        echo "All reads will be used, giving a coverage depth of ${{MAX_DEPTH}}x"
                fi
                echo "Subsetting ..."
                seqtk sample -s100 {input} ${{NUM}} > {output}
		"""

rule bbmap_sort:
	input:
		"2_{ass_type}_reads/unsorted/{sample}_{kmer}_{cov}.fasta",
	output:
		"2_{ass_type}_reads/sorted/{sample}_{kmer}_{cov}.fasta",
	log:
		"logs/bbmap_sort/{ass_type}{sample}_{kmer}_{cov}.log",
	benchmark:
		"benchmarks/bbmapsort/assemblytype_{ass_type}_prefix_{sample}_kmer_{kmer}_cov_{cov}.tsv",
	threads:
		1
	conda:
		"envs/bbmap.yaml",
	shell:
		"""
		sortbyname.sh in={input} out={output} name=f length=t ascending=f -Xmx1g
		"""

rule select_longest:
	input:
		"2_{ass_type}_reads/sorted/{sample}_{kmer}_{cov}.fasta",
	output:
		"3_{ass_type}_subset/longest/{sample}_{kmer}_{cov}_{depth}.fasta",
	log:
		"logs/select_longest/{ass_type}_{sample}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/selectlongest/assemblytype_{ass_type}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		1
	params:
		mito = LENGTH_MITOCHONDRIA,
		chlor = LENGTH_CHLOROPLAST,
		genome = LENGTH_GENOME,
	shell:
		"""
		if [ {wildcards.ass_type} == "chloroplast" ]; then
			LEN={params.chlor}
		elif [ {wildcards.ass_type} == "mitochondria" ]; then
			LEN={params.mito}
		elif [ {wildcards.ass_type} == "genome" ]; then
			LEN={params.genome}
		fi

		BP_READS="$(grep -v "^>" {input} | wc | awk "{{print \$3-\$1}}")"
		NO_READS="$(grep '>' {input} | wc -l)"
		echo "Total number of reads: ${{NO_READS}}"
		AVG="$(( ${{BP_READS}} / ${{NO_READS}} ))"
		echo "AVG length of reads: ${{AVG}}"
		NUM="$(( ${{LEN}} * {wildcards.depth} / ${{AVG}} ))"
		echo "Number of reads required to achieve depth of {wildcards.depth}: ${{NUM}}"
		MAX_DEPTH="$(( ${{BP_READS}} / ${{LEN}} ))"
		if [ ${{NUM}} -gt ${{NO_READS}} ]; then
			NUM=${{NO_READS}}
			echo "Number of reads required for requested coverage is greater than the number of reads available."
			echo "All reads will be used, giving a coverage depth of ${{MAX_DEPTH}}x"
		fi
		echo "Subsetting ..."
		awk "/^>/ {{n++}} n>${{NUM}} {{exit}} {{print}}" {input} > {output}

		"""

#Assemble with flye
rule flye:
	input:
		file = "3_{ass_type}_subset/{read_select}/{sample}_{kmer}_{cov}_{depth}.fasta",
	output:
		ass = "4_{ass_type}_assembly/flye/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	log:
		"logs/flye/{ass_type}_{sample}_{read_select}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/flye/assemblytype_{ass_type}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		5
	params:
		chlor = LENGTH_CHLOROPLAST,
		mito = LENGTH_MITOCHONDRIA,
		out = "4_{ass_type}_assembly/flye/{sample}/{read_select}_{kmer}_{cov}_{depth}",
	conda:
		"envs/flye.yaml",
	shell:
		"""
		if [ {wildcards.ass_type} == "chloroplast" ]; then
			SIZE={params.chlor}
		fi
		
		if [ {wildcards.ass_type} == "mitochondria" ]; then
			SIZE={params.mito}
		fi

		flye --pacbio-corr {input.file} --genome-size ${{SIZE}} --out-dir {params.out} --threads {threads}
		"""

#Add a check for genome size and increase time allocation for larger genomes.
rule canu:
	input:
		"3_{ass_type}_subset/{read_select}/{sample}_{kmer}_{cov}_{depth}.fasta",
	output:
		"4_{ass_type}_assembly/canu/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",		
	log:
		"logs/canu/{ass_type}_{sample}_{read_select}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/canu/assemblytype_{ass_type}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		1
	params:
		prefix = "assembly",
		dir = "4_{ass_type}_assembly/canu/{sample}/{read_select}_{kmer}_{cov}_{depth}",
		chlor = LENGTH_CHLOROPLAST,
		mito = LENGTH_MITOCHONDRIA,
		genome = LENGTH_GENOME,
		jobName = "{depth}{kmer}{cov}",
	conda:
		"envs/canu.yaml",
	shell:
		"""

		cp canuFailure.sh {params.dir}

		if [ {wildcards.ass_type} == 'chloroplast' ]; then
			SIZE={params.chlor}
			JTIME="--time=00:40:00"
			echo "Starting chloroplast assembly ..."
		fi
		if [ {wildcards.ass_type} == 'mitochondria' ]; then
			SIZE={params.mito}
			JTIME="--time=00:40:00"
			echo "Starting mitochondria assembly ..."
		fi	
		if [ {wildcards.ass_type} == 'genome' ]; then
			SIZE={params.genome}
			JTIME="--time=04:00:00"
			echo "Starting genome assembly ..."
		fi
		
		(canu -p {params.prefix} -d {params.dir} genomeSize=${{SIZE}} gridOptions="--time=00:40:00" gridOptionsJobName={params.jobName} gridOptionsExecutive="--time=00:10:00" executiveMemory=1 -pacbio-hifi {input} stopOnLowCoverage=5 onSuccess=touch onFailure=./canuFailure.sh) 2> {log}

		while :
		do
			if [ -f {params.dir}/{params.prefix} ]; then
				echo "Canu job {params.jobName} is finished!"
				mv {params.dir}/assembly.contigs.fasta {params.dir}/assembly.fasta
				break
			fi
			if [ -f {params.dir}/failure ]; then
				echo "Canu job {params.jobName} failed. :("
				break
			fi
			sleep 5m
		done

		"""
#when installing hicanu from github, get path to run from the end of the installation text
rule hicanu:
	input:
		"3_{ass_type}_subset/{read_select}/{sample}_{kmer}_{cov}_{depth}.fasta",
	output:
		"4_{ass_type}_assembly/hicanu/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	log:
		"logs/hicanu/{ass_type}_{sample}_{read_select}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/hicanu/assemblytype_{ass_type}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	params:
		hicanu = HICANU_PATH,
		prefix = "assembly",
		dir = "4_{ass_type}_assembly/hicanu/{sample}/{read_select}_{kmer}_{cov}_{depth}",
		chlor = LENGTH_CHLOROPLAST,
		mito = LENGTH_MITOCHONDRIA,
		genome = LENGTH_GENOME,
		jobName = "{depth}{kmer}{cov}",
	shell:
		"""	
		cp canuFailure.sh {params.dir}

		if [ {wildcards.ass_type} == 'chloroplast' ]; then
			SIZE={params.chlor}
			JTIME="--time=00:40:00"
			echo "Starting HiCanu chloroplast assembly ..."
		fi
		if [ {wildcards.ass_type} == 'mitochondria' ]; then
			SIZE={params.mito}
			JTIME="--time=00:40:00"
			echo "Starting HiCanu mitochondria assembly ..."
		fi	
		if [ {wildcards.ass_type} == 'genome' ]; then
			SIZE={params.genome}
			JTIME="--time=04:00:00"
			echo "Starting HiCanu genome assembly ..."
		fi
		
		({params.hicanu} -p {params.prefix} -d {params.dir} genomeSize=${{SIZE}} gridOptions=${{JTIME}} gridOptionsJobName={params.jobName} gridOptionsExecutive="--time=00:10:00" executiveMemory=1 -pacbio-hifi {input} stopOnLowCoverage=5 onSuccess=touch onFailure=./canuFailure.sh) 2> {log}

		while :
		do
			if [ -f {params.dir}/{params.prefix} ]; then
				echo "HiCanu job {params.jobName} is finished!"
				mv {params.dir}/assembly.contigs.fasta {params.dir}/assembly.fasta
				break
			fi
			if [ -f {params.dir}/failure ]; then
				echo "HiCanu job {params.jobName} failed."
				break
			fi
			sleep 5m
		done
		"""


rule hifiasm:
	input:
		"3_{ass_type}_subset/{read_select}/{sample}_{kmer}_{cov}_{depth}.fasta",
	output:
		main = "4_{ass_type}_assembly/hifiasm/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
		actg = "4_{ass_type}_assembly/hifiasm/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.a_ctg.gfa",
		pctg = "4_{ass_type}_assembly/hifiasm/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.p_ctg.gfa",
		putg = "4_{ass_type}_assembly/hifiasm/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.p_utg.gfa",
		rutg = "4_{ass_type}_assembly/hifiasm/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.r_utg.gfa",
	log:
		"logs/hifiasm/{ass_type}_{sample}_{read_select}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/hifiasm/assemblytype_{ass_type}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	params:	
		prefix = "assembly",
		hifiasm = HIFIASM_PATH,
		dir = "4_{ass_type}_assembly/hifiasm/{sample}/{read_select}_{kmer}_{cov}_{depth}",
#	shadow:
#		"shallow"
	threads:
		10
	shell:
		"""
		{params.hifiasm} -o {params.prefix} -t {threads} {input}
		mv {params.prefix}* {params.dir}/
		mv {params.dir}/{params.prefix}.ec.fa {params.dir}/{params.prefix}.fasta
		"""

rule gfa2fa:
	input:
		"4_{ass_type}_assembly/hifiasm/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.{gfa_type}.gfa",
	output:
		"4_{ass_type}_assembly/hifiasm/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.{gfa_type}.fasta",
	log:
		"logs/gfa2fa/{ass_type}_{sample}_{read_select}_{kmer}_{cov}_{depth}_{gfa_type}.log",
	benchmark:
		"benchmarks/gfa2fa/assemblytype_{ass_type}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}_gfatype_{gfa_type}.tsv"
	threads:
		1
	shell:
		"""
		awk '/^S/{{print ">" $2 "\\n" $3}}' {input} | fold > {output}
		"""

rule hifiasm_mummer:
	input:
		ref = "reference/{ass_type}.fasta",
		query = "4_{ass_type}_assembly/hifiasm/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.{gfa_type}.fasta",
	output:
		"mummer/{ass_type}/prefix_{sample}_assemblytool_hifiasmgfa_gfatype_{gfa_type}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.mums",
	log:
		"logs/mummer/{ass_type}_{sample}_hifiasmgfa_{read_select}_{kmer}_{cov}_{depth}_{gfa_type}.log",
	benchmark:
		"logs/mummer/assemblytype_{ass_type}_prefix_{sample}_assemblytool_hifiasmgfa_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}_gfatype_{gfa_type}.tsv",
	params:
		pref = "prefix_{sample}_assemblytool_hifiasmgfa_gfatype_{gfa_type}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}",
	shadow:
		"shallow"
	threads:
		1
	conda:
		"envs/mummer.yaml",
	shell:
		"""
		mummer -mum -b -c {input.ref} {input.query} > {params.pref}.mums
		mummerplot --postscript --prefix={params.pref} {params.pref}.mums
		nucmer -maxmatch -p {params.pref} {input.ref} {input.query}
		show-coords -r -c -H -d -o -T -l {params.pref}.delta > {params.pref}.coords
		show-snps -C -l -r -T -H {params.pref}.delta > {params.pref}.snps
		show-tiling {params.pref}.delta > {params.pref}.tiling
		mv {params.pref}.* mummer/{wildcards.ass_type}
		"""

rule hifiasm_quast:
	input:
		"4_{ass_type}_assembly/hifiasm/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.{gfa_type}.fasta",
	output:
		"quast_gfa/assemblytype_{ass_type}_assemblytool_hifiasmgfa_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}_gfatype_{gfa_type}/report.tsv",
	log:
		"logs/quast/{ass_type}/hifiasmgfa_{read_select}_{sample}_{kmer}_{cov}_{depth}_{gfa_type}.tsv",
	benchmark:
		"benchmarks/quast/assemblytype_{ass_type}_readselect_{read_select}_assemblytool_hifiasmgfa_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}_gfatype_{gfa_type}.tsv",
	threads:
		1
	params:
		out = "quast_gfa/assemblytype_{ass_type}_assemblytool_hifiasmgfa_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}_gfatype_{gfa_type}"
	conda:
		"envs/quast.yaml",
	shell:
		"""
		quast {input} --threads {threads} -o {params.out}
		"""

#Check quality statistics
rule quast:
	input:
		"4_{ass_type}_assembly/{tool}/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	output:
		"quast/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}/report.tsv",
	log:
		"logs/quast/{ass_type}/{tool}_{read_select}_{sample}_{kmer}_{cov}_{depth}.tsv",
	benchmark:
		"benchmarks/quast/assemblytype_{ass_type}_readselect_{read_select}_assemblytool_{tool}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		1
	params:
		out = "quast/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}"
	conda:
		"envs/quast.yaml",
	shell:
		"""
		quast {input} --threads {threads} -o {params.out}
		"""

#rule split_assemblies:
#	input:
#		
#	output:
#
#	log:
#
#	benchmark:
#
#	threads:
#
#	shell:
#		"""
#
#		"""

rule mummer:
	input:
		ref = "reference/{ass_type}.fasta",
		query = "4_{ass_type}_assembly/{tool}/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	output:
		"mummer/{ass_type}/prefix_{sample}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.mums", 
	log:
		"logs/mummer/{ass_type}_{sample}_{tool}_{read_select}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"logs/mummer/assemblytype_{ass_type}_prefix_{sample}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		1
	params:
		pref = "prefix_{sample}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}",
	shadow:
		"shallow"
	conda:
		"envs/mummer.yaml",
	shell:
		"""
		mummer -mum -b -c {input.ref} {input.query} > {params.pref}.mums
		mummerplot --postscript --prefix={params.pref} {params.pref}.mums
		nucmer -maxmatch -c 100 -p {params.pref} {input.ref} {input.query}
		show-coords -r -c -H -d -o -T -l {params.pref}.delta > {params.pref}.coords
		show-snps -C -l -r -T -H {params.pref}.delta > {params.pref}.snps
		show-tiling {params.pref}.delta > {params.pref}.tiling
		mv {params.pref}.* mummer/{wildcards.ass_type}
		"""

