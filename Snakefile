#!/bin/bash

MAX_THREADS = 48
SAMPLES = ["SRR99694"]
PREFIXES = ["SRR9969479.m64011_190430_000122", "SRR9969480.m54119U_190503_125238"]
BBDUK_KMER = ["31"]
BBDUK_COVERAGE = ["0.9"]
CHUNKS = 20
READ_DEPTH = ["10","15","20","25","30","35","40","45","49"]
#10,25,15,35,50
#"10","15","20","25","30","35","40","45","49"
#10,25,75,100,50
LENGTH_CHLOROPLAST = ["134502"]
LENGTH_MITOCHONDRIA = ["491515"]
LENGTH_GENOME = ["387500000"]
ASSEMBLY_TYPE = ["genome"]
ASSEMBLY_TOOLS = ["hifiasm","flye"]
READ_SELECTION_METHOD = ["random","longest"]
HIFIASM_PATH = "/shared/cmatthews/tools/hifiasm/hifiasm"
HICANU_PATH = "/shared/cmatthews/tools/canu/Linux-amd64/bin/canu"
GFA_TYPES = ["p_ctg"]
LTR_FINDER_PATH = "/shared/cmatthews/tools/LTR_FINDER_parallel-master/LTR_FINDER_parallel"
LTRFILES = ["rawLTR.scn","assembly.fasta.out.LAI"]
MUMMER_PATH = "/shared/cmatthews/tools/mummer-4.0.0beta2/" 
CANU_BENCHMARK_FILES = ["job_finish_state","submit_time","start_time","end_time","walltime_reserved","walltime_elapsed","max_memory","max_disk_write","max_disk_read","num_cores"]


localrules:
	all,
	canu,
	hicanu,
	generate_coverage_list,
	select_longest,

rule all:
	input:
#		expand("0_raw/{PREFIX}.subreads.bam.pbi", PREFIX = PREFIXES),
#		expand("1_subreads/chunks/{PREFIX}.ccs.{NUM}.bam", NUM = range(1,CHUNKS+1), PREFIX = PREFIXES),
#		expand("1_subreads/{PREFIX}.ccs.bam", PREFIX = PREFIXES),
#		expand("1_subreads/{SAMPLE}.ccs.bam", SAMPLE = SAMPLES),
#		expand("2_{ASS_TYPE}_reads/unsorted/{SAMPLE}_{KMER}_{COV}.fasta", ASS_TYPE = ASSEMBLY_TYPE, SAMPLE = SAMPLES, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE),
		expand("2_{ASS_TYPE}_reads/sorted/{SAMPLE}_{KMER}_{COV}.fasta", ASS_TYPE = ASSEMBLY_TYPE, SAMPLE = SAMPLES, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE),
		expand("2_{ASS_TYPE}_reads/sorted/{SAMPLE}_{KMER}_{COV}_coveragetable.txt", ASS_TYPE = ASSEMBLY_TYPE, SAMPLE = SAMPLES, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE),	
		expand("3_{ASS_TYPE}_subset/{READ_SELECT}/{SAMPLE}_{KMER}_{COV}_{DEPTH}.fasta", ASS_TYPE = ASSEMBLY_TYPE, READ_SELECT = READ_SELECTION_METHOD, SAMPLE = SAMPLES, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE, DEPTH = READ_DEPTH),
		expand("4_{ASS_TYPE}_assembly/{TOOL}/{SAMPLE}/{READ_SELECT}_{KMER}_{COV}_{DEPTH}/assembly.fasta", ASS_TYPE = ASSEMBLY_TYPE, TOOL = ASSEMBLY_TOOLS, SAMPLE = SAMPLES, READ_SELECT = READ_SELECTION_METHOD, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE, DEPTH = READ_DEPTH),
		expand("quast/assemblytype_{ASS_TYPE}_assemblytool_{TOOL}_readselect_{READ_SELECT}_prefix_{SAMPLE}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}/report.tsv", ASS_TYPE = ASSEMBLY_TYPE, TOOL = ASSEMBLY_TOOLS, READ_SELECT = READ_SELECTION_METHOD, SAMPLE = SAMPLES, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE, DEPTH = READ_DEPTH),
		expand("quastref/assemblytype_{ASS_TYPE}_assemblytool_{TOOL}_readselect_{READ_SELECT}_prefix_{SAMPLE}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}/report.tsv", ASS_TYPE = ASSEMBLY_TYPE, TOOL = ASSEMBLY_TOOLS, READ_SELECT = READ_SELECTION_METHOD, SAMPLE = SAMPLES, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE, DEPTH = READ_DEPTH),
		expand("mummer/{ASS_TYPE}/prefix_{SAMPLE}_assemblytool_{TOOL}_readselect_{READ_SELECT}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}.delta", ASS_TYPE = ASSEMBLY_TYPE, SAMPLE = SAMPLES, TOOL = ASSEMBLY_TOOLS, READ_SELECT = READ_SELECTION_METHOD, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE, DEPTH = READ_DEPTH),
		expand("busco/assemblytype_genome_assemblytool_{TOOL}_readselect_{READ_SELECT}_prefix_{SAMPLE}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}/results/run_poales_odb10/full_table.tsv", TOOL = ASSEMBLY_TOOLS, READ_SELECT = READ_SELECTION_METHOD, SAMPLE = SAMPLES, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE, DEPTH = READ_DEPTH),
		expand("ltr/harvest/assemblytype_genome_assemblytool_{TOOL}_readselect_{READ_SELECT}_prefix_{SAMPLE}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}/assembly.fa.harvest.scn",TOOL = ASSEMBLY_TOOLS, READ_SELECT = READ_SELECTION_METHOD, SAMPLE = SAMPLES, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE, DEPTH = READ_DEPTH),
		expand("ltr/finder/assemblytype_genome_assemblytool_{TOOL}_readselect_{READ_SELECT}_prefix_{SAMPLE}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}/assembly.fasta.finder.combine.scn", TOOL = ASSEMBLY_TOOLS, READ_SELECT = READ_SELECTION_METHOD, SAMPLE = SAMPLES, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE, DEPTH = READ_DEPTH),
		expand("ltr/retriever/assemblytype_genome_assemblytool_{TOOL}_readselect_{READ_SELECT}_prefix_{SAMPLE}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}/complete", TOOL = ASSEMBLY_TOOLS, READ_SELECT = READ_SELECTION_METHOD, SAMPLE = SAMPLES, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE, DEPTH = READ_DEPTH), 
		expand("benchcanu/assemblytype_{ASS_TYPE}_assemblytool_{TOOL}_prefix_{SAMPLE}_readselect_{READ_SELECT}_kmer_{KMER}_cov_{COV}_depth_{DEPTH}/{CANU_BENCHMARKS}.txt",ASS_TYPE = ASSEMBLY_TYPE, TOOL = ["canu","hicanu"], READ_SELECT = READ_SELECTION_METHOD, SAMPLE = SAMPLES, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE, DEPTH = READ_DEPTH, CANU_BENCHMARKS = CANU_BENCHMARK_FILES),
		expand("length/{ASS_TYPE}_{TOOL}_{SAMPLE}_{READ_SELECT}_{KMER}_{COV}_{DEPTH}.txt", ASS_TYPE = ASSEMBLY_TYPE, SAMPLE = SAMPLES, TOOL = ASSEMBLY_TOOLS, READ_SELECT = READ_SELECTION_METHOD, KMER = BBDUK_KMER, COV = BBDUK_COVERAGE, DEPTH = READ_DEPTH),


#creates an index file of the raw hifi reads
#rule pbindex:
#	input:
#		"0_raw/{prefix, SRR[0-9]*\.\m[0-9]*\.[0-9]*}.subreads.bam",	
#	output:
#		"0_raw/{prefix}.subreads.bam.pbi",
#	log:
#		"logs/pbindex/{prefix}.log",
#	benchmark:
#		"benchmarks/pbindex/prefix_{prefix}.tsv",
#	threads:
#		1	
#	conda:
#		"envs/pacbio.yaml"
#	shell:
#		"""
#		pbindex {input}
#		"""
#
## generate high quality consensus reads
#rule ccs:
#	input:
#		sub = "0_raw/{prefix, SRR[0-9]*\.\m[0-9]*\.[0-9]*}.subreads.bam",
#		pbi = "0_raw/{prefix}.subreads.bam.pbi",
#	output:
#		"1_subreads/chunks/{prefix}.ccs.{num}.bam",
#	log:
#		"logs/ccs/{prefix}_{num}.log"
#	benchmark:
#		"benchmarks/ccs/prefix_{prefix}_chunknum_{num}.tsv"	
#	threads:
#		MAX_THREADS
#	params:
#		chnk = CHUNKS
#	conda:
#		"envs/pacbio.yaml",
#	shell:
#		"""
#		ccs {input.sub} {output} --chunk {wildcards.num}/{params.chnk} -j {threads}
#		"""
#
##merge chunks
#rule pbmerge:
#	input:
#		expand("1_subreads/chunks/{{prefix}}.ccs.{num}.bam", prefix = PREFIXES, num = range(1,CHUNKS+1)),
#	output:
#		"1_subreads/{prefix}.ccs.bam",
#	log:
#		"logs/pbmerge/{prefix}.log",
#	benchmark:
#		"benchmarks/pbmerge/prefix_{prefix}.log",
#	threads:
#		MAX_THREADS
#	conda:
#		"envs/pacbio.yaml",
#	shell:
#		"""
#		pbmerge -o {output} {input}
#		"""
#
##merge bam files
#rule pbmerge_sample:
#	input:
#		expand("1_subreads/{prefix}.ccs.bam", prefix = PREFIXES),
#	output:
#		"1_subreads/{sample}.ccs.bam",
#	log:
#		"logs/pbmergesample/{sample}.log",
#	benchmark:
#		"benchmarks/pbmergesample/prefix_{sample}.tsv",
#	threads:
#		1
#	conda:
#		"envs/pacbio.yaml",
#	shell:
#		"""
#		pbmerge -o {output} {input}
#		"""	
#
#rule bbduk:
#	input:
#		bam = "1_subreads/{sample}.ccs.bam",
#		ref = "reference/{ass_type}.fasta",
#	output:
#		out = "2_{ass_type}_reads/unsorted/{sample}_{kmer}_{cov}.fasta",
#	log:
#		"logs/bbduk/{sample}_{ass_type}_{kmer}_{cov}.log",
#	benchmark:
#		"benchmarks/bbduk/assemblytype_{ass_type}_prefix_{sample}_kmer_{kmer}_cov_{cov}.tsv",
#	threads:
#		5
#	wildcard_constraints:
#		ass_type = "(mitochondria|chloroplast)",
#	conda:
#		"envs/bbmap.yaml",
#	shell:
#		"""
#		echo "Extracting {wildcards.ass_type} reads with bbduk..."
#		bbduk.sh in={input.bam} outm={output.out} ref={input.ref} threads={threads} k={wildcards.kmer} -Xmx1g mincovfraction={wildcards.cov}
#		"""
#
#rule bbduk_genome:
#	input:
#		seq = "1_subreads/{sample}.ccs.bam",
#		ref1 = "reference/chloroplast.fasta",
#		ref2 = "reference/mitochondria.fasta",
#	output:
#		"2_genome_reads/unsorted/{sample}_{kmer}_{cov}.fasta",
#	log:
#		"logs/bbdukgenome/assemblytype_genome_prefix_{sample}_kmer_{kmer}_cov_{cov}.log",
#	benchmark:
#		"benchmarks/bbdukgenome/assemblytype_genome_prefix_{sample}_kmer_{kmer}_cov_{cov}.tsv",
#	threads:
#		MAX_THREADS
#	conda:
#		"envs/bbmap.yaml",
#	shell:
#		"""
#		bbduk.sh in={input.seq} out={output} ref={input.ref1},{input.ref2} threads={threads} k={wildcards.kmer} -Xmx60g mincovfraction={wildcards.cov}
#		"""
#
##subset the matching reads down to specified coverage
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
	resources:
                time = lambda wildcards, input: (60 if wildcards.ass_type  == 'genome' else 1),
                mem_mb = lambda wildcards, input: (15000 if wildcards.ass_type == 'genome' else 200),
                cpu = lambda wildcards, input: (1 if wildcards.ass_type == 'genome' else 1),
	params:
		chlor = LENGTH_CHLOROPLAST,
		mito = LENGTH_MITOCHONDRIA,
		geno = LENGTH_GENOME,
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
		if [ {wildcards.ass_type} == "genome" ]; then
			LENGTH={params.geno}
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
#
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
		sortbyname.sh in={input} out={output} name=f length=t ascending=f -Xmx60g
		"""

rule generate_coverage_list:
	input:
		"2_{ass_type}_reads/sorted/{sample}_{kmer}_{cov}.fasta",
	output:
		"2_{ass_type}_reads/sorted/{sample}_{kmer}_{cov}_coveragetable.txt",
	log:
		"logs/generatecoveragelist/assemblytype_{ass_type}_prefix_{sample}_kmer_{kmer}_cov_{cov}.log",
	benchmark:
		"benchmarks/generatecoveragelist/assemblytype_{ass_type}_prefix_{sample}_kmer_{kmer}_cov_{cov}.tsv",
	threads:
		1
	params:
		length_mito = LENGTH_MITOCHONDRIA,
		length_chloro = LENGTH_CHLOROPLAST,
		length_genome = LENGTH_GENOME,
	shell:
		"""
		if [ {wildcards.ass_type} == "chloroplast" ]; then
			LEN={params.length_chloro}
		elif [ {wildcards.ass_type} == "mitochondria" ]; then
			LEN={params.length_mito}
		elif [ {wildcards.ass_type} == "genome" ]; then
			LEN={params.length_genome}
		fi
		touch length.tmp && rm length.tmp
		x=0
		awk '/^>/{{if (l!="") print l; l=0; next}}{{l+=length($0)}}' {input} >> length.tmp

		while IFS= read -r line
		do
			x=$(( ${{x}} + ${{line}} ))
			pr=$(( ${{x}} / ${{LEN}} ))
			echo "${{pr}}" >> {output}
		done < "length.tmp"

		rm length.tmp
                """

rule select_longest:
	input:
		fa = "2_{ass_type}_reads/sorted/{sample}_{kmer}_{cov}.fasta",
		list = "2_{ass_type}_reads/sorted/{sample}_{kmer}_{cov}_coveragetable.txt",
	output:
		"3_{ass_type}_subset/longest/{sample}_{kmer}_{cov}_{depth}.fasta",
	log:
		"logs/select_longest/{ass_type}/{sample}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/selectlongest/assemblytype_{ass_type}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		1
	params:
		length_mito = LENGTH_MITOCHONDRIA,
		length_chloro = LENGTH_CHLOROPLAST,
		length_genome = LENGTH_GENOME,
	shell:
		"""
		if [ {wildcards.ass_type} == "chloroplast" ]; then
			LEN={params.length_chloro}
		elif [ {wildcards.ass_type} == "mitochondria" ]; then
			LEN={params.length_mito}
		elif [ {wildcards.ass_type} == "genome" ]; then
			LEN={params.length_genome}
		fi

		NO_READS=$(grep '^>' {input.fa} | wc -l)
		NUM=$(awk '$1<{wildcards.depth}{{c++}} END{{print c+0}}' {input.list})
		NUM=$(( ${{NUM}} + 1 ))

		if [ ${{NUM}} -eq ${{NO_READS}} ]; then
			DEPTH=$( tail -n 1 {input.list})
		echo "Number of reads required for requested coverage is greater than the number of reads available."
			echo "All reads will be used, giving a coverage depth of ${{DEPTH}}x"
		fi

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
		MAX_THREADS
	params:
		chlor = LENGTH_CHLOROPLAST,
		mito = LENGTH_MITOCHONDRIA,
		genome = LENGTH_GENOME,
		out = "4_{ass_type}_assembly/flye/{sample}/{read_select}_{kmer}_{cov}_{depth}",
	resources:
		time = lambda wildcards, input: (2000 if wildcards.ass_type  == 'genome' else 10),
		mem_mb = lambda wildcards, input: (370000 if wildcards.ass_type == 'genome' else 5000),
		cpu = lambda wildcards, input: (MAX_THREADS if wildcards.ass_type == 'genome' else 5), 
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

		if [ {wildcards.ass_type} == "genome" ]; then
			SIZE={params.genome}
		fi

		flye --pacbio-hifi {input.file} --genome-size ${{SIZE}} --out-dir {params.out} --threads {threads}
		"""

#No spaces allowed in canu command between option and value
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
			ETIME="--time=00:10:00"
			echo "Starting chloroplast assembly ..."
		fi
		if [ {wildcards.ass_type} == 'mitochondria' ]; then
			SIZE={params.mito}
			JTIME="--time=00:10:00"
			ETIME="--time=00:10:00"
			echo "Starting mitochondria assembly ..."
		fi	
		if [ {wildcards.ass_type} == 'genome' ]; then
			SIZE={params.genome}
			JTIME="--time=72:00:00"
			ETIME="--time=72:00:00"
			echo "Starting genome assembly ..."
		fi
		
		(canu -p {params.prefix} -d {params.dir} genomeSize=${{SIZE}} gridOptions=${{JTIME}} gridOptionsJobName={params.jobName} gridOptionsExecutive=${{ETIME}} executiveMemory=1 correctedErrorRate=0.015 ovlMerThreshold=75 batOptions="-eg 0.01 -eM 0.01 -dg 6 -db 6 -dr 1 -ca 50 -cp 5" -pacbio-corrected {input} stopOnLowCoverage=3 onSuccess=touch onFailure=./canuFailure.sh) 2> {log}

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
	conda:
		"envs/canuv2.yaml",
	shell:
		"""	
		cp canuFailure.sh {params.dir}

		ETIME="--time=00:10:00"

		if [ {wildcards.ass_type} == 'chloroplast' ]; then
			SIZE={params.chlor}
			JTIME="--time=00:40:00"
			echo "Starting HiCanu chloroplast assembly ..."
		fi
		if [ {wildcards.ass_type} == 'mitochondria' ]; then
			SIZE={params.mito}
			JTIME="--time=00:10:00"
			echo "Starting HiCanu mitochondria assembly ..."
		fi	
		if [ {wildcards.ass_type} == 'genome' ]; then
			SIZE={params.genome}
			JTIME="--time=72:00:00"
			ETIME="--time=72:00:00"
			echo "Starting HiCanu genome assembly ..."
		fi
		
		({params.hicanu} -p {params.prefix} -d {params.dir} genomeSize=${{SIZE}} gridOptions=${{JTIME}} gridOptionsJobName={params.jobName} gridOptionsExecutive=${{ETIME}} executiveMemory=1 -pacbio-hifi {input} minInputCoverage=9 stopOnLowCoverage=3 onSuccess=touch onFailure=./canuFailure.sh) 2> {log}

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

rule benchcanu:
	input:
		"4_{ass_type}_assembly/{tool}/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta"
	output:
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/job_finish_state.txt",
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/submit_time.txt",
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/start_time.txt",
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/end_time.txt",
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/walltime_reserved.txt",
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/walltime_elapsed.txt",
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/max_memory.txt",
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/max_disk_write.txt",
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/max_disk_read.txt",
		"benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/num_cores.txt",
	log:
		"logs/benchcanu/{ass_type}/{tool}/{sample}_{read_select}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		1
	params:
		dir = "4_{ass_type}_assembly/{tool}/{sample}/{read_select}_{kmer}_{cov}_{depth}/",
		dir2 = "benchcanu/assemblytype_{ass_type}_assemblytool_{tool}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}/",
	shell:
		"""
		grep -r 'State               :' {params.dir} > {params.dir2}job_finish_state.txt
		grep -r 'Submit              : 2020' {params.dir} > {params.dir2}submit_time.txt
		grep -r 'Start               : 2020' {params.dir} > {params.dir2}start_time.txt
		grep -r 'End                 : 2020' {params.dir} > {params.dir2}end_time.txt
		grep -r 'Walltime reserved   :' {params.dir} > {params.dir2}walltime_reserved.txt
		grep -r 'Walltime elapsed (%):' {params.dir} > {params.dir2}walltime_elapsed.txt
		grep -r '% Mem used (Max)    :' {params.dir} > {params.dir2}max_memory.txt
		grep -r 'Max Disk Write      :' {params.dir} > {params.dir2}max_disk_write.txt
		grep -r 'Max Disk Read       :' {params.dir} > {params.dir2}max_disk_read.txt
		grep -r 'Cores               :' {params.dir} > {params.dir2}num_cores.txt
		"""

rule hifiasm:
	input:
		"3_{ass_type}_subset/{read_select}/{sample}_{kmer}_{cov}_{depth}.fasta",
	output:
		main = "4_{ass_type}_assembly/hifiasm/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
#		actg = "4_{ass_type}_assembly/hifiasm/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.a_ctg.gfa",
#		pctg = "4_{ass_type}_assembly/hifiasm/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.p_ctg.gfa",
#		putg = "4_{ass_type}_assembly/hifiasm/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.p_utg.gfa",
#		rutg = "4_{ass_type}_assembly/hifiasm/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.r_utg.gfa",
	log:
		"logs/hifiasm/{ass_type}_{sample}_{read_select}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/hifiasm/assemblytype_{ass_type}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	params:	
		prefix = "assembly",
		hifiasm = HIFIASM_PATH,
		dir = "4_{ass_type}_assembly/hifiasm/{sample}/{read_select}_{kmer}_{cov}_{depth}",
	shadow:
		"shallow"
	resources:
		time = lambda wildcards, input: (2000 if wildcards.ass_type  == 'genome' else 10),
		mem_mb = lambda wildcards, input: (370000 if wildcards.ass_type == 'genome' else 5000),
		cpu = lambda wildcards, input: (MAX_THREADS if wildcards.ass_type == 'genome' else 5),
	threads:
		MAX_THREADS
	shell:
		"""
		{params.hifiasm} -o {params.prefix} -t {threads} {input}
		mv {params.prefix}* {params.dir}/
		awk '/^S/{{print ">" $2 "\\n" $3}}' {params.dir}/{params.prefix}.p_ctg.gfa | fold > {output.main}
		"""

#rule gfa2fa:
#	input:
#		"4_{ass_type}_assembly/hifiasm/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.{gfa_type}.gfa",
#	output:
#		"4_{ass_type}_assembly/hifiasm/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.{gfa_type}.fasta",
#	log:
#		"logs/gfa2fa/{ass_type}_{sample}_{read_select}_{kmer}_{cov}_{depth}_{gfa_type}.log",
#	benchmark:
#		"benchmarks/gfa2fa/assemblytype_{ass_type}_prefix_{sample}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}_gfatype_{gfa_type}.tsv"
#	threads:
#		1
#	shell:
#		"""
#		awk '/^S/{{print ">" $2 "\\n" $3}}' {input} | fold > {output}
#		"""

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
	resources:
		time = lambda wildcards, input: (5 if wildcards.ass_type  == 'genome' else 1),
	params:
		out = "quast/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}"
	singularity:
		"docker://quay.io/biocontainers/quast:5.0.2--1"
	shell:
		"""
		quast {input} --threads {threads} -o {params.out}
		"""


rule quastref:
	input:
		ref = "reference/{ass_type}.fasta",
		ass = "4_{ass_type}_assembly/{tool}/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	output:
		"quastref/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}/report.tsv",
	log:
		"logs/quastref/{ass_type}/{tool}_{read_select}_{sample}_{kmer}_{cov}_{depth}.tsv",
	benchmark:
		"benchmarks/quastref/assemblytype_{ass_type}_readselect_{read_select}_assemblytool_{tool}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		MAX_THREADS
	resources:
		mem_mb = lambda wildcards, input: (15000 if wildcards.ass_type == 'genome' else 3000),
		cpu = lambda wildcards, input: (5 if wildcards.ass_type == 'genome' else 1),
		time = lambda wildcards, input: (60 if wildcards.ass_type == 'genome' else 2),
	params:
		out = "quastref/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}"
	singularity:
		"docker://quay.io/biocontainers/quast:5.0.2--1"
	shell:
		"""
		if [ {wildcards.ass_type} == 'genome' ]; then
			quast {input.ass} --eukaryote --no-icarus --no-html --large -r {input.ref} --threads {threads} -o {params.out}
		else
			quast {input.ass} --no-icarus --no-html -r {input.ref} --threads {threads} -o {params.out}
		fi
		"""

rule busco:
	input:
		fasta = "4_genome_assembly/{tool}/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	output:
		"busco/assemblytype_genome_assemblytool_{tool}_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}/results/run_poales_odb10/full_table.tsv",
	log:
		"logs/busco/genome/{tool}_{read_select}_{sample}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/busco/assemblytype_genome_readselect_{read_select}_assemblytool_{tool}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		MAX_THREADS
	shadow:
		"full",
	singularity:
		"docker://ezlabgva/busco:v4.0.5_cv1",
	params:
		outdir = "assemblytype_genome_assemblytool_{tool}_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}",
		lineage_path = "./reference/poales_odb10",
	shell:
		"""
		busco -f --in {input.fasta} --out results --lineage_dataset {params.lineage_path} --cpu {threads} --mode genome
		mv results/short_summary.specific.poales_odb10.results.txt busco/{params.outdir}/results
		mv results/blast_db busco/{params.outdir}/results
		mv results/logs busco/{params.outdir}/results
		mv results/run_poales_odb10/* busco/{params.outdir}/results/run_poales_odb10/
		"""

rule gt_ltrharvest:
	input:
		"4_{ass_type}_assembly/{tool}/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	output:
		"ltr/harvest/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}/assembly.fa.harvest.scn",
	log:
		"logs/gt_ltrharvest/{ass_type}/{tool}_{read_select}_{sample}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/gtltrharvest/assemblytype_{ass_type}_readselect_{read_select}_assemblytool_{tool}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		10
	shadow:
		"shallow",
	wildcard_constraints:
		ass_type = "genome",
	conda:
		"envs/genometools.yaml",
	params:
		dir = "ltr/harvest/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}/",
	shell:
		"""
		gt suffixerator -db {input} -indexname suffixer.fa -tis -suf -lcp -des -ssp -sds -dna
		gt ltrharvest -index suffixer.fa -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > out.scn
		mv  out.scn {output}
		mv suffixer* {params.dir}
		"""

rule ltr_finder:
	input:
		"4_{ass_type}_assembly/{tool}/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	output:
		"ltr/finder/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}/assembly.fasta.finder.combine.scn",
	log:
		"logs/ltr_finder/{ass_type}/{tool}_{read_select}_{sample}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/ltrfinder/assemblytype_{ass_type}_readselect_{read_select}_assemblytool_{tool}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		10
	shadow:
		"shallow",
	wildcard_constraints:
		ass_type = "genome",
	params:
		dir = "ltr/finder/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}/",
		finder = LTR_FINDER_PATH,
	shell:
		"""
		perl {params.finder} -seq {input} -threads {threads} -harvest_out -size 1000000 -time 300
		mv assembly* {params.dir}
		"""

rule ltr_retriever:
	input:
		genome = "4_{ass_type}_assembly/{tool}/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
		finder = "ltr/finder/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}/assembly.fasta.finder.combine.scn",
		harvest = "ltr/harvest/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}/assembly.fa.harvest.scn",
	output:
		scn = "ltr/retriever/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}/rawLTR.scn",
		lai =  "ltr/retriever/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}/assembly.fasta.out.LAI",
		dummy = "ltr/retriever/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}/complete",
	log:
		"logs/ltr_retriever/{ass_type}/{tool}_{read_select}_{sample}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"benchmarks/ltrretriever/assemblytype_{ass_type}_readselect_{read_select}_assemblytool_{tool}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		10
	shadow:
		"shallow",
	wildcard_constraints:
		ass_type = "genome",
	conda:
		"envs/ltr.yaml",
	params:
		dir = "ltr/retriever/assemblytype_{ass_type}_assemblytool_{tool}_readselect_{read_select}_prefix_{sample}_kmer_{kmer}_cov_{cov}_depth_{depth}/",
	shell:
		"""
		cat {input.finder} {input.harvest} > {output.scn}
		cp {input.genome} {params.dir}
		LTR_retriever -genome {params.dir}assembly.fasta -inharvest {output.scn} -threads {threads}
		
		if [ {wildcards.tool} == 'hifiasm' ]; then
			mv assembly.fasta.out.* {params.dir}
			rm assembly.fasta.out
			rm {params.dir}assembly.fasta
			touch {output.dummy}
		fi

		if [ {wildcards.tool} == 'canu' ]; then
			mv assembly.fasta.mod.out.* {params.dir}
			rm assembly.fasta.mod.out
			rm {params.dir}assembly.fasta
			mv {params.dir}assembly.fasta.mod.out.LAI {params.dir}assembly.fasta.out.LAI
			touch {output.dummy}
		fi

		if [ {wildcards.tool} == 'hicanu' ]; then
			mv assembly.fasta.mod.out.* {params.dir}
			rm assembly.fasta.mod.out
			rm {params.dir}assembly.fasta
			mv {params.dir}assembly.fasta.mod.out.LAI {params.dir}assembly.fasta.out.LAI
			touch {output.dummy}
		fi

		"""

rule contig_length:
	input:
		"4_{ass_type}_assembly/{tool}/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	output:
		"length/{ass_type}_{tool}_{sample}_{read_select}_{kmer}_{cov}_{depth}.txt",	
	log:
		"logs/length/{ass_type}_{sample}_{tool}_{read_select}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"logs/length/assemblytype_{ass_type}_prefix_{sample}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	shell:
		"""
		cat {input} | awk '$0 ~ ">" {{if (NR > 1) {{print c;}} c=0;printf substr($0,2,100) "\t"; }} $0 !~ ">" {{c+=length($0);}} END {{ print c; }}' > {output}
		"""
rule mummer:
	input:
		ref = "reference/{ass_type}.fasta",
		query = "4_{ass_type}_assembly/{tool}/{sample}/{read_select}_{kmer}_{cov}_{depth}/assembly.fasta",
	output:
		"mummer/{ass_type}/prefix_{sample}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.delta", 
	log:
		"logs/mummer/{ass_type}_{sample}_{tool}_{read_select}_{kmer}_{cov}_{depth}.log",
	benchmark:
		"logs/mummer/assemblytype_{ass_type}_prefix_{sample}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}.tsv",
	threads:
		MAX_THREADS
	resources:
		time = lambda wildcards, input: (90 if wildcards.ass_type == "genome" else 1),
		mem_mb = lambda wildcards, input: (40000 if wildcards.ass_type == "genome" else 200),
		cpu = lambda wildcards, input: (20 if wildcards.ass_type == "genome" else 1),
	params:
		pref = "prefix_{sample}_assemblytool_{tool}_readselect_{read_select}_kmer_{kmer}_cov_{cov}_depth_{depth}",
		mummer = MUMMER_PATH,
	shadow:
		"shallow",
	shell:
		"""
		({params.mummer}nucmer -t {threads} --prefix={params.pref} {input.ref} {input.query}) 2> {log}
		{params.mummer}show-coords -r -c -H -d -o -T -l {params.pref}.delta > {params.pref}.coords
		{params.mummer}show-snps -C -l -r -T -H {params.pref}.delta > {params.pref}.snps
		{params.mummer}show-tiling {params.pref}.delta > {params.pref}.tiling
		{params.mummer}mummerplot --postscript --prefix={params.pref} {params.pref}.delta
		mv {params.pref}.* mummer/{wildcards.ass_type}
		"""



