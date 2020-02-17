MAX_THREADS = 32

PREFIXES = ["SRR9969479.m64011_190430_000122", "SRR9969480.m54119U_190503_125238"]

CHUNKS = 20

localrules:
	all

rule all:
	input:
		expand("0_raw/{PREFIX}.subreads.bam.pbi", PREFIX=PREFIXES),
		expand("1_subreads/chunks/{PREFIX}.ccs.{NUM}.bam", NUM = range(1,CHUNKS+1), PREFIX = PREFIXES),
		expand("1_subreads/{PREFIX}.ccs.bam", PREFIX = PREFIXES),

#creates an index file of the raw hifi reads
rule pbindex:
	input:
		"0_raw/{prefix, SRR[0-9]*\.\m[0-9]*\.[0-9]*}.subreads.bam",	
	output:
		"0_raw/{prefix}.subreads.bam.pbi",
	log:
		"logs/pbindex/{prefix}.log",
	benchmark:
		"benchmarks/pbindex/{prefix}.tsv",
	threads:
		1	
	conda:
		"envs/pacbio.yaml"
	shell:
		"""
		(pbindex {input}) 2> {log}
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
		"benchmarks/ccs/{prefix}_{num}.tsv"	
	threads:
		MAX_THREADS
	params:
		chnk = CHUNKS
	conda:
		"envs/pacbio.yaml",
	shell:
		"""
		(ccs {input.sub} {output} --chunk {wildcards.num}/{params.chnk} -j {threads}) 2> {log}
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
		"benchmarks/pbmerge/{prefix}.log",
	threads:
		MAX_THREADS
	conda:
		"envs/pacbio.yaml",
	shell:
		"""
		pbmerge -o {output} {input}
		"""

#rule ccs:
#        input:
#                "0_raw/{prefix, SRR[0-9]*\.\m[0-9]*\.[0-9]*}.subreads.bam",
#        output:
#                "1_subreads/{prefix}.ccs.{num}.bam",
#        log:
#                "logs/ccs/{prefix}_{num}.log"
#        benchmark:
#                "benchmarks/ccs/{prefix}_{num}.tsv"
#        threads:
#                MAX_THREADS
#        conda:
#                "envs/pacbio.yaml",
#        shell:
#                """
#                ccs {input} {output} --chunk {wildcards.num}/CHUNKS -j {threads}
#                """

