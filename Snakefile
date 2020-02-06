MAX_THREADS = 32

PREFIXES = ["SRR9969479.m64011_190430_000122", "SRR9969480.m54119U_190503_125238"]

localrules:
	all

rule all:
	input:
		expand("0_raw/{PREFIX}.subreads.bam.pbi", PREFIX=PREFIXES),

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

