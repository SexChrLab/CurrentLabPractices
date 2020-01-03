import os

configfile: "mapping_to_GRCh38_config.json"

# Tool paths:
fastqc_path = "fastqc"
multiqc_path = "multiqc"
bwa_path = "bwa"
samtools_path = "samtools"
bbduksh_path = "bbduk.sh"
picard_path = "picard"
gatk_path = "/home/tphung3/softwares/gatk-4.1.0.0/gatk"
gatk3_path = "/home/tphung3/softwares/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar"
bcftools_path = "bcftools"

rule all:
    input: #gatk
        "PlacentaSamples/no.altaware.gatk.raw.g.vcf.gz",
        "PlacentaSamples/with.altaware.gatk.raw.g.vcf.gz"
    input: #mapping
        expand("PlacentaSamples/{sample_name}.no.altaware.sorted.bam", sample_name=config["female_placentas"]),
        expand("PlacentaSamples/{sample_name}.with.altaware.sorted.bam", sample_name=config["female_placentas"]),
        expand("PlacentaSamples/{sample_name}.no.altaware.sorted.mkdup.bam.bai", sample_name=config["female_placentas"]),
        expand("PlacentaSamples/{sample_name}.with.altaware.sorted.mkdup.bam.bai", sample_name=config["female_placentas"])
    input: #prepare reference
        "reference/GRCh38_full_analysis_set_plus_decoy_hla.XXonly.fa",
        "reference/GRCh38_full_analysis_set_plus_decoy_hla.XY.fa"

rule xyalign_create_references:
    input:
        ref = "/data/storage/SAYRES/REFERENCE_GENOMES/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        mask = "/agavescratch/tphung3/SmallProjects/MappingToGRCh38/reference/mask.bed"
    output:
        xx = "reference/GRCh38_full_analysis_set_plus_decoy_hla.XXonly.fa",
        xy = "reference/GRCh38_full_analysis_set_plus_decoy_hla.XY.fa"
    conda:
        "envs/xyalign.yml"
    shell:
        "xyalign --PREPARE_REFERENCE --ref {input.ref} "
        "--xx_ref_out {output.xx} --xy_ref_out {output.xy} "
        "--x_chromosome chrX "
        "--y_chromosome chrY "
        "--reference_mask {input.mask} "
        "--output_dir xyalign "

rule prepare_reference_female:
	input:
		"reference/GRCh38_full_analysis_set_plus_decoy_hla.XXonly.fa"
	output:
		fai = "reference/GRCh38_full_analysis_set_plus_decoy_hla.XXonly.fa.fai",
		amb = "reference/GRCh38_full_analysis_set_plus_decoy_hla.XXonly.fa.amb",
		dict = "reference/GRCh38_full_analysis_set_plus_decoy_hla.XXonly.fa.dict"
	run:
		# faidx
		shell("samtools faidx {input}")
		# .dict
		shell("samtools dict -o {output.dict} {input}")
		# bwa
		shell("bwa index {input}")

rule bwa_no_altaware:
    input:
        fq1 = os.path.join(config["trimmed_fq_path"], "{sample_name}_trimmed_R1.fastq.gz"),
        fq2 = os.path.join(config["trimmed_fq_path"], "{sample_name}_trimmed_R2.fastq.gz"),
        ref = "reference/GRCh38_full_analysis_set_plus_decoy_hla.XXonly.fa",
        fai = "reference/GRCh38_full_analysis_set_plus_decoy_hla.XXonly.fa.fai"
    output:
        "PlacentaSamples/{sample_name}.no.altaware.sorted.bam"
    params:
        id = lambda wildcards: config[wildcards.sample_name]["ID"],
        sm = lambda wildcards: config[wildcards.sample_name]["SM"],
        lb = lambda wildcards: config[wildcards.sample_name]["LB"],
        pu = lambda wildcards: config[wildcards.sample_name]["PU"],
        pl = lambda wildcards: config[wildcards.sample_name]["PL"],
        bwa = bwa_path,
        samtools = samtools_path,
        threads = 4
    threads: 4
    priority: 100
    shell:
        " {params.bwa} mem -j -R "
        "'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
        "{input.ref} {input.fq1} {input.fq2}"
        "| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
        "-O bam -o {output}"

rule bwa_with_altaware:
    input:
        fq1 = os.path.join(config["trimmed_fq_path"], "{sample_name}_trimmed_R1.fastq.gz"),
        fq2 = os.path.join(config["trimmed_fq_path"], "{sample_name}_trimmed_R2.fastq.gz"),
        ref = "reference/GRCh38_full_analysis_set_plus_decoy_hla.XXonly.fa",
        fai = "reference/GRCh38_full_analysis_set_plus_decoy_hla.XXonly.fa.fai"
    output:
        "PlacentaSamples/{sample_name}.with.altaware.sorted.bam"
    params:
        id = lambda wildcards: config[wildcards.sample_name]["ID"],
        sm = lambda wildcards: config[wildcards.sample_name]["SM"],
        lb = lambda wildcards: config[wildcards.sample_name]["LB"],
        pu = lambda wildcards: config[wildcards.sample_name]["PU"],
        pl = lambda wildcards: config[wildcards.sample_name]["PL"],
        bwa = bwa_path,
        samtools = samtools_path,
        threads = 4
    threads: 4
    priority: 100
    shell:
        " {params.bwa} mem -R "
        "'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
        "{input.ref} {input.fq1} {input.fq2}"
        "| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
        "-O bam -o {output}"

rule index_bam_no_altaware:
    input:
        "PlacentaSamples/{sample_name}.no.altaware.sorted.bam"
    output:
        "PlacentaSamples/{sample_name}.no.altaware.sorted.bam.bai"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} index {input}"

rule index_bam_with_altaware:
    input:
        "PlacentaSamples/{sample_name}.with.altaware.sorted.bam"
    output:
        "PlacentaSamples/{sample_name}.with.altaware.sorted.bam.bai"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} index {input}"

rule picard_mkdups_no_altaware:
    input:
        bam = "PlacentaSamples/{sample_name}.no.altaware.sorted.bam",
        bai = "PlacentaSamples/{sample_name}.no.altaware.sorted.bam.bai"
    output:
        bam = "PlacentaSamples/{sample_name}.no.altaware.sorted.mkdup.bam",
        metrics = "stats/{sample_name}.no.altaware.sorted.mkdup_metrics.txt"
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.bam} O={output.bam} "
        "M={output.metrics}"

rule picard_mkdups_with_altaware:
    input:
        bam = "PlacentaSamples/{sample_name}.with.altaware.sorted.bam",
        bai = "PlacentaSamples/{sample_name}.with.altaware.sorted.bam.bai"
    output:
        bam = "PlacentaSamples/{sample_name}.with.altaware.sorted.mkdup.bam",
        metrics = "stats/{sample_name}.with.altaware.sorted.mkdup_metrics.txt"
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.bam} O={output.bam} "
        "M={output.metrics}"

rule index_mkdup_bam_no_altaware:
    input:
        "PlacentaSamples/{sample_name}.no.altaware.sorted.mkdup.bam"
    output:
        "PlacentaSamples/{sample_name}.no.altaware.sorted.mkdup.bam.bai"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} index {input}"

rule index_mkdup_bam_with_altaware:
    input:
        "PlacentaSamples/{sample_name}.with.altaware.sorted.mkdup.bam"
    output:
        "PlacentaSamples/{sample_name}.with.altaware.sorted.mkdup.bam.bai"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} index {input}"

# --------------------------------
# Joint genotyoe calling with GATK
# --------------------------------
# picard CreateSequenceDictionary R=reference/GRCh38_full_analysis_set_plus_decoy_hla.XXonly.fa O=reference/GRCh38_full_analysis_set_plus_decoy_hla.XXonly.dict

rule gatk_gvcf_no_altaware:
	input:
		ref = "reference/GRCh38_full_analysis_set_plus_decoy_hla.XXonly.fa",
		bam = "PlacentaSamples/{sample_name}.no.altaware.sorted.mkdup.bam",
		bai = "PlacentaSamples/{sample_name}.no.altaware.sorted.mkdup.bam.bai"
	output:
		"PlacentaSamples/{sample_name}.no.altaware.g.vcf.gz"
	params:
		gatk = gatk_path,
	shell:
		"{params.gatk} "
		"HaplotypeCaller -R {input.ref} -I {input.bam} "
		"--emit-ref-confidence GVCF --output {output}"

rule gatk_combinegvcfs_no_altaware:
    input:
        ref = "reference/GRCh38_full_analysis_set_plus_decoy_hla.XXonly.fa",
        gvcfs = lambda wildcards: expand(
			"PlacentaSamples/{sample_name}.no.altaware.g.vcf.gz",
			sample_name=config["female_placentas"])
    params:
        gatk = gatk_path

    output:
        "PlacentaSamples/no.altaware.gatk.combinegvcf.g.vcf.gz"

    run:
        variant_files = []
        for i in input.gvcfs:
	           variant_files.append("--variant " + i)
        variant_files = " ".join(variant_files)
        print(
	       """{params.gatk} --java-options "-Xmx10g" """
	          """CombineGVCFs -R {input.ref} {variant_files} -O {output}"""
        )
        shell(
	       """{params.gatk} --java-options "-Xmx10g" """
	          """CombineGVCFs -R {input.ref} {variant_files} -O {output}""")

rule gatk_genotypegvcf_joint_no_altaware:
    input:
        ref = "reference/GRCh38_full_analysis_set_plus_decoy_hla.XXonly.fa",
        gvcf = "PlacentaSamples/no.altaware.gatk.combinegvcf.g.vcf.gz"
    output:
        "PlacentaSamples/no.altaware.gatk.raw.g.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx10g" """
        """GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output} """

#
rule gatk_gvcf_with_altaware:
	input:
		ref = "reference/GRCh38_full_analysis_set_plus_decoy_hla.XXonly.fa",
		bam = "PlacentaSamples/{sample_name}.with.altaware.sorted.mkdup.bam",
		bai = "PlacentaSamples/{sample_name}.with.altaware.sorted.mkdup.bam.bai"
	output:
		"PlacentaSamples/{sample_name}.with.altaware.g.vcf.gz"
	params:
		gatk = gatk_path,
	shell:
		"{params.gatk} "
		"HaplotypeCaller -R {input.ref} -I {input.bam} "
		"--emit-ref-confidence GVCF --output {output}"

rule gatk_combinegvcfs_with_altaware:
    input:
        ref = "reference/GRCh38_full_analysis_set_plus_decoy_hla.XXonly.fa",
        gvcfs = lambda wildcards: expand(
			"PlacentaSamples/{sample_name}.with.altaware.g.vcf.gz",
			sample_name=config["female_placentas"])
    params:
        gatk = gatk_path

    output:
        "PlacentaSamples/with.altaware.gatk.combinegvcf.g.vcf.gz"

    run:
        variant_files = []
        for i in input.gvcfs:
	           variant_files.append("--variant " + i)
        variant_files = " ".join(variant_files)
        print(
	       """{params.gatk} --java-options "-Xmx10g" """
	          """CombineGVCFs -R {input.ref} {variant_files} -O {output}"""
        )
        shell(
	       """{params.gatk} --java-options "-Xmx10g" """
	          """CombineGVCFs -R {input.ref} {variant_files} -O {output}""")

rule gatk_genotypegvcf_joint_with_altaware:
    input:
        ref = "reference/GRCh38_full_analysis_set_plus_decoy_hla.XXonly.fa",
        gvcf = "PlacentaSamples/with.altaware.gatk.combinegvcf.g.vcf.gz"
    output:
        "PlacentaSamples/with.altaware.gatk.raw.g.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx10g" """
        """GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output} """
