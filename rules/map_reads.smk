## Author: Sambit K. Mishra
## Created: 02-21-2021

## Map reads for control and treatment to human genome using STAR
rule map_reads:
    input:
        fw = "{dd}" + "/results/samples_trimmed/" + "{sample_id}" + "_R1_trimmed.fq.gz",
        rev = "{dd}" + "/results/samples_trimmed/" + "{sample_id}" + "_R2_trimmed.fq.gz",
    output:
        #outdir = expand("{dd}" + "/results/mapped_reads_control/human/",dd=config["DATADIR"]),
        samout = "{dd}" + "/results/mapped_reads/human/" + "{sample_id}Aligned.sortedByCoord.out.bam",
    params:
        gd = config["DATADIR"] + "/human_genome_indexed/",
        outdir = config["DATADIR"] + "/results/mapped_reads/human/"
    conda:
        "mapping.yml"
    threads: config["THREADS"]    
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {params.gd} --readFilesIn {input.fw} {input.rev} --outFileNamePrefix {params.outdir}{wildcards.sample_id} --sjdbOverhang 100 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --quantMode GeneCounts
        """

# ## Map reads for control to human genome using STAR
# rule map_reads_trt:
#     input:
#         ctr_fw = expand("{dd}" + "/results/samples_trimmed/" + "{sample_id}" + "_R1_trimmed.fq.gz", dd=config["DATADIR"], sample_id=config["CONTROL"]),
#         ctr_rev = expand("{dd}" + "/results/samples_trimmed/" + "{sample_id}" + "_R2_trimmed.fq.gz", dd=config["DATADIR"], sample_id=config["CONTROL"]),
#         trt_fw = expand("{dd}" + "/results/samples_trimmed/" + "{sample_id}" + "_R1_trimmed.fq.gz", dd=config["DATADIR"], sample_id=config["TREATMENT"]),
#         trt_rev = expand("{dd}" + "/results/samples_trimmed/" + "{sample_id}" + "_R2_trimmed.fq.gz", dd=config["DATADIR"], sample_id=config["TREATMENT"]),
#     output:
#         outdir = expand("{dd}" + "/results/mapped_reads/human/",dd=config["DATADIR"]),
#         samout = expand("{dd}" + "/results/mapped_reads/human/" + "{sample_id}.Aligned.out.sam",dd=config["DATADIR"], sample_id=)
#     params:
#         gd = config["DATADIR"] + "/human_genome_indexed/",
#         prefix1 = 'control',
#         prefix2 = 'treatment',
#     conda:
#         "mapping.yml"
#     threads: config["THREADS"]    
#     shell:
#         """
#         STAR --runThreadN {threads} --genomeDir {params.gd} --readFilesIn {input.ctr_fw} {input.ctr_rev} --outFileNamePrefix {output.outdir}/{params.prefix1} --sjdbOverhang 100 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --quantMode GeneCounts &> log.out
#         STAR --runThreadN {threads} --genomeDir {params.gd} --readFilesIn {input.trt_fw} {input.trt_rev} --outFileNamePrefix {output.outdir}/{params.prefix2} --sjdbOverhang 100 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --quantMode GeneCounts
#         """
