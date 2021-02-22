## Author: Sambit K. Mishra
## Created: 02-21-2021

## Map reads to human genome using STAR
rule map_reads:
    input:
        ctr_fw = expand("{dd}" + "/results/samples_trimmed/" + "{sample_id}" + "_R1_trimmed.fq.gz", dd=config["DATADIR"], sample_id=config["CONTROL"]),
        ctr_rev = expand("{dd}" + "/results/samples_trimmed/" + "{sample_id}" + "_R2_trimmed.fq.gz", dd=config["DATADIR"], sample_id=config["CONTROL"]),
        trt_fw = expand("{dd}" + "/results/samples_trimmed/" + "{sample_id}" + "_R2_trimmed.fq.gz", dd=config["DATADIR"], sample_id=config["TREATMENT"]),
        trt_rev = expand("{dd}" + "/results/samples_trimmed/" + "{sample_id}" + "_R2_trimmed.fq.gz", dd=config["DATADIR"], sample_id=config["TREATMENT"]),
    output:
        outdir = expand("{dd}" + "/results/mapped_reads/human/",dd=config["DATADIR"]),
    params:
        gd = config["DATADIR"] + "/human_genome_indexed/",
        prefix1 = 'control',
        prefix2 = 'treatment',
    conda:
        "mapping.yml"
    threads: config["THREADS"]    
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {params.gd} --readFilesIn {input.ctr_fw} {input.ctr_rev} --outFileNamePrefix {output.outdir}/{params.prefix1} --sjdbOverhang 100 --readFilesCommand gunzip -c
        STAR --runThreadN {threads} --genomeDir {params.gd} --readFilesIn {input.trt_fw} {input.trt_rev} --outFileNamePrefix {output.outdir}/{params.prefix2} --sjdbOverhang 100 --readFilesCommand gunzip -c
        """
