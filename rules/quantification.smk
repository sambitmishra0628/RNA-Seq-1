## Quantify mapped reads from STAR using RSEM
## Author: Sambit K. Mishra
## Created: 02-26-2021
rule prepare_reference_human:
    input:
        directory(config['DATADIR'] + "/reference_genomes"),
    output:
        outdir = directory(config['DATADIR'] + "/reference_genomes/rsem_ref/"),
    params:
        ref_dir = config['DATADIR'] + "/reference_genomes",
        gtf_file = config["human_gtf"],
        genome_file = config["human_genome"],
    conda:
        "quantification.yml"
    shell:
        """
        rsem-prepare-reference --gtf {ref_dir}/{gtf_file} -trusted-sources BestRefSeq,Curated\ Genomic {params.ref_dir}/{params.genome_file} {output.outdir}/human_refseq
        """

rule quantify_reads:
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