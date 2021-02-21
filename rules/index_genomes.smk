## Author: Sambit K. Mishra
## Created: 02-15-2021
## Index the genomes of human and sars-2 corona virus
rule index_genomes:
    input:
        human_genome = config["DATADIR"] +  "/reference_genomes/" + config["human_genome"],
        sars2_genome = config["DATADIR"] +  "/reference_genomes/" + config["sars2_genome"],
        human_gff = config["DATADIR"] +  "/reference_genomes/" + config["human_gff"],
        sars2_gff = config["DATADIR"] +  "/reference_genomes/" + config["sars2_gff"],
    output:
        out1 = directory(config["DATADIR"] + "/human_genome_indexed/"),
        out2 = directory(config["DATADIR"] + "/sars2_genome_indexed/"),
    conda:
        "mapping.yml"
    log:
        log1 = "human.log",
        log2 = "sars.log",    
    threads: config["THREADS"]
    params:
        human_gtf = config["DATADIR"] +  "/reference_genomes/" + config["human_gtf"],
        sars2_gtf = config["DATADIR"] +  "/reference_genomes/" + config["sars2_gtf"],
    shell:
        """
            mkdir -p {output.out1}
            mkdir -p {output.out2}
            gffread -T {input.human_gff} -o {params.human_gtf}
            gffread -T {input.sars2_gff} -o {params.sars2_gtf}
            STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.out1} --genomeFastaFiles {input.human_genome} --sjdbGTFfile {params.human_gtf} --sjdbOverhang 100 2> {output.out1}/{log.log1}
            STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.out2} --genomeFastaFiles {input.sars2_genome} --sjdbGTFfile {params.sars2_gtf} --sjdbOverhang 100 2> {output.out2}/{log.log2}
        """    
