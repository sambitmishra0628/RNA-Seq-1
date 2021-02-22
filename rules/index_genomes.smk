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
        out1 = directory(expand("{dd}" + "/human_genome_indexed/", dd=config["DATADIR"])),
        out2 = directory(expand("{dd}" + "/sars2_genome_indexed/", dd=config["DATADIR"])),
    conda:
        "mapping.yml"
    threads: config["THREADS"]
    params:
        human_gtf = config["DATADIR"] +  "/reference_genomes/" + config["human_gtf"],
        sars2_gtf = config["DATADIR"] +  "/reference_genomes/" + config["sars2_gtf"],
        prefix1 = 'human',
        prefix2 = 'sars',
    shell:
        """
        gffread -T {input.human_gff} -o {params.human_gtf}
        gffread -T {input.sars2_gff} -o {params.sars2_gtf}
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.out1} --genomeFastaFiles {input.human_genome} \
        --sjdbGTFfile {params.human_gtf} --outFileNamePrefix {output.out1}/{params.prefix1} --sjdbOverhang 100 \
        --limitGenomeGenerateRAM 25000000000 --genomeChrBinNbits=16 --genomeSAsparseD 3
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.out2} --genomeFastaFiles {input.sars2_genome} --sjdbGTFfile {params.sars2_gtf} --outFileNamePrefix {output.out2}/{params.prefix2} --sjdbOverhang 100
        """    
#             mkdir -p {output.out1}
#            mkdir -p {output.o