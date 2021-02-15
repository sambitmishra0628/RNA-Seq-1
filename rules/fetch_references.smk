## Author: Sambit K. Mishra
## Created: 02-14-2021
## Download the reference genomes and annotations for human and sars-2 corona virus
rule fetch_references:
    input:
        directory(config["DATADIR"])
    output:
        config["DATADIR"] +  "/reference_genomes/" + config["human_genome"],
        sars2_genome = config["DATADIR"] +  "/reference_genomes/" + config["sars2_genome"],
        human_gff = config["DATADIR"] +  "/reference_genomes/" + config["human_gff"],
        sars2_gff = config["DATADIR"] +  "/reference_genomes/" + config["sars2_gff"],
    params:
        hg_url = config["human_genome_url"],
        h_gff_url = config["human_gff_url"],
        cv_url = config["sars2_genome_url"],
        cv_gff_url = config["sars2_gff_url"],  
        hg_file = config["human_genome"],
        h_gff_file = config["human_gff"],
        cv_g_file = config["sars2_genome"],
        cv_gff_file = config["sars2_gff"], 
    shell:
        """
        mkdir {output}
        cd {output}
        wget {params.hg_url}
        wget {params.h_gff_url}
        wget {params.cv_url}
        wget {params.cv_gff_url}
        
        gunzip {params.hg_file}.gz
        gunzip {params.h_gff_file}.gz
        gunzip {params.cv_g_file}.gz
        gunzip {params.cv_gff_file}.gz
        
        """
