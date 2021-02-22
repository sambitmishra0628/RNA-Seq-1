## Author: Sambit K. Mishra
## Created: 02-06-2021

configfile: "config.yml"
srcdir = "rules/"

include: srcdir + "fetch_samples.smk",
include: srcdir + "quality_control.smk",
include: srcdir + "fetch_references.smk",  
include: srcdir + "index_genomes.smk",
include: srcdir + "map_reads.smk",

rule all:
    input:
        # Fetch raw RNA-seq data for control and treatment
        expand("{dd}/samples/{control_id}_R{p}.fq.gz" , dd=config["DATADIR"], control_id=config["CONTROL"], p=[1,2]),
        expand("{dd}/samples/{treatment_id}_R{p}.fq.gz" , dd=config["DATADIR"], treatment_id=config["TREATMENT"], p=[1,2]),
                
        # Run FastQC on raw data
        directory(expand("{dd}/results/fastqc_raw/", dd=config["DATADIR"])),        

        # Run MultiQC on raw data
        directory(expand("{dd}/results/multiqc_raw/",dd=config["DATADIR"])),
        
        # Run cutadapt on the raw control and treatment data
        expand("{dd}" + "/results/samples_trimmed/" + "{control_id}_R" + "{p}_trimmed.fq.gz", dd=config['DATADIR'], control_id=config['CONTROL'], p=[1,2]),
        expand("{dd}" + "/results/samples_trimmed/" + "{treatment_id}_R" + "{p}_trimmed.fq.gz", dd=config['DATADIR'], treatment_id=config['TREATMENT'], p=[1,2]),
  
        # Run FastQC on trimmed data
        directory(expand("{dd}/results/fastqc_trimmed/", dd=config["DATADIR"])),

        # Run MultiQC on trimmed data
        directory(expand("{dd}/results/multiqc_trimmed/",dd=config["DATADIR"])),
        
        # Download the reference genomes and gffs
        directory(expand("{dd}" + "/reference_genomes/", dd=config["DATADIR"])),
        expand("{dd}" + "/reference_genomes/" + "{hg}", dd=config["DATADIR"], hg=config["human_genome"]),
        expand("{dd}" + "/reference_genomes/" + "{hg_gff}", dd=config["DATADIR"], hg_gff=config["human_gff"]),
        expand("{dd}" + "/reference_genomes/" + "{sg}", dd=config["DATADIR"], sg=config["sars2_genome"]),
        expand("{dd}" + "/reference_genomes/" + "{sg_gff}", dd=config["DATADIR"], sg_gff=config["sars2_gff"]),

        # # Index the reference genomes with STAR
        directory(expand("{dd}" + "/human_genome_indexed/", dd=config["DATADIR"])),
        directory(expand("{dd}" + "/sars2_genome_indexed/", dd=config["DATADIR"])),

        # Map reads to human genome using STAR
        #directory(expand("{dd}" + "/results/mapped_reads/human/", dd=config["DATADIR"]))




        