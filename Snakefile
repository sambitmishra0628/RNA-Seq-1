## Author: Sambit K. Mishra
## Created: 02-06-2021

configfile: "config.yml"
srcdir = "rules/"

include: srcdir + "fetch_samples.smk",
include: srcdir + "quality_control.smk",    

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
        expand("{dd}" + "/results/samples_trimmed/" + "{control_id}_R1_trimmed.fq.gz", dd=config['DATADIR'], control_id=config['CONTROL']),
        expand("{dd}" + "/results/samples_trimmed/" + "{control_id}_R2_trimmed.fq.gz", dd=config['DATADIR'], control_id=config['CONTROL']),
        expand("{dd}" + "/results/samples_trimmed/" + "{treatment_id}_R1_trimmed.fq.gz", dd=config['DATADIR'], treatment_id=config['TREATMENT']),
        expand("{dd}" + "/results/samples_trimmed/" + "{treatment_id}_R2_trimmed.fq.gz", dd=config['DATADIR'], treatment_id=config['TREATMENT']),                
               
        # Run FastQC on trimmed data
        directory(expand("{dd}/results/fastqc_trimmed/", dd=config["DATADIR"])),

        # Run MultiQC on trimmed data
        directory(expand("{dd}/results/multiqc_trimmed/",dd=config["DATADIR"])),
        #directory(expand("{dd}" + "{fq}", dd=config["DATADIR"], fq=config["FQ_TRIMMED_DIR"])),
        #directory(expand("{dd}" + "{rr}/{fq}", dd=config["DATADIR"], rr=config["RESULTDIR_NAME"], fq=config["FQ_RAW"])),
        #expand("{sample_dir}" + "{control_id}" + "_R" + "{p}" + ".fq.gz.1" , sample_dir=config["SAMPLEDIR"], control_id=config["CONTROL"], p=[1,2]),
        #expand("{sample_dir}" + "{treatment_id}" + "_R" + "{p}" + ".fq.gz.1", sample_dir=config["SAMPLEDIR"], treatment_id=config["TREATMENT"], p=[1,2]),        
        #directory(expand("{sd}", sd=config["SAMPLEDIR"])),
        #expand("{sample_dir}" + "{control_id}",sample_dir=config["SAMPLEDIR"], control_id=config["CONTROL"]),
        #directory(expand("{sd}", sd=config["SAMPLEDIR"])), 
