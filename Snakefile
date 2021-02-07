## Author: Sambit K. Mishra
## Created: 02-06-2021

configfile: "config.yml"
srcdir = "snakefiles/"

include: srcdir + "fetch_samples.smk",
include: srcdir + "quality_control.smk",    

rule all:
    input:
        # Fetch raw RNA-seq data for control and treatment
        expand("{dd}" + "{sd}" + "{control_id}" + "_R" + "{p}" + ".fq.gz" , dd=config["DATADIR"], sd=config["SAMPLEDIR"], control_id=config["CONTROL"], p=[1,2]),
        expand("{dd}" + "{sd}" + "{treatment_id}" + "_R" + "{p}" + ".fq.gz", dd=config["DATADIR"], sd=config["SAMPLEDIR"], treatment_id=config["TREATMENT"], p=[1,2]),
        
        # Run FastQC on raw data
        directory(expand("{dd}" + "{fq}", dd=config["DATADIR"], fq=config["FQ_RAW_DIR"])),

        # Run MultiQC on raw data
        directory(expand("{dd}" + "{fq}", dd=config["DATADIR"], fq=config["MQ_RAW_DIR"]))
        
        
        #directory(expand("{dd}" + "{rr}/{fq}", dd=config["DATADIR"], rr=config["RESULTDIR_NAME"], fq=config["FQ_RAW"])),
        #expand("{sample_dir}" + "{control_id}" + "_R" + "{p}" + ".fq.gz.1" , sample_dir=config["SAMPLEDIR"], control_id=config["CONTROL"], p=[1,2]),
        #expand("{sample_dir}" + "{treatment_id}" + "_R" + "{p}" + ".fq.gz.1", sample_dir=config["SAMPLEDIR"], treatment_id=config["TREATMENT"], p=[1,2]),        
        #directory(expand("{sd}", sd=config["SAMPLEDIR"])),
        #expand("{sample_dir}" + "{control_id}",sample_dir=config["SAMPLEDIR"], control_id=config["CONTROL"]),
        #directory(expand("{sd}", sd=config["SAMPLEDIR"])), 
