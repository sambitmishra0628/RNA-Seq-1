## Author: Sambit K. Mishra
## Created: 02-06-2021

configfile: "config.yml"
srcdir = "snakefiles/"

include:
    srcdir + "fetch_samples.smk"

rule all:
    input:
        expand("{sample_dir}" + "{control_id}" + "_R" + "{p}" + ".fq.gz.1" , sample_dir=config["SAMPLEDIR"], control_id=config["CONTROL"], p=[1,2]),
        expand("{sample_dir}" + "{treatment_id}" + "_R" + "{p}" + ".fq.gz.1", sample_dir=config["SAMPLEDIR"], treatment_id=config["TREATMENT"], p=[1,2]),
        #directory(expand("{sd}", sd=config["SAMPLEDIR"])),
        #expand("{sample_dir}" + "{control_id}",sample_dir=config["SAMPLEDIR"], control_id=config["CONTROL"]),
        #directory(expand("{sd}", sd=config["SAMPLEDIR"])), 
