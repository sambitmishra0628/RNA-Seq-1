## Author: Sambit K. Mishra
## Created: 02-06-2021

## Download the samples for analysis
import os
import glob

rule fetch_samples:
input:
    directory(config["DATADIR"]),
output:
    expand("{dd}/samples/{control_id}_R{p}.fq.gz" , dd=config["DATADIR"], control_id=config["CONTROL"], p=[1,2]),
    expand("{dd}/samples/{treatment_id}_R{p}.fq.gz" , dd=config["DATADIR"], treatment_id=config["TREATMENT"], p=[1,2]),
run:
    # Create sample directory
    sample_dir = config["DATADIR"] + "/samples/"
    if not os.path.isdir(sample_dir):
        os.mkdir(sample_dir)        
    # Download control files
    print ("***Downloading control files ... ***")
    for srr,file_id in zip (config["SRR_CONTROL"],config["CONTROL"]):
        fw_file = config["src_url"] + srr + "/" + file_id + "_R1.fq.gz.1"
        rev_file = config["src_url"] + srr + "/" + file_id + "_R2.fq.gz.1"
        cmd1 = "wget " + fw_file + " -P " + sample_dir
        cmd2 = "wget " + rev_file + " -P " + sample_dir
        print (cmd1)
        print (cmd2)
        os.system(cmd1)
        os.system(cmd2)

    # Download treatment files
    print ("***Downloading treatment files ... ***")    
    for srr,file_id in zip (config["SRR_TREATMENT"],config["TREATMENT"]):
        fw_file = config["src_url"] + srr + "/" + file_id + "_R1.fq.gz.1"
        rev_file = config["src_url"] + srr + "/" + file_id + "_R2.fq.gz.1"
        cmd1 = "wget " + fw_file + " -P " + sample_dir
        cmd2 = "wget " + rev_file + " -P " + sample_dir
        print (cmd1)
        print (cmd2)
        os.system(cmd1)
        os.system(cmd2)            
    
    # Rename the files from *.fq.gz.1 to *.fq.gz
    curr_dir = os.getcwd()
    os.chdir(sample_dir)
    all_files = glob.glob("*.gz.1")
    for file_i in all_files:
        new_file = file_i.replace('.gz.1', '.gz')
        cmd = "mv " + file_i + " " + new_file
        os.system(cmd)
