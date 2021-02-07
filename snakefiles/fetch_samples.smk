## Author: Sambit K. Mishra
## Created: 02-06-2021

## Download the samples for analysis
import os
rule fetch_samples:
    input:
        directory(config["DATADIR"]),
    output:
        expand("{sample_dir}" + "{control_id}" + "_R" + "{p}" + ".fq.gz.1" , sample_dir=config["SAMPLEDIR"], control_id=config["CONTROL"], p=[1,2]),
        expand("{sample_dir}" + "{treatment_id}" + "_R" + "{p}" + ".fq.gz.1", sample_dir=config["SAMPLEDIR"], treatment_id=config["TREATMENT"], p=[1,2]),
    run:
        # Create sample directory
        #if not os.path.isdir(config["SAMPLEDIR"]):
        #    os.mkdir(config["SAMPLEDIR"])        
        # Download control files
        print ("***Downloading control files ... ***")
        for srr,file_id in zip (config["SRR_CONTROL"],config["CONTROL"]):
            fw_file = config["src_url"] + srr + "/" + file_id + "_R1.fq.gz.1"
            rev_file = config["src_url"] + srr + "/" + file_id + "_R2.fq.gz.1"
            cmd1 = "wget " + fw_file + " -P " + config["SAMPLEDIR"]
            cmd2 = "wget " + rev_file + " -P " + config["SAMPLEDIR"]
            print (cmd1)
            print (cmd2)
            os.system(cmd1)
            os.system(cmd2)
        # Download treatment files
        print ("***Downloading treatment files ... ***")    
        for srr,file_id in zip (config["SRR_TREATMENT"],config["TREATMENT"]):
            fw_file = config["src_url"] + srr + "/" + file_id + "_R1.fq.gz.1"
            rev_file = config["src_url"] + srr + "/" + file_id + "_R2.fq.gz.1"
            cmd1 = "wget " + fw_file + " -P " + config["SAMPLEDIR"]
            cmd2 = "wget " + rev_file + " -P " + config["SAMPLEDIR"]
            print (cmd1)
            print (cmd2)
            os.system(cmd1)
            os.system(cmd2)            

        #"""
        #wget {params.control_files[1]} -P {params.sample_dir}
        #"""
        #wget {params.src_url}{params.srr_control}/{params.control_id} -P {output.cf}
        #wget {params.src_url}{params.srr_control}/{params.control_id} -P {params.sample_dir}
        #wget {params.src_url}{params.srr_control}/{params.control_id}"_R2.fq.gz.1" -P {params.sample_dir}
        #expand("{sample_dir}" + "{control_id}" + "_R" + "{p}" + ".fq.gz.1" , sample_dir=config["SAMPLEDIR"], control_id=config["CONTROL"], p=[1,2]),
        #expand("{sample_dir}" + "{treatment_id}" + "_R" + "{p}" + ".fq.gz.1", sample_dir=config["SAMPLEDIR"], treatment_id=config["TREATMENT"], p=[1,2]),        