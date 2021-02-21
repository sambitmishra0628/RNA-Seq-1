## Author: Sambit K. Mishra
## Created: 02-06-2021

## Run quality control on raw RNA-seq data
rule fastqc_raw:
    input:
        ctr = expand("{dd}/samples/{control_id}_R{p}.fq.gz" , dd=config["DATADIR"], control_id=config["CONTROL"], p=[1,2]),
        trt = expand("{dd}/samples/{treatment_id}_R{p}.fq.gz" , dd=config["DATADIR"], treatment_id=config["TREATMENT"], p=[1,2]),
    output:
        out = directory(expand("{dd}/results/fastqc_raw/", dd=config["DATADIR"])),
    conda:
        "qc.yml"
    threads: config["THREADS"]
    shell:
        """
            mkdir -p {output.out}
            fastqc -o {output.out} -t {threads} {input.ctr} &> {output.out}/fastqc_raw_ctr.log
            fastqc -o {output.out} -t {threads} {input.trt} &> {output.out}/fastqc_raw_trt.log   
        """    
## Run MultiQC on raw RNA-seq data
rule multiqc_raw:
    input:
        inputdir = directory(expand("{dd}/results/fastqc_raw/", dd=config["DATADIR"])),
    output:
        outdir = directory(expand("{dd}/results/multiqc_raw/",dd=config["DATADIR"])),
    conda:
        "qc.yml"
    shell:
        "multiqc {input.inputdir} -o {output.outdir}"

# Remove adapters from raw control and treatment samples
rule trim_adapters:
    input:
        fw = "{dd}" + "/samples/" + "{sample_id}" + "_R1.fq.gz",
        rev = "{dd}" + "/samples/" + "{sample_id}" + "_R2.fq.gz",
    output:
        fw = "{dd}" + "/results/samples_trimmed/" + "{sample_id}" + "_R1_trimmed.fq.gz",
        rev = "{dd}" + "/results/samples_trimmed/" + "{sample_id}" + "_R2_trimmed.fq.gz",
    conda:
        "qc.yml"
    log: "{dd}" + "/results/samples_trimmed/" + "{sample_id}" + ".log"    
    threads: config["THREADS"]    
    shell:
        """
        cutadapt -q 20  -m 100 -j {threads} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o {output.fw} -p {output.rev} {input.fw} {input.rev} &> {log}
        """

# Perform fastqc quality check on trimmed data
rule fastqc_trimmed:
    input:
        ctr = expand("{dd}/results/samples_trimmed/{control_id}_R{p}_trimmed.fq.gz", dd=config["DATADIR"], control_id=config["CONTROL"], p=[1,2]),
        trt = expand("{dd}/results/samples_trimmed/{treatment_id}_R{p}_trimmed.fq.gz", dd=config["DATADIR"], treatment_id=config["TREATMENT"], p=[1,2]),
        #samples = "{dd}" + "/results/samples_trimmed/" + "{sample_id}" + "_R{p}_trimmed.fq.gz",
        #trt = expand("{dd}" + "{sd}" + "{treatment_id}" + "_R" + "{p}_trimmed" + ".fq.gz", dd=config["DATADIR"], sd=config["TRIMMED_SAMPLE_DIR"], treatment_id=config["TREATMENT"], p=[1,2]),
    output:
        out = directory(expand("{dd}/results/fastqc_trimmed/", dd=config["DATADIR"]))
    conda:
        "qc.yml"
    threads: config["THREADS"]    
    shell:
        """
        mkdir -p {output.out}
        fastqc -o {output.out} -t {threads} {input.ctr} &> {output.out}/fastqc_trimmed_ctr.log
        fastqc -o {output.out} -t {threads} {input.trt} &> {output.out}/fastqc_trimmed_trt.log
        """ 
 
# Perform multiqc quality check on trimmed reads
rule multiqc_trimmed:
    input:
        inputdir = directory(expand("{dd}/results/fastqc_trimmed/", dd=config["DATADIR"])),
    output:
        outdir = directory(expand("{dd}/results/multiqc_trimmed/",dd=config["DATADIR"])),
    conda:
        "qc.yml"
    shell:
        "multiqc {input.inputdir} -o {output.outdir}"



