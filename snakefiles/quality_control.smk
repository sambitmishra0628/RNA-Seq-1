## Author: Sambit K. Mishra
## Created: 02-06-2021

## Run quality control on raw RNA-seq data
rule fastqc_raw:
    input:
        ctr = expand("{dd}" + "{sd}" + "{control_id}" + "_R" + "{p}" + ".fq.gz" , dd=config["DATADIR"], sd=config["SAMPLEDIR"], control_id=config["CONTROL"], p=[1,2]),
        trt = expand("{dd}" + "{sd}" + "{treatment_id}" + "_R" + "{p}" + ".fq.gz", dd=config["DATADIR"], sd=config["SAMPLEDIR"], treatment_id=config["TREATMENT"], p=[1,2]),
    output:
        out = directory(expand("{dd}" + "{fq}", dd=config["DATADIR"], fq=config["FQ_RAW_DIR"]))
    conda:
        "qc.yml"
    threads: 5
    shell:
        """
            mkdir -p {output.out}
            fastqc -o {output.out} -t {threads} {input.ctr}
            fastqc -o {output.out} -t {threads} {input.trt}    
        """    
## Run MultiQC on raw RNA-seq data
rule multiqc_raw:
    input:
        inputdir = directory(expand("{dd}" + "{fq}", dd=config["DATADIR"], fq=config["FQ_RAW_DIR"]))
    output:
        outdir = directory(expand("{dd}" + "{fq}", dd=config["DATADIR"], fq=config["MQ_RAW_DIR"]))
    conda:
        "qc.yml"
    shell:
        "multiqc {input.inputdir} -o {output.outdir}"


rule trim_adapters:
    input:
        ctr_fw = expand("{dd}" + "{sd}" + "{control_id}" + "_R1"  + ".fq.gz" , dd=config["DATADIR"], sd=config["SAMPLEDIR"], control_id=config["CONTROL"]),
        ctr_rev = expand("{dd}" + "{sd}" + "{control_id}" + "_R2"  + ".fq.gz" , dd=config["DATADIR"], sd=config["SAMPLEDIR"], control_id=config["CONTROL"]),
        trt_fw = expand("{dd}" + "{sd}" + "{treatment_id}" + "_R1" + ".fq.gz", dd=config["DATADIR"], sd=config["SAMPLEDIR"], treatment_id=config["TREATMENT"]),
        trt_rev = expand("{dd}" + "{sd}" + "{treatment_id}" + "_R2" + ".fq.gz", dd=config["DATADIR"], sd=config["SAMPLEDIR"], treatment_id=config["TREATMENT"]),
    output:
        ctr_fw = expand("{dd}" + "{td}" + "{control_id}" + "_R1_trimmed"  + ".fq.gz" , dd=config["DATADIR"], td=config["TRIMMED_SAMPLE_DIR"], control_id=config["CONTROL"]),
        ctr_rev = expand("{dd}" + "{td}" + "{control_id}" + "_R2_trimmed"  + ".fq.gz" , dd=config["DATADIR"], td=config["TRIMMED_SAMPLE_DIR"], control_id=config["CONTROL"]),
        trt_fw = expand("{dd}" + "{td}" + "{treatment_id}" + "_R1_trimmed" + ".fq.gz", dd=config["DATADIR"], td=config["TRIMMED_SAMPLE_DIR"], treatment_id=config["TREATMENT"]),
        trt_rev = expand("{dd}" + "{td}" + "{treatment_id}" + "_R2_trimmed" + ".fq.gz", dd=config["DATADIR"], td=config["TRIMMED_SAMPLE_DIR"], treatment_id=config["TREATMENT"]),
    conda:
        "qc.yml"
    threads: 4    
    shell:
        """
        cutadapt -q 20 --cores {threads} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o {output.ctr_fw} -p {output.ctr_rev} {input.ctr_fw} {input.ctr_rev}
        cutadapt -q 20 --cores {threads} -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o {output.trt_fw} -p {output.trt_rev} {input.trt_fw} {input.trt_rev} 
        """

# rule fastqc_trimmed:
#     input:

#     output:

#     conda:

#     shell:



# rule multiqc_trimmed:
#     input:

#     output:

#     conda:

#     shell:


