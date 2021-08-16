rule index_cram:
  """index the subset CRAM file"""
    input:
        "../input/{pref}/cram/{sample}_mhc.cram"
    output:
        "../input/{pref}/cram/{sample}_mhc.cram.crai"
    conda:
        "../../envs/c4investigator.yaml"
    resources:
        mem_mb=2000,
        time=str(timedelta(minutes=30))
    threads:
        8
    shell:
        """samtools index -@ {threads} {input}"""


rule index_bam:
  """index the subset BAM file"""
    input:
        "../input/{pref}/bam/{sample}_mhc.bam"
    output:
        "../input/{pref}/bam/{sample}_mhc.bam.bai"
    conda:
        "../../envs/c4investigator.yaml"
    resources:
        mem_mb=2000,
        time=str(timedelta(minutes=30))
    threads:
        8
    shell:
        """samtools index -@ {threads} {input}"""



rule cram2fastq:
    """convert a cram to fastq suitable for C4Investigator"""
    input:
        ref_fq="../data/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        cramf="../input/{pref}/cram/{sample}_mhc.cram",
        crami="../input/{pref}/cram/{sample}_mhc.cram.crai",
    output:
        fq1="../input/{pref}/fastq/{sample}_MHC_1.fastq.gz",
        fq2="../input/{pref}/fastq/{sample}_MHC_2.fastq.gz"
    threads: 8
    conda:
        "../../envs/c4investigator.yaml"      
    shell:
        """bazam -Dsamjdk.reference_fasta={input.ref_fq} -Xmx4g  -bam {input.cramf} -r1 {output.fq1} -r2 {output.fq2}"""


rule bam2fastq:
    """convert a bam to fastq suitable for C4Investigator"""
    input:
        ref_fq="../data/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        cramf="../input/{pref}/bam/{sample}_mhc.bam",
        crami="../input/{pref}/bam/{sample}_mhc.bam.bai",
    output:
        fq1="../input/{pref}/fastq/{sample}_MHC_1.fastq.gz",
        fq2="../input/{pref}/fastq/{sample}_MHC_2.fastq.gz"
    threads: 8
    conda:
        "../../envs/c4investigator.yaml"      
    shell:
        """bazam -Dsamjdk.reference_fasta={input.ref_fq} -Xmx4g  -bam {input.cramf} -r1 {output.fq1} -r2 {output.fq2}"""
    
        
rule run_C4Investigator:
  """run C4Investigator"""
    input:
        fq1="../input/{pref}/fastq/{sample}_MHC_1.fastq.gz",
        fq2="../input/{pref}/fastq/{sample}_MHC_2.fastq.gz"
    output:
        "../output/{pref}/{sample}.tar.gz",
    resources:
        mem_mb=33000,
        time=str(timedelta(hours=3))
    params:
        workingDirectory="..",
        resultsDirectory="output/{pref}/{sample}",
        samplename="{sample}",
        projectName="C4Investigator",
        minDP="6",
        maxReadThreshold="50000"
    conda:
        "../../envs/c4investigator.yaml"
    log:
        "../output/{pref}/{sample}_log.txt"        
    threads:
        16
    script:
        """../scripts/C4Investigator_run.R"""
