rule star_index:
    input:
        bam = "../results/star/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "../results/star/{sample}/Aligned.sortedByCoord.out.bam.bai"
    conda:
        "../envs/invert.yaml"
    shell:
        '''
        samtools index {input} {output}
        '''

rule htseq:
    input:
        bam = "../results/star/{sample}/Aligned.sortedByCoord.out.bam",
        index = "../results/star/{sample}/Aligned.sortedByCoord.out.bam.bai"
    output:
        "../results/htseq/{sample}_counts.tsv"
    params:
        gtf=config["annotations"]["iav_only"],
    log:
        "logs/htseq/{sample}.log"
    conda:
        "../envs/htseq.yaml"
    shell:
        '''
        htseq-count {input.bam} {params.gtf} -t CDS -f bam -r pos -c {output} &> {log}
        '''