rule star:
    input:
        fq1 = get_r1,
        fq2 = get_r2
    params:
        index = config["genome_index"],
        gtf = config["annotations"]["combined"],
        outdir = '../results/star/{sample}'
    output:
        bam = '../results/star/{sample}/{sample}_Aligned.sortedByCoord.out.bam'
    log:
        "logs/star/{sample}.log"
    conda:
        "../envs/star.yaml"
    threads:
        16
    shell:
        """
        STAR --runThreadN {threads} \
            --runMode alignReads \
            --genomeDir {params.index} \
            --readFilesIn {input.fq1} {input.fq2} \
            --outSAMtype BAM SortedByCoordinate \
            --sjdbOverhang 149 \
            --outFilterType BySJout \
            --outFilterMultimapNmax 10 \
            --alignSJoverhangMin 5 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.outdir}/{wildcards.sample}_ &> {log}
        
        rm -rf ../results/star/*/*_STARtmp
        """