rule strand_separation:
    message: "Extracting forward and reverse strand-specific reads from STAR mapping step."
    input:
        bam="../results/star/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        fwd1_bam = temp("../results/invert_by_group/{group}/{sample}_fwd1.bam"),
        fwd2_bam = temp("../results/invert_by_group/{group}/{sample}_fwd2.bam"),
        fwd_bam = "../results/invert_by_group/{group}/{sample}_fwd.bam",
        rev1_bam = temp("../results/invert_by_group/{group}/{sample}_rev1.bam"),
        rev2_bam = temp("../results/invert_by_group/{group}/{sample}_rev2.bam"),
        rev_bam = ".../results/invert_by_group/{group}/{sample}_rev.bam"
    conda:
        "../envs/invert.yaml"
    shell:
        """
        #Forward strand
        # 1. alignments of the second in pair if they map to the forward strand
        samtools view -b -f 128 -F 16 {input.bam} > {output.fwd1_bam}

        # 2. alignments of the first in pair if they map to the reverse  strand 
        samtools view -b -f 80 {input.bam} > {output.fwd2_bam}
        
        # Combine alignments that originate on the forward strand.
        samtools merge -f {output.fwd_bam} {output.fwd1_bam} {output.fwd2_bam}
        samtools index {output.fwd_bam}

        # Reverse strand
        # 1. alignments of the second in pair if they map to the reverse strand
        samtools view -b -f 144 {input.bam} > {output.rev1_bam}

        # 2. alignments of the first in pair if they map to the forward strand
        samtools view -b -f 64 -F 16 {input.bam} > {output.rev2_bam}

        # Combine alignments that originate on the reverse strand.
        samtools merge -f {output.rev_bam} {output.rev1_bam} {output.rev2_bam}
        samtools index {output.rev_bam}
        """

rule segment_specific_bam:
    message: "Extracting bam files of forward reads mapping to each segment of IAV."
    input:
        bam=rules.strand_separation.output.fwd_bam
    output:
        expand("../results/invert_by_group/{{group}}/{{sample}}_{segment}.bam", segment = ['PB1','PB2','PA','NP','HA','NA','NS','M'])
    params:
        out_file_path="../results/invert_by_group/{group}"
    conda:
        "../envs/invert.yaml"
    shell:
        '''
        scripts/segment_bams.sh {input.bam} {params.out_file_path} {wildcards.sample}
        '''

rule combine_bams_by_group:
    input:
        lambda wildcards: expand("../results/invert_by_group/{group}/{sample}_{{segment}}.bam", group = wildcards.group, sample=GROUPS[wildcards.group])
    output:
        "../results/invert_by_group/{group}/{segment}_all.bam"
    conda:
        "../envs/invert.yaml"
    shell:
        '''
        samtools merge {output} {input}
        '''

rule calculate_cmrna_ratios_by_group:
    input:
        expand("../results/invert_by_group/{{group}}/{segment}_all.bam", segment = ['PB1','PB2','PA','NP','HA','NA','NS','M'])
    output:
        "../results/invert_by_group/{group}/cmrna_ratios.txt"
    params:
        in_file_dir="../results/invert_by_group/{group}"
    conda:
        "../envs/invert.yaml"
    shell:
        '''
        scripts/cmrna_count_group.sh {params.in_file_dir} {output}
        '''

rule splice_ratios:
    input:
        bam=rules.strand_separation.output.fwd_bam
    output:
        splice_counts="../results/invert_by_group/{group}/{sample}_splicing_count.txt"
    conda:
        "../envs/invert.yaml"
    shell:
        """
        mkdir -p $(dirname {output.splice_counts}/)
        scripts/cmrna_splicing.sh {input.bam} {wildcards.sample} {output.splice_counts}
        """

rule kinetics_calculations:
    input:
        ratios=lambda wildcards: expand("../results/invert_by_group/{group}/cmrna_ratios.txt", group=wildcards.group),
        splice_counts= "../results/invert_by_group/{group}/{sample}_splicing_count.txt",
        expression_levels="../results/cufflinks/{sample}/genes.fpkm_tracking"
    output:
        "../results/invert_by_group/{group}/{sample}_final_results.csv"
    shell:
        """
        python scripts/vcmRNA_calculation.py {input.ratios} {input.splice_counts} {input.expression_levels} {output}
        """