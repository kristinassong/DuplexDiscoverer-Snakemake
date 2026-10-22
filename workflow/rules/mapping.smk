###############################################################
# Map reads to the genome
###############################################################

def get_star_input(wildcards):
    """Route correct fastp output to STAR — all samples enter as SE."""
    if wildcards.sample in pe_samples:
        return "results/fastp/{sample}_merged.fastq.gz".format(sample=wildcards.sample)
    else:
        return "results/fastp/{sample}.fastq.gz".format(sample=wildcards.sample)

rule STAR_index:
    input:
        fasta = config["genome_fasta"],
        gtf = config["genome_gtf"]
    output:
        genome = "results/star/genome/Genome"
    params:
        outdir = "results/star/genome",
        overhang = config.get("sjdb_overhang", 149)
    conda:
        "../envs/star.yaml"
    threads:
        8
    log:
        "results/logs/star/index.log"
    message:
        "Build STAR genome index."
    shell:
        """
        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {params.outdir} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.overhang} \
            2> {log}
        """

rule STAR_align_PARIS2_snoKARRseq_SPLASH_LIGRseq:
    input:
        fq = get_star_input,
        idx = rules.STAR_index.output.genome
    output:
        bam = "results/star/{sample}/{sample}_Aligned.out.bam",
        chim_junc = "results/star/{sample}/{sample}_Chimeric.out.junction",
        log_final = "results/star/{sample}/{sample}_Log.final.out",
        gene_counts = "results/star/{sample}/{sample}_ReadsPerGene.out.tab",
        sj = "results/star/{sample}/{sample}_SJ.out.tab"
    params:
        output_prefix = "results/star/{sample}/{sample}_",
        genomedir = "results/star/genome"
    conda:
        "../envs/star.yaml"
    threads:
        8
    log:
        "results/logs/star/{sample}.log"
    wildcard_constraints:
        sample="|".join(paris2_samples.tolist() + snokarr_samples.tolist()
                    + splash_samples.tolist() + ligrseq_samples.tolist())
            if len(paris2_samples) + len(snokarr_samples) + len(splash_samples) + len(ligrseq_samples) > 0 else "NOMATCH"
    message:
        "Align {wildcards.sample} reads to the reference genome using STAR."
    shell:
        """
        STAR \
            --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {params.genomedir} \
            --readFilesIn {input.fq} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.output_prefix} \
            --outSAMtype BAM Unsorted \
            --outSAMattributes All \
            --outReadsUnmapped Fastx \
            --genomeLoad NoSharedMemory \
            --outFilterMultimapNmax 100 \
            --outFilterScoreMinOverLread 0 \
            --outFilterMatchNminOverLread 0 \
            --outSJfilterReads All \
            --alignIntronMin 1 \
            --alignIntronMax 10 \
            --scoreGapNoncan 0 \
            --scoreGapGCAG 0 \
            --scoreGapATAC 0 \
            --scoreGenomicLengthLog2scale -1 \
            --chimOutType Junctions \
            --chimSegmentMin 5 \
            --chimJunctionOverhangMin 5 \
            --chimOutJunctionFormat 1 \
            --chimScoreJunctionNonGTAG 0 \
            --chimScoreDropMax 80 \
            --chimMultimapNmax 200 \
            --chimMainSegmentMultNmax 200 \
            --chimFilter banGenomicN \
            --limitOutSJcollapsed 10000000 \
            --limitIObufferSize 1500000000 1500000000 \
            --quantMode GeneCounts \
            2> {log}
        """

rule STAR_align_PARIS_KARRseq:
    input:
        fq = get_star_input,
        idx = rules.STAR_index.output.genome
    output:
        bam = "results/star/{sample}/{sample}_Aligned.out.bam",
        chim_junc = "results/star/{sample}/{sample}_Chimeric.out.junction",
        log_final = "results/star/{sample}/{sample}_Log.final.out",
        gene_counts = "results/star/{sample}/{sample}_ReadsPerGene.out.tab",
        sj = "results/star/{sample}/{sample}_SJ.out.tab"
    params:
        output_prefix = "results/star/{sample}/{sample}_",
        genomedir = "results/star/genome"
    conda:
        "../envs/star.yaml"
    threads:
        8
    log:
        "results/logs/star/{sample}.log"
    wildcard_constraints:
        sample="|".join(paris_samples.tolist() + karr_samples.tolist()) \
            if len(paris_samples) + len(karr_samples) > 0 else "NOMATCH"
    message:
        "Align {wildcards.sample} reads to the reference genome using STAR."
    shell:
        """
        STAR \
            --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {params.genomedir} \
            --readFilesIn {input.fq} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.output_prefix} \
            --outSAMtype BAM Unsorted \
            --outSAMattributes All \
            --outReadsUnmapped Fastx \
            --genomeLoad NoSharedMemory \
            --outFilterMultimapNmax 100 \
            --outFilterScoreMinOverLread 0 \
            --outFilterMatchNminOverLread 0 \
            --outSJfilterReads All \
            --alignIntronMin 1 \
            --alignIntronMax 10 \
            --scoreGapNoncan -4 \
            --scoreGapGCAG -4 \
            --scoreGapATAC -4 \
            --scoreGenomicLengthLog2scale -0.25 \
            --chimOutType Junctions \
            --chimSegmentMin 15 \
            --chimJunctionOverhangMin 15 \
            --chimOutJunctionFormat 1 \
            --chimScoreJunctionNonGTAG -1 \
            --chimScoreDropMax 20 \
            --chimMultimapNmax 200 \
            --chimMainSegmentMultNmax 200 \
            --chimFilter banGenomicN \
            --limitOutSJcollapsed 10000000 \
            --limitIObufferSize 1500000000 1500000000 \
            --quantMode GeneCounts \
            2> {log}
        """

rule STAR_align_RICseq:
    input:
        fq = get_star_input,
        idx = rules.STAR_index.output.genome
    output:
        bam = "results/star/{sample}/{sample}_Aligned.out.bam",
        chim_junc = "results/star/{sample}/{sample}_Chimeric.out.junction",
        log_final = "results/star/{sample}/{sample}_Log.final.out",
        gene_counts = "results/star/{sample}/{sample}_ReadsPerGene.out.tab",
        sj = "results/star/{sample}/{sample}_SJ.out.tab"
    params:
        output_prefix = "results/star/{sample}/{sample}_",
        genomedir = "results/star/genome"
    conda:
        "../envs/star.yaml"
    threads:
        8
    log:
        "results/logs/star/{sample}.log"
    wildcard_constraints:
        sample="|".join(ricseq_samples.tolist()) if len(ricseq_samples) > 0 else "NOMATCH"
    message:
        "Align {wildcards.sample} reads to the reference genome using STAR."
    shell:
        """
        STAR \
            --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {params.genomedir} \
            --readFilesIn {input.fq} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.output_prefix} \
            --outSAMtype BAM Unsorted \
            --outSAMattributes All \
            --outReadsUnmapped Fastx \
            --genomeLoad NoSharedMemory \
            --outFilterMultimapNmax 100 \
            --outFilterScoreMinOverLread 0 \
            --outFilterMatchNminOverLread 0 \
            --outSJfilterReads All \
            --alignIntronMin 1 \
            --alignIntronMax 10 \
            --scoreGapNoncan -4 \
            --scoreGapGCAG -4 \
            --scoreGapATAC -4 \
            --scoreGenomicLengthLog2scale -0.25 \
            --chimOutType Junctions \
            --chimSegmentMin 15 \
            --chimJunctionOverhangMin 15 \
            --alignSJoverhangMin 15 \
            --alignSJDBoverhangMin 10 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --chimOutJunctionFormat 1 \
            --chimScoreJunctionNonGTAG -1 \
            --chimScoreDropMax 20 \
            --chimMultimapNmax 200 \
            --chimMainSegmentMultNmax 200 \
            --chimFilter banGenomicN \
            --limitOutSJcollapsed 10000000 \
            --limitIObufferSize 1500000000 1500000000 \
            --quantMode GeneCounts \
            2> {log}
        """

rule flip_strands:
    input:
        chim_junc = "results/star/{sample}/{sample}_Chimeric.out.junction"
    output:
        flipped = temp("results/star/{sample}/{sample}_Chimeric.flipped.junction")
    wildcard_constraints:
        sample = "|".join(sorted(REVERSE_STRANDED)) if len(REVERSE_STRANDED) > 0 else "NOMATCH"
    log:
        "results/logs/star/{sample}_flip_strands.log"
    message:
        "Flip donor and acceptor strands for {wildcards.sample} (reverse-stranded library)."
    shell:
        """
        awk 'BEGIN{{ FS = OFS = "\\t" }}
             $2 !~ /^[0-9]+$/ {{ print; next }}
             {{ $3 = ($3 == "+" ? "-" : "+"); $6 = ($6 == "+" ? "-" : "+"); print }}
        ' {input.chim_junc} > {output.flipped} 2> {log}
        """

def get_chim_junc(wildcards):
    tmpl = ("results/star/{s}/{s}_Chimeric.flipped.junction"
            if wildcards.sample in REVERSE_STRANDED
            else "results/star/{s}/{s}_Chimeric.out.junction")
    return tmpl.format(s=wildcards.sample)

rule pad_snoRNA_bed:
    input:
        bed = config["snorna_bed"]
    output:
        bed = "resources/snoRNA_padded.bed"
    params:
        pad = config.get("snorna_pad", 25)
    log:
        "results/logs/star/pad_snoRNA_bed.log"
    message:
        "Pad the input snoRNA BED by {params.pad} nt on each side for junction prefiltering."
    shell:
        """
        {{
        awk -F'\\t' -v pad={params.pad} 'BEGIN{{ OFS = "\\t" }}
            /^(#|track|browser)/ {{ next }}
            NF >= 3 {{
                s = $2 - pad; if (s < 0) s = 0;
                $2 = s; $3 = $3 + pad;
                print
            }}' {input.bed} \
            | sort -k1,1 -k2,2n > {output.bed}
        }} 2> {log}
        """

rule filter_snoRNA_junctions:
    input:
        chim_junc = get_chim_junc,
        sno_bed   = rules.pad_snoRNA_bed.output.bed
    output:
        sno_junc = "results/star/{sample}/{sample}_Chimeric.sno.junction",
        read_placements  = "results/star/{sample}/{sample}_read_placements.tsv"
    params:
        bedpe    = "results/star/{sample}/{sample}.tmp.bedpe",
        keep     = "results/star/{sample}/{sample}.tmp.keep",
        filtered = "results/star/{sample}/{sample}.tmp.filtered"
    conda:
        "../envs/bedtools.yaml"
    log:
        "results/logs/star/{sample}_filter_snoRNA_junctions.log"
    message:
        "Keep chimeric alignments with an arm on a snoRNA or scaRNA, then retain the best-scoring placement per read for {wildcards.sample}."
    shell:
        """
        {{
        awk 'BEGIN{{ OFS = "\\t" }} $2 ~ /^[0-9]+$/ {{
            print $1, $2, $2 + 1, $4, $5, $5 + 1, NR
        }}' {input.chim_junc} > {params.bedpe}

        bedtools pairtobed -a {params.bedpe} -b {input.sno_bed} -type either \
            | cut -f7 | sort -un > {params.keep}

        awk -v kf={params.keep} '
            FILENAME == kf {{ k[$1] = 1; next }}
            /^#/ {{ next }}
            $2 !~ /^[0-9]+$/ {{ print; next }}
            (FNR in k)
        ' {params.keep} {input.chim_junc} > {params.filtered}

        awk '$2 ~ /^[0-9]+$/ {{ print $10 "\\t" $15 }}' {params.filtered} \
            | sort -u > {output.read_placements}

        awk 'BEGIN{{ FS = OFS = "\\t" }}
             $2 !~ /^[0-9]+$/ {{ print; next }}
             $1 !~ /^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$/ {{ next }}
             $4 !~ /^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$/ {{ next }}
             $18 != $19 {{ next }}
             {{ n[$10]++; $10 = $10 "_" n[$10]; $15 = 1; print }}
        ' {params.filtered} > {output.sno_junc}
        }} 2> {log}

        rm -f {params.bedpe} {params.keep} {params.filtered}
        """

rule multiqc:
    input:
        fastp = expand("results/fastp/{sample}_fastp.json", sample=samples.index),
        star = expand("results/star/{sample}/{sample}_Log.final.out", sample=samples.index)
    output:
        "results/multiqc/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    log:
        "results/logs/multiqc/multiqc.log"
    message:
        "Aggregate QC reports with MultiQC."
    shell:
        """
        multiqc results/fastp/ results/star/ \
            --outdir results/multiqc/ \
            --filename multiqc_report.html \
            2> {log}
        """
