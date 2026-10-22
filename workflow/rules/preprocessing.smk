###############################################################
# Preprocessing
###############################################################

rule fastp_PARIS_trim_adapter:
    input:
        r1 = lambda wc: samples.loc[wc.sample, "r1"]
    output:
        r1 = "results/fastp/{sample}_trim_adapter.fastq.gz",
        json = "results/fastp/{sample}_trim_adapter.json",
        html = "results/fastp/{sample}_trim_adapter.html"
    threads:
        8
    conda:
        "../envs/fastp.yaml"
    log:
        "results/logs/fastp/{sample}_trim_adapter.log"
    message:
        "Running fastp QC on PARIS sample {wildcards.sample}: Trim 3' adapters (Step 1/3)."
    wildcard_constraints:
        sample="|".join(paris_samples) if len(paris_samples) > 0 else "NOMATCH"
    shell:
        """
        fastp \
            -i {input.r1} \
            -o {output.r1} \
            --disable_quality_filtering \
            --disable_length_filtering \
            -j {output.json} -h {output.html} \
            --thread {threads} \
            2> {log}
        """

rule fastp_PARIS_dedup:
    input:
        r1 = rules.fastp_PARIS_trim_adapter.output.r1
    output:
        r1 = "results/fastp/{sample}_dedup.fastq.gz",
        json = "results/fastp/{sample}_dedup.json",
        html = "results/fastp/{sample}_dedup.html"
    threads:
        8
    conda:
        "../envs/fastp.yaml"
    log:
        "results/logs/fastp/{sample}_dedup.log"
    message:
        "Running fastp QC on PARIS sample {wildcards.sample}: Deduplicating (Step 2/3)."
    wildcard_constraints:
        sample="|".join(paris_samples) if len(paris_samples) > 0 else "NOMATCH"
    shell:
        """
        fastp \
            -i {input.r1} \
            -o {output.r1} \
            --dedup \
            --disable_adapter_trimming \
            --disable_quality_filtering \
            --disable_length_filtering \
            -j {output.json} -h {output.html} \
            --thread {threads} \
            2> {log}
        """

rule fastp_PARIS_trim_umi:
    input:
        r1 = rules.fastp_PARIS_dedup.output.r1
    output:
        r1 = "results/fastp/{sample}.fastq.gz",
        json = "results/fastp/{sample}_fastp.json",
        html = "results/fastp/{sample}_fastp.html"
    params:
        min_len = config.get("min_length", 20),
        min_qual = config.get("min_quality", 20),
        umi_trim_len = config.get("paris_umi_trim_len", 11),
    threads:
        8
    conda:
        "../envs/fastp.yaml"
    log:
        "results/logs/fastp/{sample}_trim_umi.log"
    message:
        "Running fastp QC on PARIS sample {wildcards.sample}: Trim 5' UMI/barcode (Step 3/3)."
    wildcard_constraints:
        sample="|".join(paris_samples) if len(paris_samples) > 0 else "NOMATCH"
    shell:
        """
        fastp \
            -i {input.r1} \
            -o {output.r1} \
            --trim_front1 {params.umi_trim_len} \
            -q {params.min_qual} \
            --length_required {params.min_len} \
            -j {output.json} -h {output.html} \
            --thread {threads} \
            2> {log}
        """

rule fastp_PARIS2_SPLASH_LIGRseq:
    input:
        r1 = lambda wc: samples.loc[wc.sample, "r1"]
    output:
        r1 = "results/fastp/{sample}.fastq.gz",
        json = "results/fastp/{sample}_fastp.json",
        html = "results/fastp/{sample}_fastp.html"
    params:
        min_len = config.get("min_length", 20),
        min_qual = config.get("min_quality", 20)
    threads:
        8
    conda:
        "../envs/fastp.yaml"
    log:
        "results/logs/fastp/{sample}.log"
    message:
        "Running fastp QC on PARIS2 sample {wildcards.sample}."
    wildcard_constraints:
        sample="|".join(paris2_samples.tolist() + splash_samples.tolist() + ligrseq_samples.tolist()) if len(paris2_samples) + len(splash_samples) + len(ligrseq_samples) > 0 else "NOMATCH"
    shell:
        """
        fastp \
            -i {input.r1} \
            -o {output.r1} \
            -q {params.min_qual} \
            --length_required {params.min_len} \
            -j {output.json} -h {output.html} \
            --thread {threads} \
            2> {log}
        """

rule fastp_KARRseq:
    input:
        r1 = lambda wc: samples.loc[wc.sample, "r1"],
        r2 = lambda wc: samples.loc[wc.sample, "r2"]
    output:
        merged = "results/fastp/{sample}_merged.fastq.gz",
        json = "results/fastp/{sample}_fastp.json",
        html = "results/fastp/{sample}_fastp.html"
    params:
        min_len = config.get("min_length", 20),
        min_qual = config.get("min_quality", 20)
    threads:
        8
    conda:
        "../envs/fastp.yaml"
    log:
        "results/logs/fastp/{sample}.log"
    message:
        "Running fastp QC on KARR-seq sample {wildcards.sample}: merging PE reads and deduplicating."
    wildcard_constraints:
        sample="|".join(karr_samples.tolist() + snokarr_samples.tolist()) if len(karr_samples) + len(snokarr_samples) > 0 else "NOMATCH"
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            --merge \
            --merged_out {output.merged} \
            --dedup \
            -q {params.min_qual} \
            --length_required {params.min_len} \
            -j {output.json} -h {output.html} \
            --thread {threads} \
            2> {log}
        """

rule fastp_RICseq:
    input:
        r1 = lambda wc: samples.loc[wc.sample, "r1"],
        r2 = lambda wc: samples.loc[wc.sample, "r2"]
    output:
        merged = "results/fastp/{sample}_merged.fastq.gz",
        json = "results/fastp/{sample}_fastp.json",
        html = "results/fastp/{sample}_fastp.html"
    params:
        min_len = config.get("min_length", 20),
        min_qual = config.get("min_quality", 20)
    threads:
        8
    conda:
        "../envs/fastp.yaml"
    log:
        "results/logs/fastp/{sample}.log"
    message:
        "Running fastp QC on RIC-seq sample {wildcards.sample}: merging PE reads, deduplicating, clipping poly(N) tails."
    wildcard_constraints:
        sample = "|".join(ricseq_samples.tolist()) if len(ricseq_samples) > 0 else "NOMATCH"
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            --merge --merged_out {output.merged} \
            --dedup \
            --trim_poly_x \
            -q {params.min_qual} \
            --length_required {params.min_len} \
            -j {output.json} -h {output.html} \
            --thread {threads} \
            2> {log}
        """