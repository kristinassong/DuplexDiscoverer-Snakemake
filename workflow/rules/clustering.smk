###############################################################
# DuplexDiscovereR
###############################################################

rule install_duplexdiscoverer:
    output:
        touch("resources/.duplexdiscoverer_installed"),
    conda:
        "../envs/r.yaml"
    log:
        "results/logs/install_duplexdiscoverer.log",
    message:
        "Installing DuplexDiscovereR from Bioconductor git."
    shell:
        """
        rm -rf resources/DuplexDiscovereR
        git clone https://git.bioconductor.org/packages/DuplexDiscovereR resources/DuplexDiscovereR 2> {log}
        sed -i 's/GenomicAlignments::cigarWidthAlongReferenceSpace/cigarillo::cigar_extent_along_ref/g' resources/DuplexDiscovereR/R/utils_preproc.R
        sed -i 's/cigarWidthAlongReferenceSpace(/cigarillo::cigar_extent_along_ref(/g' resources/DuplexDiscovereR/R/utils_preproc.R
        sed -i 's/@importFrom GenomicAlignments cigarWidthAlongReferenceSpace/@importFrom cigarillo cigar_extent_along_ref/' resources/DuplexDiscovereR/R/DuplexDiscovereR-package.R
        R CMD INSTALL resources/DuplexDiscovereR >> {log} 2>&1
        """

rule run_duplexdiscoverer:
    input:
        chim = rules.filter_snoRNA_junctions.output.sno_junc,
        sj = "results/star/{sample}/{sample}_SJ.out.tab",
        gene_counts = "results/star/{sample}/{sample}_ReadsPerGene.out.tab",
        gtf = config["genome_gtf"],
        installed = rules.install_duplexdiscoverer.output
    output:
        dg = "results/duplexdiscoverer/{sample}/{sample}_dgs.tsv",
        stats = "results/duplexdiscoverer/{sample}/{sample}_chimeric_stats.tsv",
        rds = "results/duplexdiscoverer/{sample}/{sample}_dgs.rds"
    log:
        "results/logs/duplexdiscoverer/{sample}.log"
    wildcard_constraints:
        sample = "|".join(sorted(NORMAL_MEM)) if NORMAL_MEM else "NOMATCH"
    params:
        lib = "SE",
        collapse_n_inter = config["collapse_n_inter"],
        gap_collapse_similar = config["gap_collapse_similar"],
        compute_p_values = config["compute_p_values"],
        min_arm_len = lambda wc: 15 if (wc.sample in paris_samples or wc.sample in karr_samples or wc.sample in ricseq_samples) else 5,
        trim_alignments = config["trim_alignments"],
        trim_length = config["trim_length"]
    threads:
        4
    conda:
        "../envs/r.yaml"
    message:
        "Run DuplexDiscovereR for {wildcards.sample}."
    script:
        "../scripts/duplexdiscoverer.R"

rule run_duplexdiscoverer_large:
    input:
        chim = rules.filter_snoRNA_junctions.output.sno_junc,
        sj = "results/star/{sample}/{sample}_SJ.out.tab",
        gene_counts = "results/star/{sample}/{sample}_ReadsPerGene.out.tab",
        gtf = config["genome_gtf"],
        installed = rules.install_duplexdiscoverer.output
    output:
        dg = "results/duplexdiscoverer/{sample}/{sample}_dgs.tsv",
        stats = "results/duplexdiscoverer/{sample}/{sample}_chimeric_stats.tsv",
        rds = "results/duplexdiscoverer/{sample}/{sample}_dgs.rds"
    log:
        "results/logs/duplexdiscoverer/{sample}.log"
    wildcard_constraints:
        sample = "|".join(sorted(HIGH_MEM)) if HIGH_MEM else "NOMATCH"
    params:
        lib = "SE",
        collapse_n_inter = config["collapse_n_inter"],
        gap_collapse_similar = config["gap_collapse_similar"],
        compute_p_values = config["compute_p_values"],
        min_arm_len = lambda wc: 15 if (wc.sample in paris_samples or wc.sample in karr_samples or wc.sample in ricseq_samples) else 5,
        trim_alignments = config["trim_alignments"],
        trim_length = config["trim_length"]
    threads:
        4
    conda:
        "../envs/r.yaml"
    message:
        "Run DuplexDiscovereR for {wildcards.sample}."
    script:
        "../scripts/duplexdiscoverer.R"

rule aggregate_interactions:
    input:
        #rds = expand(rules.run_duplexdiscoverer.output.rds, sample=samples.index),
        rds = expand("results/duplexdiscoverer/{sample}/{sample}_dgs.rds", sample=samples.index),
        gtf = config["genome_gtf"],
        installed = rules.install_duplexdiscoverer.output
    output:
        merged = "results/duplexdiscoverer/merged/snoRNA_interactions_merged.tsv",
        filtered     = "results/duplexdiscoverer/merged/snoRNA_interactions_filtered.tsv",
        superset_rds = "results/duplexdiscoverer/merged/superset_gi.rds"
    params:
        samples   = list(samples.index),
        protocols = list(samples["protocol"]),
        min_ratio = config.get("merge_min_ratio", 0.1),
        minoverlap = config.get("merge_minoverlap", 5),
        maxgap    = config.get("merge_maxgap", 50),
        min_reads   = config.get("filter_min_reads", 5)
    log:
        "results/logs/duplexdiscoverer/aggregate_interactions.log"
    conda:
        "../envs/r.yaml"
    message:
        "Aggregate interactions across samples and protocols."
    script:
        "../scripts/aggregate_interactions.R"
