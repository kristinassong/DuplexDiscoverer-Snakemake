# Load libraries
library(DuplexDiscovereR)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(tibble)
library(dplyr)

# Patch col_check_rename bug present in DuplexDiscovereR up to v1.7.0
# Bug: n_col_check referenced in message() before assignment in STAR_COLNAMES14 branch
message("Applying col_check_rename patch (upstream bug, n_col_check used before assignment)")
assignInNamespace("col_check_rename", function(df, table_type = "STAR") {
    STAR_COLNAMES19 <- DuplexDiscovereR:::get_colnames_and_types_for_input("STAR_COLNAMES19")
    STAR_COLNAMES14 <- DuplexDiscovereR:::get_colnames_and_types_for_input("STAR_COLNAMES14")
    STAR_COLTYPES19 <- DuplexDiscovereR:::get_colnames_and_types_for_input("STAR_COLTYPES19")
    STAR_COLTYPES14 <- DuplexDiscovereR:::get_colnames_and_types_for_input("STAR_COLTYPES14")
    BEDPE_COLNAMES  <- DuplexDiscovereR:::get_colnames_and_types_for_input("BEDPE_COLNAMES")
    BEDPE_COLTYPES  <- DuplexDiscovereR:::get_colnames_and_types_for_input("BEDPE_COLTYPES")
    if (table_type != "STAR") {
        expected_names <- BEDPE_COLNAMES
        expected_types <- BEDPE_COLTYPES
        n_col_check <- 10
        message("Detected input type generic .bedpe dataframe, expecting 10+ columns\n")
    }
    if (table_type == "STAR") {
        if (all(STAR_COLNAMES19 %in% colnames(df))) {
            expected_names <- STAR_COLNAMES19
            expected_types <- STAR_COLTYPES19
            n_col_check <- 19
            message("Detected input type Chimeric.out.Junction with multimaps, expecting first 19 columns\n")
        } else {
            if (all(STAR_COLNAMES14 %in% colnames(df))) {
                expected_names <- STAR_COLNAMES14
                expected_types <- STAR_COLTYPES14
                n_col_check <- 14
                message("Detected input type Chimeric.out.Junction without multimaps, expecting first 14 columns\n")
            } else {
                expected_names <- STAR_COLNAMES14
                expected_types <- STAR_COLTYPES14
                n_col_check <- 14
                message("Specified format is STAR Chimeric.out.Junction, but columns are missing or unnamed\n")
                message("Will guess and validate first 14 columns.\n")
            }
        }
    }
    df <- as_tibble(df)
    testt <- data.frame(head(df))
    actual_names <- names(testt)
    actual_types <- vapply(testt %>% dplyr::select(all_of(seq_len(n_col_check))),
                           function(col) class(col), "")
    extra_columns <- setdiff(actual_names, expected_names)
    if (all(actual_types == expected_types)) {
        message("Verified the correct column types and order\n")
        if (!all(expected_names %in% colnames(testt))) {
            message("Guessing columns")
            df <- df %>% rename_with(~expected_names, seq_len(n_col_check))
            status <- data.frame(Before = actual_names[seq_len(n_col_check)],
                                 After  = names(df)[seq_len(n_col_check)])
            message(paste(status$Before, "\t-> ", status$After, collapse = "\n", sep = ""))
        } else {
            message("Column names are correct\n")
        }
    } else {
        result_guess <- tryCatch({
            df_tmp <- head(df, 5)
            df_tmp[seq_len(n_col_check)] <- lapply(seq_len(n_col_check),
                function(i) as(df_tmp[[i]], expected_types[i]))
            TRUE
        }, error = function(e) FALSE)
        if (result_guess) {
            df[seq_len(n_col_check)] <- lapply(seq_len(n_col_check),
                function(i) as(df[[i]], expected_types[i]))
            df <- df %>% rename_with(~expected_names,
                                     all_of(actual_names[seq_len(n_col_check)]))
            status <- data.frame(Before = actual_names[seq_len(n_col_check)],
                                 After  = names(df)[seq_len(n_col_check)])
            message(paste(status$Before, "\t-> ", status$After, collapse = "\n", sep = ""))
            message("Guessed column names and converted types\n")
        } else {
            stats <- ifelse(actual_types[seq_len(n_col_check)] == expected_types,
                            "ok", "<-problem")
            expected <- paste(seq_len(n_col_check), expected_names, ":\t",
                              expected_types, "\t||",
                              actual_names[seq_len(n_col_check)], ":\t",
                              actual_types[seq_len(n_col_check)], "\t:",
                              stats, collapse = "\n")
            message("Problem with the input format found")
            stop("\nColumn names, types, or order do not match the expected criteria.")
        }
    }
    return(df)
}, ns = "DuplexDiscovereR")
message("Patch applied successfully")

# Snakemake inputs and outputs
sj_file <- as.character(snakemake@input[["sj"]])
gtf_file <- as.character(snakemake@input[["gtf"]])
chim_file <- as.character(snakemake@input[["chim"]])
gc_file <- as.character(snakemake@input[["gene_counts"]])
dg_file <- as.character(snakemake@output[["dg"]])
sample <- as.character(snakemake@wildcards[["sample"]])
lib_type <- as.character(snakemake@params[["lib"]])
stats_file <- as.character(snakemake@output[["stats"]])
rds_file   <- as.character(snakemake@output[["rds"]])

# STAR Chimeric.out.junction
message("[", sample, "] Loading chimeric junction file")
ChimOutJunc <- as.data.frame(fread(chim_file, header = FALSE, fill = TRUE))

# drop a header line if STAR wrote one (col 2 is non-numeric on a header row)
if (!grepl("^[0-9]+$", as.character(ChimOutJunc[1, 2]))) {
  ChimOutJunc <- ChimOutJunc[-1, , drop = FALSE]
}

# name and type the first 19 columns from the package's own definition,
# so the multimap (19-column) branch of col_check_rename fires by construction
star19  <- DuplexDiscovereR:::get_colnames_and_types_for_input("STAR_COLNAMES19")
types19 <- DuplexDiscovereR:::get_colnames_and_types_for_input("STAR_COLTYPES19")
stopifnot("Expected 19+ columns" = ncol(ChimOutJunc) >= length(star19))
colnames(ChimOutJunc)[seq_along(star19)] <- star19
ChimOutJunc[seq_along(types19)] <- lapply(seq_along(types19), function(i)
  as(ChimOutJunc[[i]], types19[i]))

message("[", sample, "] Chimeric junction file dimensions: ", nrow(ChimOutJunc), " x ", ncol(ChimOutJunc))

# Verify columns
message("[", sample, "] Column names: ", paste(colnames(ChimOutJunc), collapse = ", "))
message("[", sample, "] Column classes: ", paste(sapply(ChimOutJunc, class), collapse = ", "))
message("[", sample, "] First row: ", paste(ChimOutJunc[1, ], collapse = ", "))

# Gene counts
message("[", sample, "] Loading gene counts")
GeneCounts <- fread(gc_file, header = FALSE, skip = 4)
GeneCounts <- as.data.frame(GeneCounts[, .(V1, V2)])

# Gene annotation
message("[", sample, "] Loading gene annotation")
gtf <- import.gff(gtf_file, feature.type = "gene")
GeneAnnoGR <- gtf
colnames(mcols(GeneAnnoGR))[colnames(mcols(GeneAnnoGR)) == "gene_biotype"] <- "gene_type"

# Splice junctions
message("[", sample, "] Loading splice junctions")
SpliceJnc <- fread(sj_file, header = FALSE)
setnames(SpliceJnc, c("chrom", "start", "end", "strand", "intron_motif", "annot",
                       "uniq_mapped_reads", "multi_mapped_reads", "max_overhang"))
SpliceJnc[, strand := fcase(strand == 0, "*",
                             strand == 1, "+",
                             strand == 2, "-")]
SpliceJncGR <- GRanges(
  seqnames = SpliceJnc$chrom,
  ranges   = IRanges(start = SpliceJnc$start, end = SpliceJnc$end),
  strand   = SpliceJnc$strand
)

# Run DuplexDiscovereR
message("[", sample, "] Running DuplexDiscovereR (lib_type=", lib_type, ")")
result <- DuplexDiscovereR::runDuplexDiscoverer(
  data = ChimOutJunc,
  junctions_gr = SpliceJncGR,
  anno_gr = GeneAnnoGR,
  df_counts = GeneCounts,
  sample_name = sample,
  lib_type = lib_type,
  table_type = "STAR",
  fafile = NULL,
  collapse_n_inter = as.integer(snakemake@params[["collapse_n_inter"]]),
  gap_collapse_similar = as.integer(snakemake@params[["gap_collapse_similar"]]),
  compute_p_values = as.logical(snakemake@params[["compute_p_values"]]),
  min_arm_len = as.integer(snakemake@params[["min_arm_len"]]),
  trim_alignments = as.logical(snakemake@params[["trim_alignments"]]),
  trim_length = as.integer(snakemake@params[["trim_length"]])
)

# Extract results using accessors
message("[", sample, "] Extracting results")
dg <- dd_get_duplex_groups(result)
stats <- dd_get_chimeric_reads_stats(result)
run_stats <- dd_get_run_stats(result)

# Print summary to log
message("[", sample, "] Run summary:")
print(run_stats)
message("[", sample, "] Chimeric read stats:")
print(stats)
message("[", sample, "] Total duplex groups: ", length(dg))

# Convert duplex groups to dataframe and add sample label
dg_df <- as.data.frame(makeDfFromGi(dg))
dg_df$sample <- sample
 
# Save outputs
message("[", sample, "] Saving outputs")

flatten_list_cols <- function(df) {
    df <- as.data.frame(df)
    bad <- vapply(df, function(x) is.list(x) || !is.atomic(x), logical(1))
    if (any(bad)) {
        message("Flattening non-atomic columns: ", paste(names(df)[bad], collapse = ", "))
        df[bad] <- lapply(df[bad], function(col)
            vapply(col, function(v) paste(unlist(v), collapse = ","), character(1)))
    }
    df
}

dg_df <- flatten_list_cols(dg_df)
fwrite(dg_df, file = dg_file, sep = "\t", row.names = FALSE)
fwrite(flatten_list_cols(stats), file = stats_file, sep = "\t", row.names = FALSE)
saveRDS(dg, file = rds_file)
message("[", sample, "] Done")
