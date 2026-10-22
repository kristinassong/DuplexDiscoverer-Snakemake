###############################################################################
# aggregate_duplex_groups.R
#
# Merges per-sample DuplexDiscovereR duplex groups into one non-redundant
# interaction table.
#
#   1. filter each sample's DGs to those with a guide arm (snoRNA/scaRNA/tRNA)
#   2. merge overlapping DGs across samples into a superset
#   3. map each sample DG back to one superset entry, pooling read support
#   4. orient each row so the highest-priority guide biotype is on the left
#   5. write a full table and a min_reads-filtered subset
#
# Interactions with no guide arm are KEPT — their coordinates are usable
# downstream. left_type %in% guide identifies rows that do have one.
#
# snakemake wiring:
#   input:  rds (list), gtf
#   output: merged, filtered, superset_rds
#   params: samples, protocols, maxgap, minoverlap, min_reads, min_ratio
###############################################################################

library(DuplexDiscovereR)
library(InteractionSet)
library(rtracklayer)
library(GenomicRanges)
library(data.table)

rds_files <- as.character(snakemake@input[["rds"]])
sample_names <- as.character(snakemake@params[["samples"]])
protocols <- as.character(snakemake@params[["protocols"]])
maxgap <- as.integer(snakemake@params[["maxgap"]])
minoverlap <- as.integer(snakemake@params[["minoverlap"]])
min_reads <- as.integer(snakemake@params[["min_reads"]])

stopifnot(
  "samples/protocols length mismatch" =
    length(rds_files) == length(sample_names) &&
    length(rds_files) == length(protocols)
)

GeneAnnoGR <- import.gff(as.character(snakemake@input[["gtf"]]), feature.type = "gene")
colnames(mcols(GeneAnnoGR))[colnames(mcols(GeneAnnoGR)) == "gene_biotype"] <- "gene_type"

# Guide biotypes, IN PRIORITY ORDER — earlier entries win the left column.
# snoDB entries in the GTF carry gene_biotype snoRNA (2145), scaRNA (33) and
# tRNA (1); tRNA is included so that single snoDB entry is not dropped.
# NOTE: including tRNA also admits every non-snoDB tRNA gene in the GTF.
guide <- c("snoRNA", "scaRNA", "tRNA")

# Lower rank = higher priority. Non-guide biotypes AND unannotated arms (NA)
# both rank last, so one-arm-annotated and both-arms-annotated cases resolve
# through the same comparison with no special casing.
rank_of <- function(x) {
  r <- match(as.character(x), guide)
  r[is.na(r)] <- length(guide) + 1L
  r
}

# fifelse() errors when yes/no differ in type. If a paired column ever comes
# through as all-NA logical instead of character, the reshape dies mid-run;
# this makes both sides character unconditionally.
ch <- function(x) as.character(x)

# ---------------------------------------------------------------------------
# 1. load and filter per-sample DGs
# ---------------------------------------------------------------------------
gi_list <- lapply(seq_along(rds_files), function(i) {
  gi <- readRDS(rds_files[i])
  a <- mcols(gi)$gene_type.A; b <- mcols(gi)$gene_type.B
  keep <- (!is.na(a) & a %in% guide) | (!is.na(b) & b %in% guide)
  message("[", sample_names[i], "] ", sum(keep), " / ", length(gi), " DGs")
  gi[keep]
})
names(gi_list) <- sample_names

message("pooled DGs entering comparison: ",
        sum(vapply(gi_list, length, integer(1))))

# ---------------------------------------------------------------------------
# 2. build the cross-sample superset
# ---------------------------------------------------------------------------
res <- compareMultipleInteractions(
  gi_list,
  min_ratio  = as.numeric(snakemake@params[["min_ratio"]]),
  minoverlap = minoverlap,
  maxgap     = maxgap,
  anno_gr    = GeneAnnoGR
)
superset <- res$gi_all
message("superset: ", length(superset), " interactions")

# ---------------------------------------------------------------------------
# 3. map each sample DG to one superset entry, carry its read support
#
# use.region = "same" requires arm A to match arm A and arm B to match arm B.
# If a superset entry's anchors are ordered differently from a contributing
# DG, that DG does not map, its reads never reach `agg`, the entry ends at
# n_reads = 0, and slim[n_reads > 0] drops it — with no error. The mapped
# count below should be at or near 100% for every sample.
# ---------------------------------------------------------------------------
message("mapping sample DGs to superset entries")
hits <- rbindlist(lapply(seq_along(gi_list), function(i) {
  gi <- gi_list[[i]]
  h <- findOverlaps(gi, superset, minoverlap = minoverlap,
                    use.region = "same")
  h <- h[!duplicated(queryHits(h))]          # one superset entry per DG
  message(sprintf("  [%-38s] mapped %8d / %8d  (%.1f%%)",
                  sample_names[i], length(h), length(gi),
                  100 * length(h) / max(1L, length(gi))))
  if (length(h) == 0) return(NULL)
  data.table(idx      = subjectHits(h),
             sample   = sample_names[i],
             protocol = protocols[i],
             n_reads  = mcols(gi)$n_reads[queryHits(h)])
}))

stopifnot("no sample DG mapped to any superset entry" = nrow(hits) > 0)

agg <- hits[, {
  # one group-by, reused for both the label and the count. as.integer() guards
  # sprintf("%d", <double>), which errors outright if n_reads is ever numeric.
  bp <- .SD[, .(v = as.integer(sum(n_reads))), by = protocol]
  .(n_reads_total     = sum(n_reads),
    n_samples         = uniqueN(sample),
    n_protocols       = uniqueN(protocol),
    protocols         = paste(sort(unique(protocol)), collapse = ","),
    reads_by_protocol = paste(sprintf("%s:%d", bp$protocol, bp$v), collapse = ";"))
}, by = idx]

# ---------------------------------------------------------------------------
# 4. flatten and merge
#
# subjectHits() indexes into `superset`; idx is the row number of
# makeDfFromGi's output. The join is only correct if those align 1:1 in order.
# ---------------------------------------------------------------------------
out <- as.data.table(makeDfFromGi(superset))
stopifnot("makeDfFromGi row count != superset length" =
            nrow(out) == length(superset))

out[, idx := .I]
out <- merge(out, agg, by = "idx", all.x = TRUE)
out[is.na(n_reads_total),
    `:=`(n_reads_total = 0L, n_samples = 0L, n_protocols = 0L,
         protocols = NA_character_, reads_by_protocol = NA_character_)]

# annotateGI does mcols(gi) <- cbind(mcols(gi)[original_columns], dt), so
# re-annotating an already-annotated object duplicates the column names and
# makeDfFromGi's data.frame() conversion suffixes the SECOND occurrence.
# Unsuffixed = inherited per-sample annotation; ".1" = fresh annotation from
# compareMultipleInteractions. The reshape below reads ".1".
stopifnot("expected .1 annotation columns from compareMultipleInteractions" =
            all(c("gene_type.A.1", "gene_type.B.1", "gene_name.A.1",
                  "gene_name.B.1", "gene_id.A.1", "gene_id.B.1") %in% names(out)))

# If an arm carries multiple annotations the column can come through as a
# List, and match() inside rank_of() then silently matches nothing — every arm
# ranks last, left_is_A is TRUE everywhere, and the priority ordering quietly
# stops working on a complete-looking table.
message("gene_type.A.1 class: ", class(out$gene_type.A.1)[1],
        " | gene_type.B.1 class: ", class(out$gene_type.B.1)[1])
message("A   : ", paste(head(ch(out$gene_type.A), 6), collapse = " | "))
message("A.1 : ", paste(head(ch(out$gene_type.A.1), 6), collapse = " | "))

# ---------------------------------------------------------------------------
# 5. reshape: highest-priority guide biotype goes in the left column
#    snoRNA beats scaRNA beats tRNA beats everything else (including NA).
#    Ties (both arms the same rank) leave arm A on the left.
# ---------------------------------------------------------------------------
rank_A <- rank_of(out$gene_type.A.1)
rank_B <- rank_of(out$gene_type.B.1)
left_is_A <- rank_A <= rank_B

a_is_guide <- rank_A <= length(guide)
b_is_guide <- rank_B <= length(guide)

message("superset entries with no guide arm: ",
        sum(!a_is_guide & !b_is_guide), " / ", nrow(out))

slim <- data.table(
  left_name     = fifelse(left_is_A, ch(out$gene_name.A.1), ch(out$gene_name.B.1)),
  left_id       = fifelse(left_is_A, ch(out$gene_id.A.1),   ch(out$gene_id.B.1)),
  left_type     = fifelse(left_is_A, ch(out$gene_type.A.1), ch(out$gene_type.B.1)),
  left_chrom    = fifelse(left_is_A, ch(out$chromA),  ch(out$chromB)),
  left_start    = fifelse(left_is_A, out$startA,  out$startB),
  left_end      = fifelse(left_is_A, out$endA,    out$endB),
  left_strand   = fifelse(left_is_A, ch(out$strandA), ch(out$strandB)),
  right_name    = fifelse(left_is_A, ch(out$gene_name.B.1), ch(out$gene_name.A.1)),
  right_id      = fifelse(left_is_A, ch(out$gene_id.B.1),   ch(out$gene_id.A.1)),
  right_type    = fifelse(left_is_A, ch(out$gene_type.B.1), ch(out$gene_type.A.1)),
  right_chrom   = fifelse(left_is_A, ch(out$chromB), ch(out$chromA)),
  right_start   = fifelse(left_is_A, out$startB, out$startA),
  right_end     = fifelse(left_is_A, out$endB,   out$endA),
  right_strand  = fifelse(left_is_A, ch(out$strandB), ch(out$strandA)),
  n_reads       = out$n_reads_total,
  n_samples     = out$n_samples,
  n_protocols   = out$n_protocols,
  protocols     = out$protocols,
  reads_by_protocol = out$reads_by_protocol
)

slim <- slim[n_reads > 0]
setorder(slim, -n_reads)

message("interactions with read support: ", nrow(slim))
message("left arm biotypes:")
print(table(slim$left_type, useNA = "ifany"))
message("right arm biotypes (top 12):")
print(head(sort(table(slim$right_type, useNA = "ifany"), decreasing = TRUE), 12))

fwrite(slim, snakemake@output[["merged"]], sep = "\t")
saveRDS(superset, snakemake@output[["superset_rds"]])

filt <- slim[n_reads >= min_reads]
message("filtered: ", nrow(filt), " of ", nrow(slim), " (n_reads >= ", min_reads, ")")
fwrite(filt, snakemake@output[["filtered"]], sep = "\t")

message("Done")