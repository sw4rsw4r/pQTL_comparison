library(dplyr)
library(rtracklayer)

load_liftOver_hg19ToHg38 <- function() {
  filepath_liftOver <- "data/liftOver/hg19ToHg38.over.chain"
  if (!file.exists(filepath_liftOver)) {
    check_dir(dirname(filepath_liftOver))
    cmd1 <- "wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz -O data/liftOver/hg19ToHg38.over.chain.gz"
    cmd2 <- "gunzip data/liftOver/hg19ToHg38.over.chain.gz"
    system(cmd1)
    system(cmd2)
  }
  chain <- rtracklayer::import.chain(filepath_liftOver)
  return(chain)
}


load_CASE1 <- function(gene_of_interest, risk_factor, data_dir, window_size, data_type = "quant", case_prop = NA) {
  file_pQTL <- paste0(data_dir, "/", gene_of_interest, "_", risk_factor, ".protein_level.glm.linear")
  file_freq <- paste0(data_dir, "/", gene_of_interest, "_", risk_factor, ".afreq")

  df_pQTL <- read.delim(file_pQTL)
  df_freq <- read.delim(file_freq)
  merged_pQTL <- df_pQTL %>% left_join(df_freq, by = c("X.CHROM", "ID", "REF", "ALT"))

  gr <- with(merged_pQTL, GenomicRanges::GRanges(
    seqnames = paste0("chr", X.CHROM),
    beta = BETA,
    se = SE,
    varbeta = SE^2,
    snp = ID,
    effect = A1,
    other = ifelse(A1 == ALT, REF, ALT),
    chrom = X.CHROM,
    effect_AF = ALT_FREQS,
    type = data_type,
    pval = P,
    nlog10P = -log10(P),
    nsample = mean(OBS_CT.x),
    ranges = IRanges(start = POS - 1, end = POS)
  ))
  chain_hg19ToHg38 <- load_liftOver_hg19ToHg38()
  lifted_over <- rtracklayer::liftOver(gr, chain_hg19ToHg38)

  res <- with(as.data.frame(lifted_over), dplyr::tibble(
    beta,
    se,
    varbeta,
    snp,
    ID = paste0("chr", seqnames, "_", end, "_", other, "_", effect),
    effect,
    other,
    chrom = sub("chr", "", seqnames),
    position = end,
    effect_AF,
    type,
    s = case_prop,
    pval,
    nlog10P,
    nsample
  ))

  # remove duplicated rsIDs
  res_sorted <- res[order(abs(res$beta), decreasing = T), ]
  res_filt <- res_sorted[!duplicated(res_sorted$snp), ]
  res_final <- res_filt[order(res_filt$position), ]

  return(res_final)
}

