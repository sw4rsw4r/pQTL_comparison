options(scipen = -1)
library(Rmpfr)
library(dplyr)
library(ggplot2)
library(tidyr)
# library(biomaRt)

check_dir <- function(dirpath) {
  if (!file.exists(dirpath)) {
    dir.create(dirpath, recursive = T)
  }
}

# get_gene_region_GRCh38 <- function(gene_of_interest, window_size) {
#   this_chrom <- read.delim("data/FinnGen/pQTL/Olink/probe_map.tsv") %>%
#     filter(geneName == gene_of_interest) %>%
#     pull(chr)
#   this_chrom <- sub("chr", "", this_chrom)

#   mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = "https://May2024.archive.ensembl.org/")
#   dataset <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl", mart = mart)
#   resultTable <- biomaRt::getBM(
#     attributes = c("external_gene_name", "ensembl_transcript_id", "chromosome_name", "transcript_start", "transcript_end", "gene_biotype", "description"),
#     filters = "external_gene_name",
#     values = gene_of_interest,
#     mart = dataset
#   )
#   if (nrow(resultTable) == 0) {
#     return(NULL)
#   }
#   resultTable <- resultTable %>%
#     mutate(size = transcript_end - transcript_start) %>%
#     group_by(external_gene_name, chromosome_name) %>%
#     summarise(
#       start = min(transcript_start),
#       end = max(transcript_end),
#       .groups = "drop"
#     ) %>%
#     subset(chromosome_name == this_chrom)

#   res <- with(
#     resultTable,
#     list(
#       CHR = chromosome_name,
#       locus_lower = max(start - window_size, 0),
#       locus_upper = end + window_size,
#       gene_name = gene_of_interest
#     )
#   )
#   return(res)
# }


get_gene_region_GRCh38_UKBB <- function(gene_of_interest, window_size) {
  df_anno <- vroom("data/UKBB/Metadata/Protein_annotation/olink_protein_map_3k_v1.tsv")

  if (gene_of_interest == "MYLPF") {
    gene_of_interest0 <- "MYL11"
  } else {
    gene_of_interest0 <- gene_of_interest
  }
  resultTable <- df_anno %>%
    filter(HGNC.symbol == gene_of_interest0) %>%
    head(1)

  if (gene_of_interest == "ANP32C") {
    resultTable$gene_end <- 164197711
  }
  res <- with(
    resultTable,
    list(
      CHR = chr,
      locus_lower = max(gene_start - window_size, 0),
      locus_upper = gene_end + window_size,
      gene_name = gene_of_interest
    )
  )
  return(res)
}

p_value_filtering <- function(list_data, pval_cutoff = 1e-08) {
  names_data <- names(list_data)

  snps_filtered <- list()
  for (data_type in names_data) {
    df <- list_data[[data_type]]
    snps_filtered[[data_type]] <- setdiff(with(subset(df, pval <= pval_cutoff), paste0(chrom, ":", position)), ":")
  }
  common_snps <- Reduce(intersect, snps_filtered)
  res <- length(common_snps) > 0
  return(res)
}

compute_pval <- function(beta, se, prec = 100) {
  nlog10pval <- NA
  while (any(is.na(nlog10pval)) || any(nlog10pval == Inf)) {
    nlog10pval <- asNumeric(
      round(-log10(
        2 * (1 - pnorm(abs(mpfr(beta, prec) / mpfr(se, prec))))
      ), 5)
    )
    if (any(nlog10pval == Inf)) {
      prec <- prec + 5000
    } else {
      pval <- format((2 * (1 - pnorm(abs(mpfr(beta, prec) / mpfr(se, prec))))), digits = 3, scientific = TRUE)
    }
  }
  return(list(pval = pval, nlog10pval = nlog10pval))
}

load_pQTL_deCODE <- function(gene_of_interest, risk_factor, data_dir, window_size, data_type = "quant", case_prop = NA) {
  # https://www.decode.com/summarydata/
  # Grímur Hjörleifsson Eldjarn, Egil Ferkingstad et al. Large-scale plasma proteomics comparisons through genetics and disease associations

  filepath <- list.files(data_dir, pattern = paste0("_", gene_of_interest, "_"), full.names = T)
  if (risk_factor == "deCODE_SMP") {
    filepath <- grep(".txt.gz$", grep("Proteomics_SMP_PC0", filepath, value = T), value = T)
  } else if (risk_factor == "UKBB") {
    filepath <- grep(".txt.gz$", grep("GBR_UKB_OLINK2_", filepath, value = T), value = T)
  }

  region <- get_gene_region_GRCh38(gene_of_interest, window_size)
  df <- vroom::vroom(filepath)
  df_filt0 <- df %>%
    dplyr::filter(
      Chrom == paste0("chr", region$CHR),
      Pos >= region$locus_lower,
      Pos <= region$locus_upper,
      effectAllele != otherAllele
    )
  list_multiallelic <- df_filt0$Pos[duplicated(df_filt0$Pos)]
  df_filt <- df_filt0 %>%
    dplyr::filter(
      effectAllele != "*",
      otherAllele != "*",
      nchar(effectAllele) == 1,
      nchar(otherAllele) == 1
    )
  # group_by(rsids) %>%
  # slice_max(order_by = abs(Beta), n = 1, with_ties = F) %>%
  # ungroup()
  n_total <- length(unique(df_filt0$Pos))
  n_indel <- length(unique(subset(df_filt0, !Pos %in% df_filt$Pos)$Pos))
  n_multi <- length(unique(list_multiallelic))
  ratio_indel <- n_indel / n_total
  ratio_multiallelic <- n_multi / n_total
  ratio_intersect <- length(intersect(list_multiallelic, subset(df_filt0, !Pos %in% df_filt$Pos)$Pos)) / n_total
  df_filt <- df_filt %>%
    dplyr::filter(
      !Pos %in% list_multiallelic,
      rsids != "."
    )

  computed_P <- with(df_filt, compute_pval(Beta, SE))

  res <- with(df_filt, dplyr::tibble(
    beta = Beta,
    se = SE,
    varbeta = SE^2,
    snp = sapply(strsplit(df_filt$rsids, ","), function(x) x[1]),
    effect = toupper(effectAllele),
    other = toupper(otherAllele),
    chrom = Chrom,
    position = Pos,
    MAF = ImpMAF,
    type = data_type,
    s = case_prop,
    pval = computed_P$pval,
    nlog10P = computed_P$nlog10pval,
    nsample = N,
  ))
  return(list(res = res, ratio_indel = ratio_indel, ratio_multiallelic = ratio_multiallelic, ratio_intersect = ratio_intersect))
}

download_files_FinnGen <- function(gene_of_interest, risk_factor) {
  if (risk_factor == "Olink") {
    addr <- "https://storage.googleapis.com/finngen-public-data-r10/omics/proteomics/release_2023_03_02/data/Olink/pQTL/"
    filename <- paste0("Olink_Batch1_", gene_of_interest, ".txt.gz")
  } else if (risk_factor == "Somascan") {
    addr <- "https://storage.googleapis.com/finngen-public-data-r10/omics/proteomics/release_2023_03_02/data/Somascan/pQTL/"

    probe_somascan <- "data/FinnGen/pQTL/Somascan/probe_map.tsv"
    df_probe_somascan <- read.delim(probe_somascan)
    seqid <- subset(df_probe_somascan, geneName == gene_of_interest)$AptName[1]
    filename <- paste0("SomaScan_Batch2_", seqid, ".txt.gz")
  }
  dir_download <- "data/FinnGen/pQTL/downloaded/"
  file_downloaded <- paste0(risk_factor, "_", gene_of_interest, ".txt.gz")
  cmd1 <- paste0("wget ", addr, filename, " -O ", dir_download, file_downloaded)
  cmd2 <- paste0("wget ", addr, filename, ".tbi -O ", dir_download, file_downloaded, ".tbi")
  system(cmd1)
  system(cmd2)
}


get_rsID_info <- function(region) {
  # filepath_dbSNP <- "https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz"
  filepath_dbSNP <- "/home/seongwonhwang/Desktop/projects/git/pQTL_comparison/data/dbSNP/common_all_20180418.vcf.gz"
  tabix_file <- TabixFile(filepath_dbSNP)
  region$CHR <- sub(23, "X", region$CHR)
  tabix_out <- scanTabix(tabix_file, param = with(region, GRanges(CHR, IRanges(locus_lower, locus_upper))))[[1]]
  df_dbSNP <- tabix_out %>%
    strsplit("\t") %>%
    do.call(rbind, .) %>%
    as.data.frame()
  if (nrow(df_dbSNP) == 0) {
    return(NULL)
  }
  col_names <- c("CHR", "POS", "rsID", "REF", "ALT")
  colnames(df_dbSNP)[1:5] <- col_names
  df_dbSNP <- df_dbSNP %>%
    separate_rows(ALT, sep = ",")
  df_dbSNP$POS <- as.numeric(df_dbSNP$POS)
  return(df_dbSNP[, col_names])
}

complement_allele <- function(allele) {
  comp_map <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  return(comp_map[allele])
}

read_tabix <- function(filepath, region) {
  tabix_file <- TabixFile(filepath)
  this_chrom <- sub("X", 23, region$CHR)
  tabix_out <- tryCatch(
    {
      scanTabix(tabix_file, param = with(region, GRanges(this_chrom, IRanges(locus_lower, locus_upper))))[[1]]
    },
    error = function(x) {
      return(NULL)
    }
  )
  if (length(tabix_out) == 0) {
    return(NULL)
  }
  df_pQTL <- tabix_out %>%
    strsplit("\t") %>%
    do.call(rbind, .) %>%
    as.data.frame()
  colnames(df_pQTL) <- c("CHR", "POS", "ID", "REF", "ALT", "ALT_FREQ", "BETA", "SE", "T_STAT", "P", "N")

  df_pQTL <- df_pQTL %>%
    mutate(
      ALT_FREQ = as.numeric(ALT_FREQ),
      BETA = as.numeric(BETA),
      SE = as.numeric(SE),
      P = as.numeric(P),
      N = as.numeric(N),
      POS = as.numeric(POS)
    )
  # df_dbSNP <- get_rsID_info(region)
  df_UKBBrsIDmap <- vroom(paste0("data/UKBB/Metadata/SNP_RSID_maps/olink_rsid_map_mac5_info03_b0_7_chr", region$CHR, "_patched_v2.tsv.gz")) %>%
    mutate(
      CHR = this_chrom,
      POS = POS38
    ) %>%
    dplyr::filter(CHR == region$CHR, POS >= region$locus_lower, POS <= region$locus_upper)

  df_merged <- df_pQTL %>%
    inner_join(df_UKBBrsIDmap, by = c("CHR", "POS", "REF", "ALT"))

  if (nrow(df_UKBBrsIDmap) != 0 && nrow(df_merged) == 0) {
    df_UKBBrsIDmap$REF <- sapply(df_UKBBrsIDmap$REF, complement_allele)
    df_UKBBrsIDmap$ALT <- sapply(df_UKBBrsIDmap$ALT, complement_allele)
    df_merged <- df_pQTL %>%
      inner_join(df_UKBBrsIDmap, by = c("CHR", "POS", "REF", "ALT"))
  }
  return(df_merged)
}

delete_files <- function(filepath) {
  file.remove(filepath)
  file.remove(paste0(filepath, ".tbi"))
}

load_pQTL_Finngen <- function(gene_of_interest, risk_factor, data_dir, window_size, data_type = "quant", case_prop = NA) {
  filepath <- paste0("data/FinnGen/pQTL/downloaded/", risk_factor, "_", gene_of_interest, ".txt.gz")
  if (!file.exists(filepath)) download_files_FinnGen(gene_of_interest, risk_factor)

  region <- get_gene_region_GRCh38_UKBB(gene_of_interest, window_size)
  df_pQTL <- read_tabix(filepath, region)
  delete_files(filepath)

  if (any(as.numeric(df_pQTL$P) == 0)) stop("zero p-value")

  res <- with(df_pQTL, dplyr::tibble(
    beta = BETA,
    se = SE,
    varbeta = SE^2,
    snp = rsid,
    ID = ID.x,
    effect = ALT,
    other = REF,
    chrom = CHR,
    position = POS,
    effect_AF = ALT_FREQ,
    type = data_type,
    s = case_prop,
    pval = P,
    nlog10P = -log10(P),
    nsample = N,
  ))

  # remove duplicated rsIDs
  res_sorted <- res[order(abs(res$beta), decreasing = T), ]
  res_filt <- res_sorted[!duplicated(res_sorted$snp), ]
  res_final <- res_filt[order(res_filt$position), ]

  return(res_final)
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


load_pQTL_UKBB <- function(target_id, filename_tar, gene_of_interest, window_size, risk_factor = "Olink", data_type = "quant", case_prop = NA) {
  # download files
  path_download <- "data/UKBB/UKB_PPP_pGWAS_summary_statistics"
  synGet(target_id, downloadLocation = path_download)
  list_files <- untar(file.path(path_download, filename_tar), list = T)

  region <- get_gene_region_GRCh38_UKBB(gene_of_interest, window_size)
  selected <- grep(paste0("_chr", region$CHR, "_"), list_files, value = T)
  untar(file.path(path_download, filename_tar), files = selected)

  # read the file
  # ID is based on hg19 but GENPOS is based on GRCh38
  df_pQTL <- vroom::vroom(selected) %>%
    dplyr::filter(
      GENPOS >= region$locus_lower,
      GENPOS <= region$locus_upper
    ) %>%
    mutate(CHROM = as.character(CHROM))

  df_UKBBrsIDmap <- vroom(paste0("data/UKBB/Metadata/SNP_RSID_maps/olink_rsid_map_mac5_info03_b0_7_chr", region$CHR, "_patched_v2.tsv.gz"))

  # This rsID is going to be used to match IDs in the LD matrix (https://registry.opendata.aws/ukbb-ld/)
  df_merged <- df_pQTL %>%
    inner_join(df_UKBBrsIDmap, by = "ID")

  res <- with(df_merged, dplyr::tibble(
    beta = BETA,
    se = SE,
    varbeta = SE^2,
    snp = rsid,
    ID = paste0("chr", CHROM, "_", GENPOS, "_", ALLELE0, "_", ALLELE1),
    effect = ALLELE1,
    other = ALLELE0,
    chrom = paste0("chr", CHROM),
    position = GENPOS,
    effect_AF = A1FREQ,
    type = data_type,
    s = case_prop,
    pval = 10**(-LOG10P),
    nlog10P = LOG10P,
    nsample = N,
  ))

  # remove duplicated rsIDs
  res_sorted <- res[order(abs(res$beta), decreasing = T), ]
  res_filt <- res_sorted[!duplicated(res_sorted$snp), ]
  res_final <- res_filt[order(res_filt$position), ]

  # delete files
  file.remove(file.path(path_download, filename_tar))
  unlink(dirname(selected), recursive = T)
  return(res_final)
}

load_liftOver_hg38ToHg19 <- function() {
  filepath_liftOver <- "data/liftOver/hg38ToHg19.over.chain"
  if (!file.exists(filepath_liftOver)) {
    check_dir(dirname(filepath_liftOver))
    cmd1 <- "wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -O data/liftOver/hg38ToHg19.over.chain.gz"
    cmd2 <- "gunzip data/liftOver/hg38ToHg19.over.chain.gz"
    system(cmd1)
    system(cmd2)
  }
  chain <- rtracklayer::import.chain(filepath_liftOver)
  return(chain)
}


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


get_liftOver_hg38Tohg19 <- function(chrom, position) {
  this_chrom <- ifelse(grepl("^chr", chrom), chrom, paste0("chr", chrom))
  gr <- GenomicRanges::GRanges(
    seqnames = this_chrom,
    ranges = IRanges(start = position - 1, end = position)
  )
  chain_hg38Tohg19 <- load_liftOver_hg38ToHg19()
  lifted_over <- rtracklayer::liftOver(gr, chain_hg38Tohg19)
  return(as.data.frame(lifted_over)$end)
}


get_liftOver_hg19ToHg38 <- function(chrom, position) {
  this_chrom <- ifelse(grepl("^chr", chrom), chrom, paste0("chr", chrom))
  gr <- GenomicRanges::GRanges(
    seqnames = this_chrom,
    ranges = IRanges(start = position - 1, end = position)
  )
  chain_hg19ToHg38 <- load_liftOver_hg38ToHg19()
  lifted_over <- rtracklayer::liftOver(gr, chain_hg19ToHg38)
  return(as.data.frame(lifted_over)$end)
}

get_ld_path_UKBB <- function(gene_of_interest, window_size) {
  region <- get_gene_region_GRCh38_UKBB(gene_of_interest, window_size)

  hg19_start <- with(region, get_liftOver_hg38Tohg19(CHR, locus_lower))
  hg19_end <- with(region, get_liftOver_hg38Tohg19(CHR, locus_upper))

  LD_files <- read.table("UKBB_ld_file_list.txt")$V4
  LD_files <- grep(".npz$", LD_files, value = T)
  LD_words <- strsplit(LD_files, "_")

  filename_ld <- tibble(
    filename = LD_files,
    chrom = sub("chr", "", sapply(LD_words, function(x) x[1])),
    start = as.numeric(sapply(LD_words, function(x) x[2])),
    end = as.numeric(sub(".npz", "", sapply(LD_words, function(x) x[3])))
  ) %>%
    dplyr::filter(
      chrom == region$CHR,
      start <= hg19_start
    )
  if (length(hg19_end) == 0) hg19_end <- hg19_start + 2e+06
  out <- filename_ld %>%
    dplyr::filter(
      end >= hg19_end
    ) %>%
    pull(filename) %>%
    head(1)
  return(out)
}

download_files_LD_UKBB <- function(filepath_ld) {
  dir_LD <- dirname(filepath_ld)
  fname_LD <- sub(".npz", "", basename(filepath_ld))
  cmd <- paste0("aws s3 cp s3://broad-alkesgroup-ukbb-ld/UKBB_LD/ ", dir_LD, ' --recursive --no-sign-request --exclude "*" --include "', fname_LD, '*"')
  system(cmd)
}

load_ld_mat_UKBB <- function(gene_of_interest, window_size, rsids_to_keep, dir_output) {
  dir_LD <- "data/UKBB/ld_matrix"
  filename_LD <- file.path(dir_output, "LD", "UKBB", paste0(gene_of_interest, "_window_", window_size, "_LD.RDS"))
  check_dir(dirname(filename_LD))
  if (file.exists(filename_LD)) {
    LD <- readRDS(filename_LD)
  } else {
    filename_ld <- get_ld_path_UKBB(gene_of_interest, window_size)
    filepath_ld_mat <- file.path(dir_LD, filename_ld)
    filepath_ld_meta <- sub(".npz$", ".gz", filepath_ld_mat)

    if (!file.exists(filepath_ld_mat)) download_files_LD_UKBB(filepath_ld_mat)
    ld_meta <- vroom::vroom(filepath_ld_meta, delim = "\t", comment = "#")
    ld_meta <- with(ld_meta, data.frame(snp = rsid, effect = allele2, other = allele1, idx = 1:nrow(ld_meta)))

    fname_save <- sub(".npz$", ".RDS", filepath_ld_mat)
    if (file.exists(fname_save)) {
      mat <- readRDS(fname_save)
    } else {
      np <- reticulate::import("numpy")
      npz_file <- np$load(filepath_ld_mat)
      i <- as.numeric(npz_file$f[["row"]])
      j <- as.numeric(npz_file$f[["col"]])
      v <- as.numeric(npz_file$f[["data"]])
      dims <- as.numeric(npz_file$f[["shape"]])
      mat <- Matrix::sparseMatrix(i, j, x = v, index1 = FALSE, dims = dims)
      saveRDS(mat, file = fname_save)
    }

    ld_meta_filt <- subset(ld_meta, snp %in% rsids_to_keep)
    mat_filt <- as.matrix(mat[ld_meta_filt$idx, ld_meta_filt$idx])
    mat_t <- t(mat_filt)
    lower_tri <- lower.tri(mat_filt, diag = T)
    mat_t[lower_tri] <- mat_t[lower_tri] + mat_filt[lower_tri]
    rownames(mat_t) <- colnames(mat_t) <- ld_meta_filt$snp

    LD <- list(meta = ld_meta_filt, mat = mat_t)
    saveRDS(LD, file = filename_LD)
  }
  return(LD)
}


harmonize_ld <- function(ld_mat, ld_merged) {
  flip <- with(ld_merged, ifelse(effect_1 == effect_2, 1, -1))
  ld_mat_corrected <- ld_mat * flip %o% flip
  return(ld_mat_corrected)
}

load_ld_mat_FinnGen <- function(gene_of_interest, window_size, IDs_to_keep, dir_output) {
  dir_LD <- "data/FinnGen/ld_matrix/"
  filename_LD <- file.path(dir_output, "LD", "FinnGen", paste0(gene_of_interest, "_window_", window_size, "_LD.RDS"))
  check_dir(dirname(filename_LD))
  if (file.exists(filename_LD)) {
    LD <- readRDS(filename_LD)
  } else {
    region <- get_gene_region_GRCh38_UKBB(gene_of_interest, window_size)

    filename_ld <- paste0(dir_LD, "finngen_r12_chr", region$CHR, "_ld.tsv.gz")
    tabix_file <- TabixFile(filename_ld)
    tabix_out <- scanTabix(tabix_file, param = with(region, GRanges(region$CHR, IRanges(locus_lower, locus_upper))))[[1]]

    # These IDs are based on GRCh38
    df_LD <- tabix_out %>%
      strsplit("\t") %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      setNames(c("chrom", "position", "ID1", "ID2", "R", "R2")) %>%
      mutate(
        R2 = as.numeric(R2)
      )
    if (region$CHR == "X") IDs_to_keep <- gsub("chr23", "chrX", IDs_to_keep)
    IDs_final <- Reduce(intersect, list(df_LD$ID1, df_LD$ID2, IDs_to_keep))
    df_LD_filt <- df_LD %>% dplyr::filter(ID1 %in% IDs_final, ID2 %in% IDs_final)

    df_UKBBrsIDmap <- vroom(paste0("data/UKBB/Metadata/SNP_RSID_maps/olink_rsid_map_mac5_info03_b0_7_chr", region$CHR, "_patched_v2.tsv.gz")) %>%
      mutate(
        CHR = sub(":.*", "", ID),
        POS = POS38,
        ID = paste0("chr", CHR, "_", POS, "_", REF, "_", ALT)
      ) %>%
      dplyr::filter(CHR == region$CHR, POS >= region$locus_lower, POS <= region$locus_upper)

    df_merged <- tryCatch(
      {
        df_LD_filt %>%
          inner_join(df_UKBBrsIDmap, by = c("ID1" = "ID")) %>%
          dplyr::rename(snp = rsid) %>%
          inner_join(df_UKBBrsIDmap[, c("ID", "rsid")], by = c("ID2" = "ID")) %>%
          dplyr::rename(snp2 = rsid) %>%
          dplyr::select(-ID1, -ID2) %>%
          distinct()
      },
      error = function(e) {
        return(data.frame())
      }
    )
    if (nrow(df_merged) == 0) {
      df_UKBBrsIDmap$REF <- sapply(df_UKBBrsIDmap$REF, complement_allele)
      df_UKBBrsIDmap$ALT <- sapply(df_UKBBrsIDmap$ALT, complement_allele)
      df_UKBBrsIDmap <- df_UKBBrsIDmap %>% mutate(ID = paste0("chr", CHR, "_", POS, "_", REF, "_", ALT))
      df_merged <- df_LD_filt %>%
        inner_join(df_UKBBrsIDmap, by = c("ID1" = "ID")) %>%
        dplyr::rename(snp = rsid) %>%
        inner_join(df_UKBBrsIDmap[, c("ID", "rsid")], by = c("ID2" = "ID")) %>%
        dplyr::rename(snp2 = rsid) %>%
        dplyr::select(-ID1, -ID2) %>%
        distinct()
    }

    ld_meta <- df_merged %>%
      dplyr::rename(effect = ALT, other = REF) %>%
      dplyr::select(chrom, position, snp, effect, other) %>%
      arrange(position) %>%
      distinct()

    ld_mat <- acast(df_merged, snp ~ snp2, mean, value.var = "R2", fill = 0)
    diag(ld_mat) <- 1

    LD <- list(meta = ld_meta, mat = ld_mat)
    saveRDS(LD, file = filename_LD)
  }
  return(LD)
}

harmonize <- function(runID, gene_of_interest, window_size, lst_data, LD_type, dir_output) {
  names_risk_factor <- names(lst_data$risk_factor)
  names_outcome <- names(lst_data$outcome)
  fname_harmonize <- file.path(dir_output, "harmonize", "harmonize.RDS")

  if (file.exists(fname_harmonize)) {
    res <- readRDS(fname_harmonize)
  } else {
    check_dir(dirname(fname_harmonize))
    lst_data_combined <- c(lst_data$risk_factor, lst_data$outcome)

    rsids_to_keep <- Reduce("intersect", lapply(lst_data_combined, function(x) x$snp))

    n_samples <- sapply(lst_data_combined, function(x) mean(x$nsample, na.rm = T))
    max_sample <- max(n_samples)
    selected_factor <- sort(names(which(n_samples == max_sample)))[1]

    res <- list()
    for (this_factor in setdiff(names(lst_data_combined), selected_factor)) {
      merged <- lst_data_combined[[selected_factor]] %>%
        dplyr::inner_join(lst_data_combined[[this_factor]], by = "snp", suffix = c("_1", "_2")) %>%
        dplyr::filter(snp %in% rsids_to_keep) %>%
        dplyr::mutate(
          effect_2_ori = effect_2,
          beta_2 = ifelse(other_2 == effect_1, -beta_2, beta_2),
          effect_AF_2 = ifelse(other_2 == effect_1, 1 - effect_AF_2, effect_AF_2),
          effect_2 = ifelse(other_2 == effect_1, other_2, effect_2),
          other_2 = ifelse(other_2 == effect_1, effect_2_ori, other_2)
        ) %>%
        dplyr::filter(effect_1 == effect_2)

      if (length(res) == 0) {
        res[[1]] <- merged %>% dplyr::select(snp, grep("_1", colnames(merged), value = T))
        colnames(res[[1]]) <- sub("_1$", "", colnames(res[[1]]))
      }
      dat <- merged %>% dplyr::select(snp, grep("_2", colnames(merged), value = T))
      colnames(dat) <- sub("_2$", "", colnames(dat))
      res[[length(res) + 1]] <- dat
    }
    names(res) <- c(selected_factor, setdiff(names(lst_data_combined), selected_factor))

    if (length(LD_type) == 1) {
      if (LD_type == "FinnGen") {
        IDs_to_keep <- res[[selected_factor]]$ID
        LD <- load_ld_mat_FinnGen(gene_of_interest, window_size, IDs_to_keep, dir_output)
      } else if (LD_type == "UKBB") {
        LD <- load_ld_mat_UKBB(gene_of_interest, window_size, rsids_to_keep, dir_output)
      }
      ld_merged <- as_tibble(res[[selected_factor]]) %>% left_join(LD$meta, by = "snp", suffix = c("_1", "_2"))
      ld_merged <- ld_merged %>% dplyr::filter((effect_1 == effect_2 & other_1 == other_2) | (effect_1 == other_2 & other_1 == effect_2))
      ld_mat_harmonised <- harmonize_ld(LD$mat[ld_merged$snp, ld_merged$snp], ld_merged)
      list_LD <- list()
      for (id in names(res)) list_LD[[id]] <- ld_mat_harmonised
    } else {
      IDs_to_keep <- res[[selected_factor]]$ID
      LD <- load_ld_mat_FinnGen(gene_of_interest, window_size, IDs_to_keep, dir_output)
      ld_merged <- as_tibble(res[[selected_factor]]) %>% inner_join(LD$meta, by = "snp", suffix = c("_1", "_2"))
      ld_merged <- ld_merged %>% dplyr::filter((effect_1 == effect_2 & other_1 == other_2) | (effect_1 == other_2 & other_1 == effect_2))

      LD2 <- load_ld_mat_UKBB(gene_of_interest, window_size, rsids_to_keep, dir_output)
      ld_merged2 <- as_tibble(res[[selected_factor]]) %>% inner_join(LD2$meta, by = "snp", suffix = c("_1", "_2"))
      ld_merged2 <- ld_merged2 %>% dplyr::filter((effect_1 == effect_2 & other_1 == other_2) | (effect_1 == other_2 & other_1 == effect_2))

      ld_merged <- ld_merged %>% dplyr::filter(snp %in% ld_merged2$snp)
      ld_merged2 <- ld_merged2[match(ld_merged$snp, ld_merged2$snp), ]

      list_LD <- list()
      list_LD[["FinnGen"]] <- harmonize_ld(LD$mat[ld_merged$snp, ld_merged$snp], ld_merged)
      list_LD[["UKBB"]] <- harmonize_ld(LD2$mat[ld_merged$snp, ld_merged$snp], ld_merged2)
    }

    for (id in names(res)) {
      res[[id]] <- as.list(res[[id]][match(ld_merged$snp, res[[id]]$snp), ])
      res[[id]]$type <- setdiff(unique(res[[id]]$type), NA)
      res[[id]]$N <- mean(res[[id]]$nsample)
      res[[id]]$LD <- list_LD[[id]]
      if (all(is.na(res[[id]]$s))) {
        res[[id]]$s <- NULL
      } else {
        res[[id]]$s <- unique(res[[id]]$s)
      }
    }

    res$names$risk_factors <- names_risk_factor
    res$names$outcome <- names_outcome
    saveRDS(res, file = fname_harmonize)
  }
  return(res)
}


run_colocProp <- function(res, dir_results, this_prune = 0.4, this_J = 10) {
  RF1 <- res$names$risk_factors[1]
  RF2 <- res$names$risk_factors[2]

  fname_prop.coloc.res1 <- file.path(dir_results, "propcoloc", "prop.coloc.RDS")
  fname_prop.coloc.res2 <- file.path(dir_results, "propcoloc", "prop.coloc.txt")
  check_dir(dirname(fname_prop.coloc.res1))
  # if (file.exists(fname_prop.coloc.res1) & file.exists(fname_prop.coloc.res2)) {
  #   return()
  # }

  prop.coloc.res <- NULL
  n_try <- 0
  while (is.null(prop.coloc.res)) {
    prop.coloc.res <- tryCatch(
      {
        prop.coloc::prop.coloc(
          b1 = res[[RF1]]$beta, se1 = res[[RF1]]$se, b2 = res[[RF2]]$beta, se2 = res[[RF2]]$se,
          n = c(res[[RF1]]$N, res[[RF2]]$N),
          ld = res[[RF1]]$LD, figs = TRUE, traits = c(RF1, RF2),
          prune = this_prune, J = this_J
        )
      },
      error = function(e) {
        return(NULL)
      }
    )
    if (is.null(prop.coloc.res) & n_try == 0) {
      this_prune <- 0.2
      this_J <- 5
      n_try <- n_try + 1
    } else if (is.null(prop.coloc.res) & n_try == 1) {
      prop.coloc.res <- list(p_full = NA, p_cond = NA, LM_full = NA, LM_cond = NA)
    }
  }
  prop.coloc.res$J <- this_J
  prop.coloc.res$prune <- this_prune
  prop.coloc.res$rsIDs <- res[[1]]$snp

  saveRDS(prop.coloc.res, file = fname_prop.coloc.res1)
  write.table(as.data.frame(prop.coloc.res[sapply(prop.coloc.res, length) == 1]),
    fname_prop.coloc.res2,
    quote = F, row.names = F, col.names = T, sep = "\t"
  )
}


run_coloc <- function(res, dir_results) {
  RF1 <- res$names$risk_factors[1]
  RF2 <- res$names$risk_factors[2]
  fname_output1 <- file.path(dir_results, "coloc", "coloc_sensitivity.pdf")
  fname_output2 <- file.path(dir_results, "coloc", "coloc.txt")
  check_dir(dirname(fname_output1))
  # if (file.exists(fname_output1) & file.exists(fname_output2)) {
  #   return()
  # }
  res[[RF1]]$MAF <- with(res[[RF1]], ifelse(effect_AF < .5, effect_AF, 1 - effect_AF))
  res[[RF2]]$MAF <- with(res[[RF2]], ifelse(effect_AF < .5, effect_AF, 1 - effect_AF))

  # estimated_sdY <- with(res[[RF1]], sqrt(varbeta * (N * 2 * MAF * (1 - MAF))))
  # print(estimated_sdY)

  pdf(fname_output1, width = 7, height = 5)
  coloc.res <- coloc.abf(dataset1 = res[[RF1]], dataset2 = res[[RF2]])
  sensitivity(coloc.res, rule = "H4 > 0.5")
  dev.off()

  write.table(data.frame(coloc.res$summary),
    fname_output2,
    quote = F, row.names = T, col.names = T, sep = "\t"
  )
}



run_susie <- function(res, dir_results) {
  RF1 <- res$names$risk_factors[1]
  RF2 <- res$names$risk_factors[2]

  fname_output1 <- file.path(dir_results, "susie", "susie.txt")
  fname_output2 <- file.path(dir_results, "susie", "susie_sensitivity.pdf")
  check_dir(dirname(fname_output1))
  # if (file.exists(fname_output1) & file.exists(fname_output2)) {
  #   return()
  # }
  res[[RF1]]$MAF <- with(res[[RF1]], ifelse(effect_AF < .5, effect_AF, 1 - effect_AF))
  res[[RF2]]$MAF <- with(res[[RF2]], ifelse(effect_AF < .5, effect_AF, 1 - effect_AF))

  selected <- c("beta", "varbeta", "position", "MAF", "snp", "type", "N", "LD")
  dat1 <- res[[RF1]][selected]
  dat2 <- res[[RF2]][selected]
  names(dat1$beta) <- names(dat1$varbeta) <- names(dat1$MAF) <- names(dat1$position) <- names(dat2$beta) <- names(dat2$varbeta) <- names(dat2$MAF) <- names(dat2$position) <- dat1$snp
  dat1$Z <- with(dat1, beta / sqrt(varbeta))
  dat2$Z <- with(dat2, beta / sqrt(varbeta))

  susie.res <- NULL
  tryCatch(
    {
      out <- coloc::coloc.susie(dat1, dat2)
      susie.res <- out$summary
    },
    error = function(cond) cond
  )
  if (is.null(susie.res)) {
    return()
  }
  susie.res$ori_nsnps <- length(dat1$snp)

  write.table(susie.res,
    fname_output1,
    quote = F, row.names = F, col.names = T, sep = "\t"
  )
  # pdf(fname_output2, width = 8, height = 5)
  # sensitivity(out, "H3 > 0.5 & H4 < 0.5", row = 1, dataset1 = dat1, dataset2 = dat2)
  # dev.off()
  return(susie.res)
}


run_propcoloc_Wallace <- function(res, dir_results) {
  RF1 <- res$names$risk_factors[1]
  RF2 <- res$names$risk_factors[2]

  fname_output <- file.path(dir_results, "propcoloc_Wallace", "propcoloc_Wallace.txt")
  check_dir(dirname(fname_output))
  # if (file.exists(fname_output)) {
  #   return()
  # }
  res_proptests <- run_proptests(res[[RF1]], res[[RF2]], LD = res[[RF1]]$LD)
  write.table(res_proptests,
    fname_output,
    quote = F, row.names = F, col.names = T, sep = "\t"
  )
}



get_propcoloc_res <- function(dir_results, list_factors) {
  RF1 <- list_factors[1]
  RF2 <- list_factors[2]

  fname_propcoloc <- file.path(dir_results, "propcoloc", "prop.coloc.RDS")
  res_propcoloc <- NULL
  if (file.exists(fname_propcoloc)) {
    df_tmp <- readRDS(fname_propcoloc)
    this_p_cond <- ifelse(df_tmp$p_cond == F, 1, df_tmp$p_cond)
    res_propcoloc <- with(
      df_tmp,
      data.frame(RF1 = RF1, RF2 = RF2, p_full, p_cond = this_p_cond, LM_full, LM_cond)
    )
  }

  return(res_propcoloc)
}



get_susie_res <- function(dir_results, list_factors) {
  RF1 <- list_factors[1]
  RF2 <- list_factors[2]

  fname_susie <- file.path(dir_results, "susie", "susie.txt")
  res_susie <- NULL
  if (file.exists(fname_susie)) {
    res_susie <- with(
      read.delim(fname_susie),
      data.frame(RF1, RF2, nsnps, H0 = PP.H0.abf, H1 = PP.H1.abf, H2 = PP.H2.abf, H3 = PP.H3.abf, H4 = PP.H4.abf)
    )
  }
  return(res_susie)
}

get_coloc_res <- function(dir_results, list_factors) {
  RF1 <- list_factors[1]
  RF2 <- list_factors[2]
  fname_coloc <- file.path(dir_results, "coloc", "coloc.txt")
  res_coloc <- NULL
  if (file.exists(fname_coloc)) {
    res_coloc <- with(
      as.data.frame(t(read.delim(fname_coloc))),
      data.frame(RF1 = RF1, RF2 = RF2, nsnps, H0 = PP.H0.abf, H1 = PP.H1.abf, H2 = PP.H2.abf, H3 = PP.H3.abf, H4 = PP.H4.abf)
    )
  }
  return(res_coloc)
}

get_colocPropTest_res <- function(dir_results, list_factors) {
  res_colocPropTest <- NULL
  RF1 <- list_factors[1]
  RF2 <- list_factors[2]

  fname_wallace <- file.path(dir_results, "propcoloc_Wallace", "propcoloc_Wallace.txt")
  if (file.exists(fname_wallace) && file.info(fname_wallace)$size > 1) {
    df1 <- read.delim(fname_wallace)
    res_colocPropTest <- data.frame(RF1 = RF1, RF2 = RF2, min_p = min(df1$p), min_fdr = min(df1$fdr))
  }
  return(res_colocPropTest)
}


get_all_results <- function(runID, names_risk_factor, window_size) {
  lst <- list.files(file.path("results", runID, paste0("window_", window_size)))
  merged <- NULL
  for (gene_of_interest in lst) {
    dir_results <- file.path("results", runID, paste0("window_", window_size), gene_of_interest)
    df_propcoloc <- get_propcoloc_res(dir_results, names_risk_factor)
    df_susie <- get_susie_res(dir_results, names_risk_factor)
    df_coloc <- get_coloc_res(dir_results, names_risk_factor)
    df_wallace <- get_colocPropTest_res(dir_results, names_risk_factor)

    merged <- rbind(merged, data.frame(
      gene = gene_of_interest,
      p_cond = ifelse(is.null(df_propcoloc$p_cond), NA, df_propcoloc$p_cond),
      LM_cond = ifelse(is.null(df_propcoloc$LM_cond), NA, df_propcoloc$LM_cond),
      coloc_H4 = ifelse(is.null(df_coloc$H4), NA, df_coloc$H4),
      susie_H4 = ifelse(is.null(df_susie$H4), NA, max(df_susie$H4)),
      fdr = ifelse(is.null(df_wallace$min_fdr), NA, df_wallace$min_fdr)
    ))
  }
  merged <- merged %>% mutate(
    g_coloc_H4 = factor(ifelse(coloc_H4 <= .5, "not", ">0.5"), levels = c("not", ">0.5")),
    g_susie_H4 = factor(ifelse(susie_H4 <= .5, "not", ">0.5"), levels = c("not", ">0.5")),
    g_p_cond = factor(ifelse(p_cond < .05, "sig", "not"), levels = c("not", "sig")),
    g_FDR = factor(ifelse(fdr < .05, "sig", "not"), levels = c("not", "sig"))
  )
  return(merged)
}


print_summary_tab <- function(CASE_filt, print_percent = F) {
  if (print_percent) {
    print(with(CASE_filt, round(prop.table(table(g_coloc_H4, g_susie_H4, useNA = "always")), 2)))
    print(with(CASE_filt, round(prop.table(table(g_coloc_H4, g_p_cond, LM_cond > .05, useNA = "always")), 2)))
    print(with(CASE_filt, round(prop.table(table(g_coloc_H4, g_FDR, useNA = "always")), 2)))

    print(with(CASE_filt, round(prop.table(table(g_susie_H4, g_coloc_H4, useNA = "always")), 2)))
    print(with(CASE_filt, round(prop.table(table(g_susie_H4, g_p_cond, LM_cond > .05, useNA = "always")), 2)))
    print(with(CASE_filt, round(prop.table(table(g_susie_H4, g_FDR, useNA = "always")), 2)))

    print(with(CASE_filt, round(prop.table(table(g_p_cond, g_coloc_H4, LM_cond > .05, useNA = "always")), 2)))
    print(with(CASE_filt, round(prop.table(table(g_p_cond, g_susie_H4, LM_cond > .05, useNA = "always")), 2)))
    print(with(CASE_filt, round(prop.table(table(g_p_cond, g_FDR, LM_cond > .05, useNA = "always")), 2)))

    print(with(CASE_filt, round(prop.table(table(g_FDR, g_coloc_H4, useNA = "always")), 2)))
    print(with(CASE_filt, round(prop.table(table(g_FDR, g_susie_H4, useNA = "always")), 2)))
    print(with(CASE_filt, round(prop.table(table(g_FDR, g_p_cond, LM_cond > .05, useNA = "always")), 2)))
  } else {
    print(with(CASE_filt, table(g_coloc_H4, g_susie_H4, useNA = "always")))
    print(with(subset(CASE_filt, LM_cond <= .05), table(g_coloc_H4, g_p_cond, useNA = "always")))
    print(with(subset(CASE_filt, LM_cond > .05), table(g_coloc_H4, g_p_cond, useNA = "always")))
    print(with(CASE_filt, table(g_coloc_H4, g_FDR, useNA = "always")))

    print(with(CASE_filt, table(g_susie_H4, g_coloc_H4, useNA = "always")))
    print(with(subset(CASE_filt, LM_cond <= .05), table(g_susie_H4, g_p_cond, useNA = "always")))
    print(with(subset(CASE_filt, LM_cond > .05), table(g_susie_H4, g_p_cond, useNA = "always")))
    print(with(CASE_filt, table(g_susie_H4, g_FDR, useNA = "always")))

    print(with(subset(CASE_filt, LM_cond <= .05), table(g_p_cond, g_coloc_H4, useNA = "always")))
    print(with(subset(CASE_filt, LM_cond <= .05), table(g_p_cond, g_susie_H4, useNA = "always")))
    print(with(subset(CASE_filt, LM_cond <= .05), table(g_p_cond, g_FDR, useNA = "always")))

    print(with(subset(CASE_filt, LM_cond > .05), table(g_p_cond, g_coloc_H4, useNA = "always")))
    print(with(subset(CASE_filt, LM_cond > .05), table(g_p_cond, g_susie_H4, useNA = "always")))
    print(with(subset(CASE_filt, LM_cond > .05), table(g_p_cond, g_FDR, useNA = "always")))

    print(with(CASE_filt, table(g_FDR, g_coloc_H4, useNA = "always")))
    print(with(CASE_filt, table(g_FDR, g_susie_H4, useNA = "always")))
    print(with(subset(CASE_filt, LM_cond <= .05), table(g_FDR, g_p_cond, useNA = "always")))
    print(with(subset(CASE_filt, LM_cond > .05), table(g_FDR, g_p_cond, useNA = "always")))
  }
}
