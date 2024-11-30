get_gene_region_GRCh37 <- function(gene_of_interest, window_size) {
  if (gene_of_interest == "HMGCR") {
    ### target: HMGCR gene with REF: GRCh37
    # https://www.ncbi.nlm.nih.gov/gene/3156
    CHR <- "5"
    gene_start <- 74632993
    gene_end <- 74657941
    gene_name <- "ENSG00000113161"
  } else {
    stop("Gene information (CHR, gene_start, gene_end, gene_name) is required.")
  }
  return(list(CHR = CHR, locus_lower = max(gene_start - window_size, 0), locus_upper = gene_end + window_size, gene_name = gene_name))
}

get_gene_region_GRCh38 <- function(gene_of_interest, window_size) {
  if (gene_of_interest == "HMGCR") {
    CHR <- "5"
    gene_start <- 75336529
    gene_end <- 75362116
    gene_name <- "ENSG00000113161"
  } else {
    stop("Gene information (CHR, gene_start, gene_end, gene_name) is required.")
  }
  return(list(CHR = CHR, locus_lower = max(gene_start - window_size, 0), locus_upper = gene_end + window_size, gene_name = gene_name))
}


load_T2D <- function(filepath_T2D, gene_of_interest, window_size) {
  region <- get_gene_region_GRCh37(gene_of_interest, window_size)

  N <- 251739.50 # 51.15% This is accurate
  case_prop <- 80154 / (80154 + 853816)

  df <- vroom::vroom(filepath_T2D)
  res <- with(df, dplyr::tibble(
    beta = ifelse(effect_allele_frequency < .5, `Fixed-effects_beta`, -`Fixed-effects_beta`),
    se = `Fixed-effects_SE`,
    varbeta = `Fixed-effects_SE`^2,
    snp = rsID,
    effect = ifelse(effect_allele_frequency < .5, toupper(effect_allele), toupper(other_allele)),
    other = ifelse(effect_allele_frequency < .5, toupper(other_allele), toupper(effect_allele)),
    chrom = `chromosome(b37)`,
    position = `position(b37)`,
    MAF = ifelse(effect_allele_frequency < .5, effect_allele_frequency, 1 - effect_allele_frequency),
    type = "cc",
    s = case_prop,
    pval = `Fixed-effects_p-value`,
    nsample = N
  )) %>% dplyr::filter(!is.na(snp), chrom == region$CHR, position >= region$locus_lower, position <= region$locus_upper)
  return(res)
}

check_dir <- function(dirpath) {
  if (!file.exists(dirpath)) {
    dir.create(dirpath, recursive = T)
  }
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

get_rsID_info <- function(gene_of_interest, window_size) {
  fname_tabix_input <- paste0("data/GWAS/", gene_of_interest, "_window_", window_size, "_region_for_tabix_input.txt")
  fname_tabix_output <- paste0("data/GWAS/", gene_of_interest, "_window_", window_size, "_region_for_tabix_output.txt")

  if (!file.exists(fname_tabix_output)) {
    region <- get_gene_region_GRCh38(gene_of_interest, window_size)
    df_tabix_input <- with(region, data.frame(CHR, locus_lower, locus_upper))
    write.table(df_tabix_input, fname_tabix_input, quote = F, row.names = F, col.names = F, sep = "\t")

    filepath_dbSNP <- "data/dbSNP/common_all_20180418.vcf.gz"
    if (!file.exists(filepath_dbSNP)) {
      check_dir(dirname(filepath_dbSNP))
      cmd1 <- "wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz -O data/dbSNP/common_all_20180418.vcf.gz"
      cmd2 <- "wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz.tbi -O data/dbSNP/common_all_20180418.vcf.gz.tbi"
      system(cmd1)
      system(cmd2)
    }
    cmd <- paste0("tabix ", filepath_dbSNP, " -R ", fname_tabix_input, " > ", fname_tabix_output)
    system(cmd)
  }

  df_tabix_output <- vroom(fname_tabix_output, col_names = F) %>%
    mutate(CHR = X1, POS = X2, rsID = X3, REF = X4, ALT = X5) %>%
    filter(nchar(REF) == 1, nchar(ALT) == 1) %>%
    select(CHR, POS, rsID)

  return(df_tabix_output)
}

load_BMI <- function(gene_of_interest, window_size) {
  filepath_BMI <- "data/GWAS/fat-distn.giant.ukbb.meta-analysis.bmi.combined.tbl.gz"

  if (!file.exists(filepath_BMI)) {
    # http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009004/readme.txt
    check_dir(dirname(filepath_BMI))
    cmd <- "wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009004/fat-distn.giant.ukbb.meta-analysis.bmi.combined.tbl.gz -O data/GWAS/fat-distn.giant.ukbb.meta-analysis.bmi.combined.tbl.gz"
    system(cmd)
  }
  nsample <- 806834
  df <- vroom::vroom(filepath_BMI)
  df_tabix_output <- get_rsID_info(gene_of_interest, window_size)
  rsIDs <- sapply(strsplit(df$SNP, ":"), function(x) x[1])
  df_filt <- df %>%
    mutate(rsID = rsIDs) %>%
    inner_join(df_tabix_output, by = "rsID")

  gr <- with(df_filt, GenomicRanges::GRanges(
    seqnames = paste0("chr", CHR),
    beta = ifelse(Freq1 < .5, Effect, -Effect),
    se = StdErr,
    snp = rsID,
    effect = ifelse(Freq1 < .5, toupper(Allele1), toupper(Allele2)),
    other = ifelse(Freq1 < .5, toupper(Allele2), toupper(Allele1)),
    MAF = ifelse(Freq1 < .5, Freq1, 1 - Freq1),
    pval = `P-value`,
    ranges = IRanges(start = POS - 1, end = POS)
  ))
  chain_hg38Tohg19 <- load_liftOver_hg38ToHg19()
  lifted_over <- rtracklayer::liftOver(gr, chain_hg38Tohg19)

  res <- with(as.data.frame(lifted_over), dplyr::tibble(
    beta,
    se,
    varbeta = se^2,
    snp,
    effect = effect,
    other = other,
    chrom = sub("chr", "", seqnames),
    position = end,
    MAF,
    type = "quant",
    s = NA,
    pval,
    nsample
  ))
  return(res)
}


get_rsID_map <- function(gene_of_interest, window_size) {
  fname_rsID <- paste0("data/dbSNP/", gene_of_interest, "_window_", window_size, "_rsID_map.RDS")
  if (file.exists(fname_rsID)) {
    ID_map <- readRDS(fname_rsID)
  } else {
    df_tabix_output <- get_rsID_info(gene_of_interest, window_size)

    gr <- with(df_tabix_output, GenomicRanges::GRanges(
      seqnames = paste0("chr", CHR),
      snp = rsID,
      ranges = IRanges(start = POS - 1, end = POS)
    ))
    chain_hg38Tohg19 <- load_liftOver_hg38ToHg19()
    lifted_over <- rtracklayer::liftOver(gr, chain_hg38Tohg19)
    ID_map <- list()
    for (idx in 1:length(lifted_over)) {
      CHR <- sub("chr", "", as.data.frame(lifted_over[[idx]])$seqnames)
      POS <- as.data.frame(lifted_over[[idx]])$end
      rsID <- as.data.frame(lifted_over[[idx]])$snp
      ID_map[[paste0(CHR, ":", POS)]] <- rsID
    }
    saveRDS(ID_map, file = fname_rsID)
  }
  return(ID_map)
}

load_GWAS <- function(pheno, gene_of_interest, window_size, data_type = "quant", case_prop = NA) {
  check_dir("data/GWAS/")

  # https://www.ebi.ac.uk/gwas/
  if (pheno == "Acute Insulin response") {
    fname <- "data/GWAS/28490609-GCST004575-EFO_0006831.h.tsv.gz"
    download_addr <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004575/harmonised/28490609-GCST004575-EFO_0006831.h.tsv.gz"
    nsample <- 4765
  }
  if (pheno == "Fasting Insulin") {
    fname <- "data/GWAS/22581228-GCST005185-EFO_0004466.h.tsv.gz"
    download_addr <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005185/harmonised/22581228-GCST005185-EFO_0004466.h.tsv.gz"
    nsample <- 51750
  }
  if (pheno == "Fasting Glucose") {
    fname <- "data/GWAS/22581228-GCST005186-EFO_0004465.h.tsv.gz"
    download_addr <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005186/harmonised/22581228-GCST005186-EFO_0004465.h.tsv.gz"
    nsample <- 58074
  }
  if (pheno == "Triglyceride") {
    fname <- "data/GWAS/GCST90239664_buildGRCh37.tsv"
    download_addr <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90239001-GCST90240000/GCST90239664/GCST90239664_buildGRCh37.tsv"
    nsample <- 1320016
  }
  if (pheno == "Leptin") {
    fname <- "data/GWAS/32917775-GCST90007310-EFO_0005000-Build37.f.tsv.gz"
    download_addr <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90007001-GCST90008000/GCST90007310/harmonised/32917775-GCST90007310-EFO_0005000-Build37.f.tsv.gz"
    nsample <- 49909
  }
  if (pheno == "Sterol") {
    fname <- "data/GWAS/34503513-GCST90060133-EFO_0010231-Build37.f.tsv.gz"
    download_addr <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90060001-GCST90061000/GCST90060133/harmonised/34503513-GCST90060133-EFO_0010231-Build37.f.tsv.gz"
    nsample <- 13814
  }
  if (pheno == "Cortisol") {
    fname <- "data/GWAS/GCST90200378_buildGRCh38.tsv.gz"
    download_addr <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90200001-GCST90201000/GCST90200378/GCST90200378_buildGRCh38.tsv.gz"
    nsample <- 8193
  }
  if (pheno == "Estradiol") {
    fname <- "data/GWAS/34255042-GCST90020092-EFO_0004697-Build38.f.tsv.gz"
    download_addr <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90020001-GCST90021000/GCST90020092/harmonised/34255042-GCST90020092-EFO_0004697-Build38.f.tsv.gz"
    nsample <- 163985
    data_type <- "cc"
    case_prop <- 37461 / (37461 + 126524)
  }
  if (pheno == "Vitamin D") {
    fname <- "data/GWAS/GCST90162562_buildGRCh37.tsv"
    download_addr <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90162001-GCST90163000/GCST90162562/GCST90162562_buildGRCh37.tsv"
    nsample <- 64988
  }
  if (pheno == "Bile acid") {
    fname <- "data/GWAS/34503513-GCST90060135-EFO_0010231-Build37.f.tsv.gz"
    download_addr <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90060001-GCST90061000/GCST90060135/harmonised/34503513-GCST90060135-EFO_0010231-Build37.f.tsv.gz"
    nsample <- 13814
  }
  if (pheno == "Aldosterone") {
    fname <- "data/GWAS/GCST90012609_buildGRCh37.tsv.gz"
    download_addr <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90012001-GCST90013000/GCST90012609/GCST90012609_buildGRCh37.tsv.gz"
    nsample <- 1128
  }
  if (pheno == "CRP") {
    fname <- "data/GWAS/GCST90309897.h.tsv.gz"
    download_addr <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90309001-GCST90310000/GCST90309897/harmonised/GCST90309897.h.tsv.gz"
    nsample <- 174488
  }
  if (pheno == "Ubiquinone") {
    fname <- "data/GWAS/35668104-GCST90024608-EFO_0021486-Build37.f.tsv.gz"
    download_addr <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90024001-GCST90025000/GCST90024608/harmonised/35668104-GCST90024608-EFO_0021486-Build37.f.tsv.gz"
    nsample <- 4492
  }
  if (pheno == "CAD") {
    fname <- "data/GWAS/GCST90132314.h.tsv.gz"
    download_addr <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132314/harmonised/GCST90132314.h.tsv.gz"
    nsample <- 1165690
    data_type <- "cc"
    case_prop <- 181522 / (181522 + 984168)
  }
  if (pheno == "Testosterone") {
    fname <- "data/GWAS/GCST90319620.h.tsv.gz"
    download_addr <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90319001-GCST90320000/GCST90319620/harmonised/GCST90319620.h.tsv.gz"
    nsample <- 197921
  }
  if (pheno == "LDL-C") {
    fname <- "data/GWAS/GCST90239658_buildGRCh37.tsv"
    download_addr <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90239001-GCST90240000/GCST90239658/GCST90239658_buildGRCh37.tsv"
    nsample <- 1320016
  }
  if (pheno == "HDL-C") {
    fname <- "data/GWAS/GCST90239652_buildGRCh37.tsv"
    download_addr <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90239001-GCST90240000/GCST90239652/GCST90239652_buildGRCh37.tsv"
    nsample <- 1320016
  }

  if (!file.exists(fname)) {
    cmd1 <- paste0("wget ", download_addr, " -O ", fname)
    system(cmd1)
  }
  if (pheno == "HDL-C") {
    # correct rsID using LDL-C data
    cmd2 <- paste0("Rscript correct_HDLC_rsID.R")
    system(cmd2)
    fname <- "data/GWAS/GCST90239652_buildGRCh37_corrected_rsID.tsv"
  }

  df <- vroom(fname)
  genomic_version <- ifelse(grepl("GRCh37", fname) || grepl("[Bb]uild37", fname), "GRCh37", "GRCh38")

  if (fname == "data/GWAS/34255042-GCST90020092-EFO_0004697-Build38.f.tsv.gz") genomic_version <- "GRCh37"

  if (genomic_version == "GRCh38") {
    region <- get_gene_region_GRCh38(gene_of_interest, window_size)
    df_filt <- df %>% dplyr::filter(
      chromosome == region$CHR,
      base_pair_location >= region$locus_lower,
      base_pair_location <= region$locus_upper
    )
    if (fname == "data/GWAS/GCST90309897.h.tsv.gz") df_filt <- df_filt %>% mutate(variant_id = rsid)
    if (fname == "data/GWAS/GCST90319620.h.tsv.gz") df_filt <- df_filt %>% mutate(variant_id = rsid)

    gr <- with(df_filt, GenomicRanges::GRanges(
      seqnames = paste0("chr", chromosome),
      beta = ifelse(effect_allele_frequency < .5, beta, -beta),
      se = standard_error,
      snp = variant_id,
      effect = ifelse(effect_allele_frequency < .5, effect_allele, other_allele),
      other = ifelse(effect_allele_frequency < .5, other_allele, effect_allele),
      MAF = ifelse(effect_allele_frequency < .5, effect_allele_frequency, 1 - effect_allele_frequency),
      pval = p_value,
      ranges = IRanges(start = base_pair_location - 1, end = base_pair_location)
    ))

    chain_hg38Tohg19 <- load_liftOver_hg38ToHg19()
    lifted_over <- rtracklayer::liftOver(gr, chain_hg38Tohg19)

    res <- with(as.data.frame(lifted_over), dplyr::tibble(
      beta,
      se,
      varbeta = se^2,
      snp,
      effect = toupper(effect),
      other = toupper(other),
      chrom = sub("chr", "", seqnames),
      position = end,
      MAF,
      type = data_type,
      s = case_prop,
      pval,
      posID = paste0(chrom, ":", position),
      nsample = nsample
    ))
  } else if (genomic_version == "GRCh37") {
    region <- get_gene_region_GRCh37(gene_of_interest, window_size)
    df_filt <- df %>% dplyr::filter(
      chromosome == region$CHR,
      base_pair_location >= region$locus_lower,
      base_pair_location <= region$locus_upper
    )
    if (fname == "data/GWAS/GCST90012609_buildGRCh37.tsv.gz") {
      df_filt <- df_filt %>% rename(BETA = "beta", SE = "standard_error", A1 = "effect_allele", A2 = "other_allele", AF1 = "effect_allele_frequency")
    }
    if (fname == "data/GWAS/GCST90162562_buildGRCh37.tsv") {
      df_filt <- df_filt %>% rename(SNP = "variant_id")
    }
    res <- with(df_filt, dplyr::tibble(
      beta = ifelse(effect_allele_frequency < 0.5, beta, -beta),
      se = standard_error,
      varbeta = standard_error^2,
      snp = variant_id,
      effect = toupper(ifelse(effect_allele_frequency < 0.5, effect_allele, other_allele)),
      other = toupper(ifelse(effect_allele_frequency < 0.5, other_allele, effect_allele)),
      chrom = chromosome,
      position = base_pair_location,
      MAF = ifelse(effect_allele_frequency < 0.5, effect_allele_frequency, 1 - effect_allele_frequency),
      type = data_type,
      s = case_prop,
      pval = p_value,
      posID = paste0(chrom, ":", position),
      nsample = nsample
    ))
  }

  ID_map <- get_rsID_map(gene_of_interest, window_size)
  rsIDs_from_dbSNP <- sapply(
    res$posID,
    function(x) ifelse(is.null(ID_map[[x]]), NA, ID_map[[x]])
  )
  res <- res %>%
    dplyr::mutate(snp = ifelse(is.na(snp), rsIDs_from_dbSNP, snp)) %>%
    dplyr::filter(!is.na(beta) & !is.na(effect) & !is.na(MAF)) %>%
    distinct()
  return(res)
}


load_ld_mat_UKBB <- function(filepath_ld_mat, ld_merged) {
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

  filt <- sort(ld_merged$idx)
  mat_filt <- as.matrix(mat[filt, filt])
  mat_t <- t(mat_filt)
  lower_tri <- lower.tri(mat_filt, diag = T)
  mat_t[lower_tri] <- mat_t[lower_tri] + mat_filt[lower_tri]
  rownames(mat_t) <- colnames(mat_t) <- ld_merged$snp
  return(mat_t)
}


harmonize_ld <- function(ld_mat, ld_merged) {
  flip <- with(ld_merged, ifelse(effect_1 == effect_2, 1, -1))
  ld_mat_corrected <- ld_mat * flip %o% flip
  return(ld_mat_corrected)
}

harmonize <- function(gene_of_interest, window_size, lst_data, filepath_ld_mat, filepath_ld_meta, dir_output) {
  names_risk_factor <- names(lst_data$risk_factor)
  names_outcome <- names(lst_data$outcome)
  ID <- paste(sort(unique(c(names_risk_factor, names_outcome))), collapse = "_")
  fname_harmonize <- paste0(dir_output, "/harmonize/", ID, "/harmonize_", ID, ".RDS")


  if (file.exists(fname_harmonize)) {
    res <- readRDS(fname_harmonize)
  } else {
    check_dir(dirname(fname_harmonize))
    lst_data_combined <- c(lst_data$risk_factor, lst_data$outcome)

    # if (any(sapply(lst_data_combined, nrow) < 2)) {
    #   return(NULL)
    # }

    ld_meta <- vroom::vroom(filepath_ld_meta, delim = "\t", comment = "#")
    ld_meta <- with(ld_meta, data.frame(snp = rsid, effect = allele1, other = allele2, idx = 1:nrow(ld_meta)))

    ids_to_keep <- Reduce("intersect", lapply(lst_data_combined, function(x) x$snp))

    n_samples <- sapply(lst_data_combined, function(x) mean(x$nsample))
    max_sample <- max(n_samples)
    selected_factor <- sort(names(which(n_samples == max_sample)))[1]

    res <- list()
    ids_do_not_match <- NULL
    for (this_factor in setdiff(names(lst_data_combined), selected_factor)) {
      merged <- lst_data_combined[[selected_factor]] %>%
        dplyr::inner_join(lst_data_combined[[this_factor]], by = "snp", suffix = c("_1", "_2")) %>%
        dplyr::filter(snp %in% ids_to_keep) %>%
        dplyr::mutate(
          effect_2_ori = effect_2,
          beta_2 = ifelse(other_2 == effect_1, -beta_2, beta_2),
          MAF_2 = ifelse(other_2 == effect_1, 1 - MAF_2, MAF_2),
          effect_2 = ifelse(other_2 == effect_1, other_2, effect_2),
          other_2 = ifelse(other_2 == effect_1, effect_2_ori, other_2)
        )

      ids_do_not_match <- unique(c(ids_do_not_match, merged %>% dplyr::filter(effect_1 != effect_2) %>% pull(snp)))

      if (length(res) == 0) {
        res[[1]] <- merged %>% dplyr::select(snp, grep("_1", colnames(merged), value = T))
        colnames(res[[1]]) <- sub("_1$", "", colnames(res[[1]]))
      }
      dat <- merged %>% dplyr::select(snp, grep("_2", colnames(merged), value = T))
      colnames(dat) <- sub("_2$", "", colnames(dat))
      res[[length(res) + 1]] <- dat
    }
    names(res) <- c(selected_factor, setdiff(names(lst_data_combined), selected_factor))

    ld_merged <- as_tibble(res[[selected_factor]]) %>% left_join(ld_meta, by = "snp", suffix = c("_1", "_2"))
    ld_merged <- ld_merged %>% dplyr::filter(!is.na(effect_2) & (effect_1 == effect_2 | effect_1 == other_2))
    ids_do_not_match_with_ld <- setdiff(ids_to_keep, ld_merged$snp)

    ids_to_remove <- unique(c(ids_do_not_match_with_ld, ids_do_not_match))

    ld_merged <- ld_merged %>% dplyr::filter(!snp %in% ids_to_remove)
    ld_merged <- ld_merged[!duplicated(ld_merged$snp), ]
    for (id in names(res)) {
      res[[id]] <- res[[id]] %>%
        dplyr::filter(!duplicated(snp)) %>%
        as_tibble() %>%
        dplyr::filter(!snp %in% ids_to_remove)
    }

    message(paste0("SNPs removed: ", paste(ids_to_remove, collapse = ",")))

    # read and correct LD
    ld_mat <- load_ld_mat_UKBB(filepath_ld_mat, ld_merged)
    res$ld <- harmonize_ld(ld_mat, ld_merged)

    SNPs_final <- colnames(res$ld)

    for (id in names(res)) {
      if (id == "ld") {
        res$ld <- res$ld[SNPs_final, SNPs_final]
      } else {
        res[[id]] <- res[[id]] %>%
          dplyr::filter(snp %in% SNPs_final) %>%
          as.list()
        res[[id]]$type <- unique(res[[id]]$type)
        res[[id]]$N <- mean(res[[id]]$nsample)
        if (all(is.na(res[[id]]$s))) {
          res[[id]]$s <- NULL
        } else {
          res[[id]]$s <- unique(res[[id]]$s)
        }
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

  dir_output <- file.path(dir_results, "propcoloc", paste(res$names$risk_factors, collapse = "_"))
  check_dir(dir_output)
  fname_prop.coloc.res1 <- paste0(dir_output, "/03_prop.coloc_", RF1, "_", RF2, ".RDS")
  fname_prop.coloc.res2 <- paste0(dir_output, "/03_prop.coloc_", RF1, "_", RF2, ".txt")
  if (file.exists(fname_prop.coloc.res1) & file.exists(fname_prop.coloc.res2)) {
    return()
  }

  lst_lipids <- c("LDL-C", "HDL-C", "Triglyceride")

  prop.coloc.res <- NULL
  n_try <- 0
  while (is.null(prop.coloc.res)) {
    if (RF1 %in% lst_lipids & RF2 %in% lst_lipids) {
      prop.coloc.res <- tryCatch(
        {
          prop.coloc::prop.coloc(
            b1 = res[[RF1]]$beta, se1 = res[[RF1]]$se, b2 = res[[RF2]]$beta, se2 = res[[RF2]]$se,
            n = res[[RF1]]$N,
            ld = res$ld, figs = TRUE, traits = c(RF1, RF2),
            prune = this_prune, J = this_J
          )
        },
        error = function(e) {
          return(NULL)
        }
      )
    } else {
      prop.coloc.res <- tryCatch(
        {
          prop.coloc::prop.coloc(
            b1 = res[[RF1]]$beta, se1 = res[[RF1]]$se, b2 = res[[RF2]]$beta, se2 = res[[RF2]]$se,
            n = c(res[[RF1]]$N, res[[RF2]]$N),
            ld = res$ld, figs = TRUE, traits = c(RF1, RF2),
            prune = this_prune, J = this_J
          )
        },
        error = function(e) {
          return(NULL)
        }
      )
    }
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
  dir_output <- file.path(dir_results, "coloc", paste(res$names$risk_factors, collapse = "_"))
  check_dir(dir_output)
  fname_output1 <- paste0(dir_output, "/01_coloc_sensitivity_", RF1, "_", RF2, ".pdf")
  fname_output2 <- paste0(dir_output, "/01_coloc_", RF1, "_", RF2, ".txt")
  if (file.exists(fname_output1) & file.exists(fname_output2)) {
    return()
  }

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
  dir_output <- file.path(dir_results, "susie", paste(res$names$risk_factors, collapse = "_"))
  check_dir(dir_output)
  fname_output1 <- file.path(dir_output, "susie.txt")
  fname_output2 <- file.path(dir_output, "susie_sensitivity.pdf")
  if (fname_output1 == 'results/window_10000/HMGCR/susie/LDL-C_T2D/susie.txt') {
    return(read.delim(fname_output1))
  } else{
    return()
  }

  selected <- c("beta", "varbeta", "position", "MAF", "snp", "type", "N")
  dat1 <- res[[RF1]][selected]
  dat2 <- res[[RF2]][selected]
  dat1$LD <- dat2$LD <- res$ld
  names(dat1$beta) <- names(dat1$varbeta) <- names(dat1$MAF) <- names(dat1$position) <- names(dat2$beta) <- names(dat2$varbeta) <- names(dat2$MAF) <- names(dat2$position) <- colnames(res$ld) <- rownames(res$ld) <- dat1$snp
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
  dir_output <- file.path(dir_results, "propcoloc_Wallace", paste(res$names$risk_factors, collapse = "_"))
  check_dir(dir_output)
  fname_output <- paste0(dir_output, "/03_propcoloc_Wallace_", RF1, "_", RF2, ".txt")
  if (file.exists(fname_output)) {
    return()
  }

  res <- run_proptests(res[[RF1]], res[[RF2]], LD = res$ld)
  write.table(res,
    fname_output,
    quote = F, row.names = F, col.names = T, sep = "\t"
  )
}


get_propcoloc_res <- function(dir_results, list_factors) {
  res <- NULL
  for (idx1 in 1:length(list_factors)) {
    for (idx2 in idx1:length(list_factors)) {
      fac1 <- list_factors[idx1]
      fac2 <- list_factors[idx2]
      this_dir <- file.path(dir_results, "propcoloc", paste0(fac1, "_", fac2))

      fname_propcoloc <- paste0(this_dir, "/03_prop.coloc_", fac1, "_", fac2, ".RDS")
      if (!file.exists(fname_propcoloc)) next

      fac1_ <- paste(unlist(strsplit(fac1, " ")), collapse = "\n")
      fac2_ <- paste(unlist(strsplit(fac2, " ")), collapse = "\n")
      df_tmp <- readRDS(fname_propcoloc)
      this_p_cond <- ifelse(df_tmp$p_cond == F, 1, df_tmp$p_cond)
      df <- with(
        df_tmp,
        data.frame(group = paste0(fac1, "_", fac2), group0 = paste0(fac2, "_", fac1), fac1 = fac1_, fac2 = fac2_, p_full, p_cond = this_p_cond, LM_full, LM_cond)
      )
      res <- rbind(res, df)
    }
  }
  res_filtered <- res %>% filter(p_cond < .05 & LM_cond < .05)
  return(list(full = res, filtered = res_filtered))
}

get_susie_res <- function(dir_results, list_factors) {
  res <- NULL
  for (idx1 in 1:length(list_factors)) {
    for (idx2 in idx1:length(list_factors)) {
      fac1 <- list_factors[idx1]
      fac2 <- list_factors[idx2]
      this_dir <- file.path(dir_results, "susie", paste0(fac1, "_", fac2))

      fname_susie <- file.path(this_dir, "susie.txt")
      if (file.exists(fname_susie)) {
        fac1_ <- paste(unlist(strsplit(fac1, " ")), collapse = "\n")
        fac2_ <- paste(unlist(strsplit(fac2, " ")), collapse = "\n")
        df <- with(
          read.delim(fname_susie)[1, ],
          data.frame(group = paste0(fac1, "_", fac2), group0 = paste0(fac2, "_", fac1), fac1 = fac1_, fac2 = fac2_, nsnps, H0 = PP.H0.abf, H1 = PP.H1.abf, H2 = PP.H2.abf, H3 = PP.H3.abf, H4 = PP.H4.abf)
        )
        res <- rbind(res, df)
      }
    }
  }
  res_filtered <- res %>% filter(H3 >= .5 & H4 <= .5)
  return(list(full = res, filtered = res_filtered))
}

get_coloc_res <- function(dir_results, list_factors) {
  res <- NULL
  for (idx1 in 1:length(list_factors)) {
    for (idx2 in idx1:length(list_factors)) {
      fac1 <- list_factors[idx1]
      fac2 <- list_factors[idx2]
      this_dir <- file.path(dir_results, "coloc", paste0(fac1, "_", fac2))
      fname_coloc <- paste0(this_dir, "/01_coloc_", fac1, "_", fac2, ".txt")
      if (file.exists(fname_coloc)) {
        fac1_ <- paste(unlist(strsplit(fac1, " ")), collapse = "\n")
        fac2_ <- paste(unlist(strsplit(fac2, " ")), collapse = "\n")
        df <- with(
          as.data.frame(t(read.delim(fname_coloc))),
          data.frame(group = paste0(fac1, "_", fac2), group0 = paste0(fac2, "_", fac1), fac1 = fac1_, fac2 = fac2_, nsnps, H0 = PP.H0.abf, H1 = PP.H1.abf, H2 = PP.H2.abf, H3 = PP.H3.abf, H4 = PP.H4.abf)
        )
        res <- rbind(res, df)
      }
    }
  }
  res_filtered <- res %>% filter(H3 >= .5 & H4 <= .5)
  return(list(full = res, filtered = res_filtered))
}

get_colocPropTest_res <- function(dir_results, list_factors) {
  res <- NULL
  for (idx1 in 1:length(list_factors)) {
    for (idx2 in idx1:length(list_factors)) {
      fac1 <- list_factors[idx1]
      fac2 <- list_factors[idx2]
      this_dir <- file.path(dir_results, "propcoloc_Wallace", paste0(fac1, "_", fac2))

      fname_wallace <- paste0(this_dir, "/03_propcoloc_Wallace_", fac1, "_", fac2, ".txt")
      if (file.exists(fname_wallace)) {
        fac1_ <- paste(unlist(strsplit(fac1, " ")), collapse = "\n")
        fac2_ <- paste(unlist(strsplit(fac2, " ")), collapse = "\n")
        df1 <- read.delim(fname_wallace)
        df2 <- data.frame(group = paste0(fac1, "_", fac2), group0 = paste0(fac2, "_", fac1), fac1 = fac1_, fac2 = fac2_, min_p = min(df1$p), min_fdr = min(df1$fdr))
        res <- rbind(res, df2)
      }
    }
  }
  res_filtered <- res %>% filter(min_fdr < .05)
  return(list(full = res, filtered = res_filtered))
}

plot_propcoloc_and_susie_barplots_pairwise <- function(gene_of_interest, list_factors, dir_results, fig_name) {
  df_propcoloc <- get_propcoloc_res(dir_results, list_factors)
  df_susie <- get_susie_res(dir_results, list_factors)
  df_coloc <- get_coloc_res(dir_results, list_factors)

  p_propcoloc <- df_propcoloc$full %>%
    mutate(
      GROUP = group,
      RF1 = fac1,
      RF2 = fac2,
      val1 = ifelse(-log10(p_cond) > 10, 10, -log10(p_cond)),
      val2 = ifelse(-log10(LM_cond) > 10, 10, -log10(LM_cond))
    ) %>%
    select(GROUP, RF1, RF2, "val1", "val2") %>%
    reshape2::melt() %>%
    mutate(
      Color = factor(ifelse(value < -log10(0.05), "Not significant", as.character(variable)),
        levels = c("val1", "val2", "Not significant"),
        labels = c("Proportionality P < 0.05", "LM P < 0.05", "Not significant")
      ),
      position = "top"
    )

  p_susie <- df_susie$full %>%
    mutate(GROUP = group0, RF1 = fac2, RF2 = fac1, val1 = H3, val2 = H4) %>%
    select(GROUP, RF1, RF2, "val1", "val2") %>%
    reshape2::melt() %>%
    mutate(
      value = value * 10,
      Color = factor(
        ifelse(value < 5, "Not significant", as.character(variable)),
        levels = c("val1", "val2", "Not significant"),
        labels = c("H3 > 0.5", "H4 > 0.5", "Not significant")
      ),
      position = "bottom"
    )

  p_coloc <- df_coloc$full %>%
    filter(group %in% setdiff(group, df_susie$full$group)) %>%
    mutate(GROUP = group0, RF1 = fac2, RF2 = fac1, val1 = H3, val2 = H4) %>%
    select(GROUP, RF1, RF2, "val1", "val2") %>%
    reshape2::melt() %>%
    mutate(
      value = value * 10,
      Color = factor(
        ifelse(value < 5, "Not significant", as.character(variable)),
        levels = c("val1", "val2", "Not significant"),
        labels = c("H3 > 0.5", "H4 > 0.5", "Not significant")
      ),
      position = "bottom"
    )

  df_p <- Reduce("rbind", list(p_propcoloc, p_susie, p_coloc))

  levels_RF <- levels(factor(df_p$RF1))
  p <- easyGgplot2::ggplot2.barplot(
    data = df_p,
    xName = "variable", yName = "value",
    facetingVarNames = c("RF1", "RF2"),
    groupName = "Color",
    faceting = TRUE
  ) +
    scale_fill_manual(values = c(`Proportionality P < 0.05` = "#00c5c5", `LM P < 0.05` = "#19bb04", `Not significant` = "#888888", `H3 > 0.5` = "#d49409", `H4 > 0.5` = "#cb4cd6")) +
    theme_bw() +
    geom_hline(aes(yintercept = ifelse(position == "bottom", 5, -log10(0.05))), color = "red", lty = 2) +
    ggtitle(gene_of_interest) + ylab("") + xlab("") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    geom_rect(
      aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
      data = data.frame(RF1 = levels_RF, RF2 = levels_RF, variable = "val1", value = 0),
      fill = "#bdbdbd"
    ) +
    geom_rect(
      data = (p_susie),
      fill = "gold", alpha = .1, colour = NA, xmin = -Inf, xmax = Inf,
      ymin = -Inf, ymax = Inf
    ) +
    geom_rect(
      data = (p_propcoloc %>% filter(GROUP %in% df_propcoloc$filtered$group)),
      fill = NA, colour = "black", xmin = -Inf, xmax = Inf,
      lwd = 2, ymin = -Inf, ymax = Inf
    ) +
    geom_rect(
      data = (p_susie %>% filter(GROUP %in% df_susie$filtered$group0)),
      fill = NA, colour = "black", xmin = -Inf, xmax = Inf,
      lwd = 2, ymin = -Inf, ymax = Inf
    ) +
    geom_rect(
      data = (p_coloc %>% filter(GROUP %in% df_coloc$filtered$group0)),
      fill = NA, colour = "black", xmin = -Inf, xmax = Inf,
      lwd = 2, ymin = -Inf, ymax = Inf
    ) +
    geom_rect(
      data = (p_propcoloc %>% filter(GROUP %in% intersect(df_propcoloc$filtered$group, c(df_susie$filtered$group, df_coloc$filtered$group)))),
      fill = NA, colour = "red", xmin = -Inf, xmax = Inf,
      lwd = 2, ymin = -Inf, ymax = Inf
    ) +
    geom_rect(
      data = (p_susie %>% filter(GROUP %in% intersect(df_susie$filtered$group0, df_propcoloc$filtered$group0))),
      fill = NA, colour = "red", xmin = -Inf, xmax = Inf,
      lwd = 2, ymin = -Inf, ymax = Inf
    ) +
    geom_rect(
      data = (p_coloc %>% filter(GROUP %in% intersect(df_coloc$filtered$group0, df_propcoloc$filtered$group0))),
      fill = NA, colour = "red", xmin = -Inf, xmax = Inf,
      lwd = 2, ymin = -Inf, ymax = Inf
    )
  this_width <- max(5, length(list_factors) * 1.5)
  this_height <- max(5, length(list_factors) * 1.5) - 1.5
  pdf(fig_name, width = this_width, height = this_height)
  plot(p)
  dev.off()
  plot(p)
}


plot_propcoloc_and_wallace_barplots_pairwise <- function(gene_of_interest, list_factors, dir_results, fig_name) {
  df_propcoloc <- get_propcoloc_res(dir_results, list_factors)
  df_wallace <- get_colocPropTest_res(dir_results, list_factors)

  p_propcoloc <- df_propcoloc$full %>%
    mutate(
      GROUP = group,
      RF1 = fac1,
      RF2 = fac2,
      val1 = ifelse(-log10(p_cond) > 10, 10, -log10(p_cond)),
      val2 = ifelse(-log10(LM_cond) > 10, 10, -log10(LM_cond))
    ) %>%
    select(GROUP, RF1, RF2, "val1", "val2") %>%
    reshape2::melt() %>%
    mutate(
      Color = factor(ifelse(value < -log10(0.05), "Not significant", as.character(variable)),
        levels = c("val1", "val2", "Not significant"),
        labels = c("Proportionality P < 0.05", "LM P < 0.05", "Not significant")
      )
    )

  p_wallace <- df_wallace$full %>%
    mutate(GROUP = group0, RF1 = fac2, RF2 = fac1, val1 = min_p, val2 = min_fdr) %>%
    select(GROUP, RF1, RF2, "val1", "val2") %>%
    reshape2::melt() %>%
    mutate(
      value = -log10(value),
      Color = factor(
        ifelse(value < -log10(0.05), "Not significant", as.character(variable)),
        levels = c("val1", "val2", "Not significant"),
        labels = c("Uncorrected P < 0.05", "FDR < 0.05", "Not significant")
      )
    )

  df_p <- Reduce("rbind", list(p_propcoloc, p_wallace))

  levels_RF <- levels(factor(df_p$RF1))
  p <- easyGgplot2::ggplot2.barplot(
    data = df_p,
    xName = "variable", yName = "value",
    facetingVarNames = c("RF1", "RF2"),
    groupName = "Color",
    faceting = TRUE
  ) +
    scale_fill_manual(values = c(`Proportionality P < 0.05` = "#00c5c5", `LM P < 0.05` = "#19bb04", `Not significant` = "#888888", `Uncorrected P < 0.05` = "#d49409", `FDR < 0.05` = "#cb4cd6")) +
    theme_bw() +
    geom_hline(yintercept = -log10(0.05), color = "red", lty = 2) +
    ggtitle(gene_of_interest) + ylab("") + xlab("") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    geom_rect(
      aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
      data = data.frame(RF1 = levels_RF, RF2 = levels_RF, variable = "val1", value = 0),
      fill = "#bdbdbd"
    ) +
    geom_rect(
      data = (p_propcoloc %>% filter(GROUP %in% df_propcoloc$filtered$group)),
      fill = NA, colour = "black", xmin = -Inf, xmax = Inf,
      lwd = 2, ymin = -Inf, ymax = Inf
    ) +
    geom_rect(
      data = (p_wallace %>% filter(GROUP %in% df_wallace$filtered$group0)),
      fill = NA, colour = "black", xmin = -Inf, xmax = Inf,
      lwd = 2, ymin = -Inf, ymax = Inf
    ) +
    geom_rect(
      data = (p_propcoloc %>% filter(GROUP %in% intersect(df_propcoloc$filtered$group, df_wallace$filtered$group))),
      fill = NA, colour = "red", xmin = -Inf, xmax = Inf,
      lwd = 2, ymin = -Inf, ymax = Inf
    ) +
    geom_rect(
      data = (p_wallace %>% filter(GROUP %in% intersect(df_wallace$filtered$group0, df_propcoloc$filtered$group0))),
      fill = NA, colour = "red", xmin = -Inf, xmax = Inf,
      lwd = 2, ymin = -Inf, ymax = Inf
    )

  this_width <- max(5, length(list_factors) * 1.5)
  this_height <- max(5, length(list_factors) * 1.5) - 1.5
  pdf(fig_name, width = this_width, height = this_height)
  plot(p)
  dev.off()
  plot(p)
}


select_phenotypic_heterogeneity <- function(dir_results, list_factors) {
  df_propcoloc <- get_propcoloc_res(dir_results, list_factors)
  df_susie <- get_susie_res(dir_results, list_factors)
  df_coloc <- get_coloc_res(dir_results, list_factors)
  df_wallace <- get_colocPropTest_res(dir_results, list_factors)

  return(Reduce("intersect", list(
    df_propcoloc$filtered$group, df_wallace$filtered$group,
    c(df_susie$filtered$group, df_coloc$filtered$group)
  )))
}


get_table_from_pca_res <- function(pca_res, n_risk_factors) {
  ci <- with(pca_res, data.frame(
    estimate = liml,
    lower = liml - qnorm(1 - 0.05 / 2) * se.liml,
    upper = liml + qnorm(1 - 0.05 / 2) * se.liml,
    `p-value` = 2 * (1 - pnorm(abs(liml / se.liml))),
    PCs = rep(factors, n_risk_factors),
    `OID test p-value` = rep(1 - pchisq(Q, factors - (n_risk_factors + 1)), n_risk_factors)
  ))
  return(round(ci, 3))
}

compute_Fstat <- function(MRInputObj, nx, ny, n_PCs, multivariate) {
  if (multivariate) {
    mr_pcgmm_Fstat <- tryCatch(
      {
        mr_mvpcgmm(MRInputObj, nx = nx, ny = ny, r = n_PCs, robust = T)@CondFstat
      },
      error = function(e) {
        return(NULL)
      }
    )
  } else {
    mr_pcgmm_Fstat <- tryCatch(
      {
        mr_pcgmm(MRInputObj, nx = nx, ny = ny, r = n_PCs, robust = T)@Fstat
      },
      error = function(e) {
        return(NULL)
      }
    )
  }
  return(mr_pcgmm_Fstat)
}

pca.no <- function(res, thres) {
  bx <- sapply(res[res$names$risk_factors], function(x) x$beta)
  sx <- sapply(res[res$names$risk_factors], function(x) x$se)
  nx <- sapply(res[res$names$risk_factors], function(x) mean(x$nsample))

  a <- 1 / ((nx[1] * sx[, 1]^2) + bx[, 1]^2)
  A <- (sqrt(a) %*% t(sqrt(a))) * res$ld
  return(which(cumsum(prcomp(A, scale = FALSE)$sdev^2 / sum((prcomp(A, scale = FALSE)$sdev^2))) > thres)[1])
}

run_PCA_liml <- function(res, dir_results, n_PCs = NULL, cor.x = NULL) {
  names_risk_factors <- res$names$risk_factors
  names_outcome <- res$names$outcome
  ld <- res$ld

  dir_output <- file.path(dir_results, "MVMR", names_outcome, paste(sort(names_risk_factor), collapse = "_"))
  check_dir(dir_output)

  bx <- sapply(res[names_risk_factors], function(x) x$beta)
  sx <- sapply(res[names_risk_factors], function(x) x$se)
  nx <- sapply(res[names_risk_factors], function(x) mean(x$nsample))

  by <- sapply(res[names_outcome], function(x) x$beta)
  sy <- sapply(res[names_outcome], function(x) x$se)
  ny <- sapply(res[names_outcome], function(x) mean(x$nsample))

  if (is.null(n_PCs)) {
    thres <- 0.99
    n_PCs <- pca.no(res, thres)
    fname_output <- paste0("MVMR_PCA_liml_thres_", thres, ".txt")
  } else {
    fname_output <- paste0("MVMR_PCA_liml_nPCs_", n_PCs, ".txt")
  }

  # unconditional correlation between exposures set to 0
  if (is.null(cor.x)) cor.x <- diag(length(names_risk_factors))
  get_pca_res <- function(n_PCs) {
    PCA_liml(
      bx = bx, sx = sx,
      by = by, sy = sy,
      rho0 = ld,
      cor.x = cor.x,
      nx = nx, ny = ny,
      r = n_PCs
    )
  }

  if (ncol(bx) == 1) {
    multivariate <- F
    MRInputObj <- MendelianRandomization::mr_input(bx = as.vector(bx), bxse = as.vector(sx), by = as.vector(by), byse = as.vector(sy), correlation = ld)
  } else {
    multivariate <- T
    MRInputObj <- MendelianRandomization::mr_mvinput(bx = bx, bxse = sx, by = as.vector(by), byse = as.vector(sy), correlation = ld)
  }

  pca_res <- tryCatch(
    {
      get_pca_res(n_PCs)
    },
    error = function(e) {
      return(NULL)
    }
  )

  if (!is.null(pca_res)) {
    # robust PCA-GMM results using principal components explaining 99.99% of genetic variation
    res_ci <- get_table_from_pca_res(pca_res, n_risk_factors = length(names_risk_factors))
    res_ci$Fstat <- compute_Fstat(MRInputObj, nx, ny, n_PCs, multivariate)
    res_ci$factor <- names_risk_factors
    rownames(res_ci) <- names_risk_factors
  } else {
    stop("Something wrong in PCA_liml.")
  }

  write.table(res_ci,
    file.path(dir_output, fname_output),
    row.names = F, quote = F, col.names = T, sep = "\t"
  )
  return(res_ci)
}


forestplot_OR <- function(list_of_outcomes, dir_results, fig_name) {
  p_df <- NULL
  for (name_outcome in list_of_outcomes) {
    DIR <- file.path(dir_results, "MVMR", name_outcome)
    list_files <- list.files(DIR, full.names = T)
    list_files <- list_files[base::sapply(list_files, function(x) file.info(x)$isdir)]

    for (this_file in list_files) {
      this_pair <- basename(this_file)
      names_risk_factors <- sort(unlist(strsplit(this_pair, "_")))

      if (name_outcome %in% names_risk_factors) next
      filename_MVMR <- file.path(this_file, "MVMR_PCA_liml_thres_0.99.txt")
      if (!file.exists(filename_MVMR)) next

      df <- read.delim(filename_MVMR) %>%
        dplyr::mutate(
          outcome = name_outcome,
          estimate = as.numeric(estimate),
          lower = as.numeric(lower),
          upper = as.numeric(upper),
          p.value = as.numeric(p.value),
          nlogP = -log10(p.value),
          Fstat = ifelse(Fstat > 10, 10, Fstat),
          Direction = factor(ifelse(estimate > 0, "Pos", "Neg"), levels = c("Pos", "Neg")),
          Pairs = this_pair,
          factor = factor
        )
      p_df <- rbind(p_df, df)
    }
  }
  p_df_filt <- p_df %>%
    mutate(
      sig = p.value < 0.05,
      col = ifelse(p.value < .05, "sig", "Not significant")
    ) %>%
    group_by(outcome, Pairs) %>%
    filter(sum(sig) > 0) %>%
    group_by(Pairs) %>%
    filter(length(unique(factor[sig])) == 2)

  p1 <-
    ggplot(p_df_filt, aes(x = factor, y = exp(estimate), colour = col)) +
    geom_hline(yintercept = 1, lty = 2, col = "#bdb9b9") +
    geom_point(shape = 15, size = 1.25, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymin = exp(lower), ymax = exp(upper)), position = position_dodge(0.5), width = 0.2, size = 0.5) +
    theme_classic(base_size = 16) +
    ggh4x::facet_grid2(Pairs ~ outcome, scales = "free") +
    scale_color_manual(values = c(`Not significant` = "black", sig = "#aa0202")) +
    theme(legend.position = "none", strip.text.y = element_blank(), text = element_text(face = "bold")) +
    xlab("") +
    ylab("Estimate (OR)") +
    scale_y_continuous(trans = "log10") +
    coord_flip()

  pdf(fig_name, width = 7.5, height = 3.5)
  plot(p1)
  dev.off()
  # write.table(p_df_filt, sub(".pdf$", ".txt", fig_name), quote = F, row.names = F, col.names = T, sep = "\t")
  plot(p1)
}


calc_lambda <- function(dat, ld, n_PCs, normalise = T) {
  a <- with(dat, 1 / ((mean(nsample) * se^2) + beta^2))
  A <- (sqrt(a) %*% t(sqrt(a))) * ld
  lambda <- sqrt(nrow(ld)) * prcomp(A, scale = FALSE)$rotation[, 1:n_PCs]
  if (normalise) {
    evec <- eigen((t(lambda) %*% lambda))$vectors
    eval <- eigen((t(lambda) %*% lambda))$values
    lambda <- lambda %*% (solve(evec %*% diag(sqrt(eval)) %*% t(evec)))
  }
  return(lambda)
}


# transform to PCA-standardised multivariable betas
pca.transform <- function(beta, se, n, ld, lambda, n_PCs) {
  a <- 1 / ((n * se^2) + beta^2)
  A <- (sqrt(a) %*% t(sqrt(a))) * ld
  B <- a * beta
  A.f <- t(lambda) %*% A %*% lambda
  B.f <- as.vector(t(lambda) %*% B)
  pca.var <- solve(A.f) * (1 - as.numeric(t(B.f) %*% solve(A.f) %*% B.f))
  pca.var <- pca.var * (1 / (n - n_PCs + 1))
  inv.pca.var <- solve(pca.var)
  evec <- eigen(inv.pca.var)$vectors
  eval <- eigen(inv.pca.var)$values
  inv.pca.var.sq <- evec %*% diag(sqrt(eval)) %*% t(evec)
  pca.beta <- inv.pca.var.sq %*% as.vector(solve(A.f) %*% B.f)
  return(pca.beta)
}

get_comb_names <- function(names_risk_factor) {
  res_idx <- res_names <- NULL
  n_RF <- length(names_risk_factor)
  lst_idx <- sapply(1:n_RF, function(x) combn(1:n_RF, x, simplify = F))
  lst_names <- sapply(1:n_RF, function(x) combn(names_risk_factor, x, simplify = F))
  for (n_fact in 1:length(lst_names)) {
    comb_idx <- lst_idx[[n_fact]]
    comb_names <- lst_names[[n_fact]]
    for (idx in 1:length(comb_names)) {
      res_idx <- c(res_idx, paste(comb_idx[[idx]], collapse = ","))
      res_names <- c(res_names, paste(comb_names[[idx]], collapse = ","))
    }
  }
  return(list(idx = res_idx, names = res_names))
}

run_BMA <- function(res, dir_results) {
  names_risk_factor <- res$names$risk_factor
  names_outcome <- res$names$outcome[1]
  ld <- res$ld

  dir_output <- file.path(dir_results, "BMA", names_outcome, paste(sort(names_risk_factor), collapse = "_"))
  check_dir(dir_output)

  bx <- sapply(res[names_risk_factor], function(x) x$beta)
  sx <- sapply(res[names_risk_factor], function(x) x$se)
  nx <- sapply(res[names_risk_factor], function(x) mean(x$nsample))

  n_PCs <- pca.no(res, thres = 0.99)
  lambda <- calc_lambda(res[[1]], ld, n_PCs, normalise = F)

  pca.bx <- sapply(res[names_risk_factor], function(x) {
    pca.transform(x$beta, x$se, mean(x$nsample), ld, lambda, n_PCs)
  })
  pca.by <- with(
    res[[names_outcome]],
    pca.transform(beta, se, mean(nsample), ld, lambda, n_PCs)
  )

  amd_nmr_input <- new("mvMRInput",
    betaX = pca.bx,
    betaY = pca.by,
    snps = rownames(ld),
    exposure = colnames(pca.bx),
    outcome = names_outcome
  )
  BMA_output <- summarymvMR_SSS(amd_nmr_input,
    kmin = 1,
    kmax = 10,
    prior_prob = 0.1,
    max_iter = 100
  )

  out_comb <- get_comb_names(names_risk_factor)
  names(BMA_output@pp) <- factor(names(BMA_output@pp), levels = out_comb$idx, labels = out_comb$names)
  names(BMA_output@pp_marginal) <- names_risk_factor
  BMA_pp <- sort(BMA_output@pp, decreasing = TRUE)
  BMA_marginal <- sort(BMA_output@pp_marginal, decreasing = TRUE)

  write.table(data.frame(BMA_pp), file.path(dir_output, "03_BMA_pp.txt"), quote = F, row.names = T, col.names = T, sep = "\t")
  write.table(data.frame(BMA_marginal), file.path(dir_output, "03_BMA_pp_marginal.txt"), quote = F, row.names = T, col.names = T, sep = "\t")

  df_BMA <- tibble(BMA_pp, group = names(BMA_pp)) %>%
    arrange(desc(BMA_pp)) %>%
    dplyr::mutate(group = factor(group, levels = group)) %>%
    head(20)

  p <- ggplot(df_BMA, aes(x = group, y = BMA_pp)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(names_outcome)
  pdf(file.path(dir_output, "03_BMA_pp.pdf"), width = 8, height = 4.5)
  plot(p)
  dev.off()

  df_BMA_pp_marginal <- tibble(group = names((BMA_output@pp_marginal)), pp_marginal = BMA_output@pp_marginal) %>%
    arrange(desc(pp_marginal)) %>%
    dplyr::mutate(group = factor(group, levels = group))
  p1 <- ggplot(df_BMA_pp_marginal, aes(x = group, y = pp_marginal)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(names_outcome)
  pdf(file.path(dir_output, "03_BMA_pp_marginal.pdf"), width = 8, height = 4.5)
  plot(p1)
  dev.off()
  return(data.frame(BMA_pp))
}
