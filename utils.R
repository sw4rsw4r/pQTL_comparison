options(scipen = -1)
library(Rmpfr)
library(dplyr)
library(ggplot2)

check_dir <- function(dirpath) {
  if (!file.exists(dirpath)) {
    dir.create(dirpath, recursive = T)
  }
}

get_gene_region_GRCh38 <- function(gene, window_size = 100000) {
  # fname_genes <- "data/deCODE/2023_Large_scale_plasma/pQTL_random_seed12345_10genes_biomaRt_GRCh38.txt"
  fname_genes <- "data/deCODE/2023_Large_scale_plasma/pQTL_intersected_genes_biomaRt_GRCh38.txt"
  df_genes <- read.delim(fname_genes)
  df_genes_subset <- subset(df_genes, external_gene_name == gene)
  CHR <- df_genes_subset$chromosome_name
  locus_lower <- df_genes_subset$start
  locus_upper <- df_genes_subset$end
  gene_name <- df_genes_subset$ensembl_gene_id

  return(list(CHR = CHR, locus_lower = max(locus_lower - window_size, 0), locus_upper = locus_upper + window_size, gene_name = gene_name))
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

load_pQTL <- function(gene_of_interest, risk_factor, data_dir, window_size, data_type = "quant", case_prop = NA) {
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
  df_filt <- df %>% dplyr::filter(
    Chrom == paste0("chr", region$CHR),
    Pos >= region$locus_lower,
    Pos <= region$locus_upper,
    rsids != ".",
    effectAllele != otherAllele,
    effectAllele != "*",
    otherAllele != "*",
    nchar(effectAllele) == 1,
    nchar(otherAllele) == 1
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
