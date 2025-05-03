options(scipen = -1)
suppressMessages(library(Rmpfr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(vroom))
suppressMessages(library(rtracklayer))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(curl))
suppressMessages(library(MendelianRandomization))
suppressMessages(library(coloc))
suppressMessages(library(susieR))
suppressMessages(library(Rsamtools))
suppressMessages(library(stringr))
suppressMessages(library(combinat))
suppressMessages(library(hash))
suppressMessages(library(jsonlite))

ensure_dir <- function(dirpath) {
  if (!file.exists(dirpath)) {
    dir.create(dirpath, recursive = T)
  }
}

get_gene_region_GRCh38_UKBB <- function(gene_of_interest, WINDOW_SIZE) {
  # This file was downloaded from the following URL: https://www.synapse.org/Synapse:syn51396728
  df_anno <- vroom("data/UKBB/Metadata/Protein_annotation/olink_protein_map_3k_v1.tsv", show_col_types = FALSE)

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
      locus_lower = max(gene_start - WINDOW_SIZE, 0),
      locus_upper = gene_end + WINDOW_SIZE,
      gene_name = gene_of_interest
    )
  )
  return(res)
}

download_metadata <- function(){
  # Modify your config.json file to include your Synapse authentication token
  # The config.json file should look like this:
  # {
  #   "synapse_token": "your_token_here"
  # }
  # Load the configuration from config.json file
  config <- fromJSON("config.json")
  # Retrieve the Synapse authentication token from the configuration file
  auth_token <- config$synapse_token
  # Login to Synapse using the token
  synLogin(authToken = auth_token)

  path_SNP_RSID_maps = "data/UKBB/Metadata/SNP_RSID_maps"
  path_Protein_annotation = "data/UKBB/Metadata/Protein_annotation"
  # Download the rsID mapping files from Synapse
  synapserutils::syncFromSynapse("syn51396727", path = path_SNP_RSID_maps)
  # Download the protein annotation file from Synapse
  synapserutils::syncFromSynapse("syn51396728", path = path_Protein_annotation)
}

download_files_FinnGen <- function(gene_of_interest, risk_factor, dir_download) {
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
  file_downloaded <- paste0(risk_factor, "_", gene_of_interest, ".txt.gz")
  cmd1 <- paste0("wget ", addr, filename, " -O ", dir_download, file_downloaded)
  cmd2 <- paste0("wget ", addr, filename, ".tbi -O ", dir_download, file_downloaded, ".tbi")
  system(cmd1)
  system(cmd2)
}


# complement_allele <- function(allele) {
#   comp_map <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
#   return(comp_map[allele])
# }

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
  df_UKBBrsIDmap <- vroom(paste0("data/UKBB/Metadata/SNP_RSID_maps/olink_rsid_map_mac5_info03_b0_7_chr", region$CHR, "_patched_v2.tsv.gz")) %>%
    mutate(
      CHR = this_chrom,
      POS = POS38
    ) %>%
    dplyr::filter(CHR == region$CHR, POS >= region$locus_lower, POS <= region$locus_upper)

  df_merged <- df_pQTL %>%
    inner_join(df_UKBBrsIDmap, by = c("CHR", "POS", "REF", "ALT"))

  # if (nrow(df_UKBBrsIDmap) != 0 && nrow(df_merged) == 0) {
  #   df_UKBBrsIDmap$REF <- sapply(df_UKBBrsIDmap$REF, complement_allele)
  #   df_UKBBrsIDmap$ALT <- sapply(df_UKBBrsIDmap$ALT, complement_allele)
  #   df_merged <- df_pQTL %>%
  #     inner_join(df_UKBBrsIDmap, by = c("CHR", "POS", "REF", "ALT"))
  # }
  return(df_merged)
}

delete_files <- function(filepath) {
  file.remove(filepath)
  file.remove(paste0(filepath, ".tbi"))
}

load_pQTL_Finngen <- function(gene_of_interest, risk_factor, WINDOW_SIZE) {
  dir_download <- "data/FinnGen/pQTL/downloaded/"
  ensure_dir(dir_download)
  filepath <- paste0(dir_download, "/", risk_factor, "_", gene_of_interest, ".txt.gz")
  if (!file.exists(filepath)) download_files_FinnGen(gene_of_interest, risk_factor, dir_download)

  region <- get_gene_region_GRCh38_UKBB(gene_of_interest, WINDOW_SIZE)
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
    type = "quant",
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


load_pQTL_UKBB <- function(target_id, filename_tar, gene_of_interest, WINDOW_SIZE) {
  # download files
  path_download <- "data/UKBB/UKB_PPP_pGWAS_summary_statistics"
  synGet(target_id, downloadLocation = path_download)
  list_files <- untar(file.path(path_download, filename_tar), list = T)

  region <- get_gene_region_GRCh38_UKBB(gene_of_interest, WINDOW_SIZE)
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
    type = "quant",
    pval = 10**(-LOG10P),
    nlog10P = LOG10P,
    nsample = N
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

get_variant_data_with_retry <- function(batch, genomes = "GRCh38", max_retries = 5, wait_time = 5) {
  retries <- 0

  while (retries < max_retries) {
    tryCatch(
      {
        batch_result <- getVariantPopData(varids = batch, genomes = genomes)
        batch_result <- Filter(Negate(is.null), batch_result)
        if (length(batch_result) == 0) {
          return(NULL)
        }
        if (is.data.frame(batch_result[[1]])) {
          return(batch_result)
        } else {
          stop("Retry")
        }
      },
      error = function(e) {
        if (grepl("Retry", e$message)) {
          retries <- retries + 1
          message(sprintf("Retrying in %d seconds... (Attempt %d/%d)", wait_time, retries, max_retries))
          Sys.sleep(wait_time)
        } else {
          stop(e)
        }
      }
    )
  }
  stop("Failed after maximum retries")
}

add_AF <- function(res) {
  snp_list <- res %>%
    mutate(varid = paste0(chrom, "-", as.integer(position), "-", other, "-", effect)) %>%
    pull(varid)

  batch_size <- 10
  snp_batches <- split(snp_list, ceiling(seq_along(snp_list) / batch_size))

  result_list <- list()
  for (batch in snp_batches) {
    batch_result <- get_variant_data_with_retry(batch)
    temp <- lapply(batch_result, function(df) {
      if (is.null(df)) {
        return(data.frame())
      } else {
        return(with(subset(df, id == "nfe"), data.frame(varid, chrom, pos, AF = ac / an)))
      }
    })
    result_list <- c(result_list, temp)
  }
  final_df <- bind_rows(result_list, .id = "varid")
  final_df$varid <- NULL
  res$effect_AF <- NULL
  colnames(final_df) <- c("chrom", "position", "effect_AF")
  merged <- merge(res, final_df, by = c("chrom", "position"))
  return(merged)
}


load_pQTL_EGA <- function(gene_of_interest, filepath, WINDOW_SIZE) {
  region <- get_gene_region_GRCh38_UKBB(gene_of_interest, WINDOW_SIZE)
  filename <- paste0("data/EGA/", gene_of_interest, ".txt")
  if (!file.exists(filename)) {
    cmd <- paste0("wget ", filepath, " -O ", filename)
    system(cmd)
  }

  df_pQTL <- vroom::vroom(filename, delim = "\t") %>%
    filter(
      chromosome == region$CHR,
      base_pair_location > region$locus_lower,
      base_pair_location < region$locus_upper
    )
  if (nrow(df_pQTL) == 0) {
    return(NULL)
  }

  cols_reqd <- c("beta", "standard_error", "variant_id", "effect_allele", "other_allele", "chromosome", "base_pair_location", "effect_allele_frequency", "p_value", "log(P)")

  if (all(cols_reqd %in% names(df_pQTL))) {
    df_pQTL_subset <- df_pQTL[, cols_reqd]
    if (all(is.na(df_pQTL$variant_id))) df_pQTL_subset$variant_id <- df_pQTL$rsid

    res <- with(df_pQTL_subset, dplyr::tibble(
      beta = beta,
      se = standard_error,
      varbeta = standard_error^2,
      snp = variant_id,
      ID = paste0("chr", chromosome, "_", base_pair_location, "_", toupper(other_allele), "_", toupper(effect_allele)),
      effect = toupper(effect_allele),
      other = toupper(other_allele),
      chrom = chromosome,
      position = base_pair_location,
      effect_AF = effect_allele_frequency,
      type = "quant",
      pval = p_value,
      nlog10P = `log(P)`,
      nsample = 3301
    ))
  }
  file.remove(filename)
  if (all(is.na(res$effect_AF))) res <- add_AF(res)

  # remove duplicated rsIDs
  res_sorted <- res[order(abs(res$beta), decreasing = T), ]
  res_filt <- res_sorted[!duplicated(res_sorted$snp), ]
  res_final <- res_filt[order(res_filt$position), ]

  return(res_final)
}

check_significant_hits <- function(dataset, gene, platform) {
  nrow(subset(dataset[[platform]], pval < pval_cutoff)) > 0
}


load_liftOver_hg38ToHg19 <- function() {
  filepath_liftOver <- "data/liftOver/hg38ToHg19.over.chain"
  if (!file.exists(filepath_liftOver)) {
    ensure_dir(dirname(filepath_liftOver))
    cmd1 <- "wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -O data/liftOver/hg38ToHg19.over.chain.gz"
    cmd2 <- "gunzip data/liftOver/hg38ToHg19.over.chain.gz"
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

get_ld_path_UKBB <- function(gene_of_interest, WINDOW_SIZE) {
  region <- get_gene_region_GRCh38_UKBB(gene_of_interest, WINDOW_SIZE)

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

load_ld_mat_UKBB <- function(gene_of_interest, WINDOW_SIZE, rsids_to_keep, dir_output) {
  dir_LD <- "data/UKBB/ld_matrix"
  filename_LD <- file.path(dir_output, "LD", "UKBB", paste0(gene_of_interest, "_LD.RDS"))
  ensure_dir(dirname(filename_LD))
  if (file.exists(filename_LD)) {
    LD <- readRDS(filename_LD)
  } else {
    filename_ld <- get_ld_path_UKBB(gene_of_interest, WINDOW_SIZE)
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
    ld_meta_filt <- ld_meta_filt[!duplicated(ld_meta_filt$snp), ]
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

load_ld_mat_FinnGen <- function(gene_of_interest, WINDOW_SIZE, IDs_to_keep, dir_output) {
  dir_LD <- "data/FinnGen/ld_matrix/"
  filename_LD <- file.path(dir_output, "LD", "FinnGen", paste0(gene_of_interest, "_LD.RDS"))
  ensure_dir(dirname(filename_LD))
  if (file.exists(filename_LD)) {
    LD <- readRDS(filename_LD)
  } else {
    region <- get_gene_region_GRCh38_UKBB(gene_of_interest, WINDOW_SIZE)

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
    if (length(IDs_final) == 0) {
      IDs_to_keep <- sapply(strsplit(IDs_to_keep, "_"), function(x) paste0(x[1], "_", x[2], "_", x[4], "_", x[3]))
      IDs_final <- Reduce(intersect, list(df_LD$ID1, df_LD$ID2, IDs_to_keep))
    }
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
    # if (nrow(df_merged) == 0) {
    #   df_UKBBrsIDmap$REF <- sapply(df_UKBBrsIDmap$REF, complement_allele)
    #   df_UKBBrsIDmap$ALT <- sapply(df_UKBBrsIDmap$ALT, complement_allele)
    #   df_UKBBrsIDmap <- df_UKBBrsIDmap %>% mutate(ID = paste0("chr", CHR, "_", POS, "_", REF, "_", ALT))
    #   df_merged <- df_LD_filt %>%
    #     inner_join(df_UKBBrsIDmap, by = c("ID1" = "ID")) %>%
    #     dplyr::rename(snp = rsid) %>%
    #     inner_join(df_UKBBrsIDmap[, c("ID", "rsid")], by = c("ID2" = "ID")) %>%
    #     dplyr::rename(snp2 = rsid) %>%
    #     dplyr::select(-ID1, -ID2) %>%
    #     distinct()
    # }

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

harmonize <- function(runID, gene_of_interest, WINDOW_SIZE, list_data, LD_type, dir_output) {
  names_risk_factor <- names(list_data)
  fname_harmonize <- file.path(dir_output, "harmonize", "harmonize.RDS")

  if (file.exists(fname_harmonize)) {
    res <- readRDS(fname_harmonize)
  } else {
    ensure_dir(dirname(fname_harmonize))

    rsids_to_keep <- Reduce("intersect", lapply(list_data, function(x) x$snp))

    n_samples <- sapply(list_data, function(x) mean(x$nsample, na.rm = T))
    max_sample <- max(n_samples)
    selected_factor <- sort(names(which(n_samples == max_sample)))[1]

    res <- list()
    for (this_factor in setdiff(names(list_data), selected_factor)) {
      merged <- list_data[[selected_factor]] %>%
        dplyr::inner_join(list_data[[this_factor]], by = "snp", suffix = c("_1", "_2")) %>%
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
    names(res) <- c(selected_factor, setdiff(names(list_data), selected_factor))

    if (length(LD_type) == 1) {
      if (LD_type == "FinnGen") {
        IDs_to_keep <- res[[selected_factor]]$ID
        LD <- load_ld_mat_FinnGen(gene_of_interest, WINDOW_SIZE, IDs_to_keep, dir_output)
      } else if (LD_type == "UKBB") {
        LD <- load_ld_mat_UKBB(gene_of_interest, WINDOW_SIZE, rsids_to_keep, dir_output)
      }
      ld_merged <- as_tibble(res[[selected_factor]]) %>% left_join(LD$meta, by = "snp", suffix = c("_1", "_2"))
      ld_merged <- ld_merged %>% dplyr::filter((effect_1 == effect_2 & other_1 == other_2) | (effect_1 == other_2 & other_1 == effect_2))
      ld_mat_harmonised <- harmonize_ld(LD$mat[ld_merged$snp, ld_merged$snp], ld_merged)
      list_LD <- list()
      for (id in names(res)) list_LD[[id]] <- ld_mat_harmonised
    } else {
      IDs_to_keep <- res[[selected_factor]]$ID
      LD <- load_ld_mat_FinnGen(gene_of_interest, WINDOW_SIZE, IDs_to_keep, dir_output)
      ld_merged <- as_tibble(res[[selected_factor]]) %>% inner_join(LD$meta, by = "snp", suffix = c("_1", "_2"))
      ld_merged <- ld_merged %>% dplyr::filter((effect_1 == effect_2 & other_1 == other_2) | (effect_1 == other_2 & other_1 == effect_2))

      LD2 <- load_ld_mat_UKBB(gene_of_interest, WINDOW_SIZE, rsids_to_keep, dir_output)
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
      res[[id]]$LD <- list_LD[[ifelse(id == "EGA", "UKBB", id)]]
      if (all(is.na(res[[id]]$s))) {
        res[[id]]$s <- NULL
      } else {
        res[[id]]$s <- unique(res[[id]]$s)
      }
    }

    res$names$risk_factors <- names_risk_factor
    saveRDS(res, file = fname_harmonize)
  }
  return(res)
}


run_propcoloc <- function(res, dir_results, this_J = 10) {
  RF1 <- res$names$risk_factors[1]
  RF2 <- res$names$risk_factors[2]

  fname_prop.coloc.res1 <- file.path(dir_results, "propcoloc", "prop.coloc.RDS")
  fname_prop.coloc.res2 <- file.path(dir_results, "propcoloc", "prop.coloc.txt")
  ensure_dir(dirname(fname_prop.coloc.res1))
  # if (file.exists(fname_prop.coloc.res1) & file.exists(fname_prop.coloc.res2)) {
  #   return()
  # }

  prop.coloc.res <- NULL
  prune_values <- c(0.4, 0.3, 0.2, 0.1, 0.01, 0.001, 0.0001, 0.00001)
  current_prune_index <- 1
  while (is.null(prop.coloc.res) && current_prune_index <= length(prune_values)) {
    this_prune <- prune_values[current_prune_index]
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
    if (is.null(prop.coloc.res)) {
      current_prune_index <- current_prune_index + 1
    } else {
      break
    }
  }
  if (is.null(prop.coloc.res)) prop.coloc.res <- list(p_full = NA, p_cond = NA, LM_full = NA, LM_cond = NA)

  prop.coloc.res <- tryCatch(
    {
      prop.coloc::prop.coloc(
        b1 = res[[RF1]]$beta, se1 = res[[RF1]]$se, b2 = res[[RF2]]$beta, se2 = res[[RF2]]$se,
        n = c(res[[RF1]]$N, res[[RF2]]$N),
        ld = res[[RF1]]$LD, figs = F, traits = c(RF1, RF2),
        prune = this_prune, J = this_J
      )
    },
    error = function(cond) {
      writeLines(capture.output(cond), paste0(fname_prop.coloc.res2, "_err"))
      return(cond)
    }
  )
  if (any(class(prop.coloc.res) == "error")) {
    return()
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
  ensure_dir(dirname(fname_output1))
  # if (file.exists(fname_output1) & file.exists(fname_output2)) {
  #   return()
  # }
  res[[RF1]]$MAF <- with(res[[RF1]], ifelse(effect_AF < .5, effect_AF, 1 - effect_AF))
  res[[RF2]]$MAF <- with(res[[RF2]], ifelse(effect_AF < .5, effect_AF, 1 - effect_AF))
  res[[RF1]]$MAF <- with(res[[RF1]], ifelse(MAF == 0, 1 / res[[RF1]]$N, MAF))
  res[[RF2]]$MAF <- with(res[[RF2]], ifelse(MAF == 0, 1 / res[[RF2]]$N, MAF))

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
  ensure_dir(dirname(fname_output1))
  # if (file.exists(fname_output1) & file.exists(fname_output2)) {
  #   return()
  # }
  res[[RF1]]$MAF <- with(res[[RF1]], ifelse(effect_AF < .5, effect_AF, 1 - effect_AF))
  res[[RF2]]$MAF <- with(res[[RF2]], ifelse(effect_AF < .5, effect_AF, 1 - effect_AF))
  res[[RF1]]$MAF <- with(res[[RF1]], ifelse(MAF == 0, 1 / res[[RF1]]$N, MAF))
  res[[RF2]]$MAF <- with(res[[RF2]], ifelse(MAF == 0, 1 / res[[RF2]]$N, MAF))

  selected <- c("beta", "varbeta", "position", "MAF", "snp", "type", "N", "LD")
  dat1 <- res[[RF1]][selected]
  dat2 <- res[[RF2]][selected]
  names(dat1$beta) <- names(dat1$varbeta) <- names(dat1$MAF) <- names(dat1$position) <- names(dat2$beta) <- names(dat2$varbeta) <- names(dat2$MAF) <- names(dat2$position) <- dat1$snp
  dat1$Z <- with(dat1, beta / sqrt(varbeta))
  dat2$Z <- with(dat2, beta / sqrt(varbeta))

  susie.res <- tryCatch(
    {
      out <- coloc::coloc.susie(dat1, dat2)
      if (is.null(out$summary) || is.null(out$summary$PP.H4.abf)) {
        stop("No result")
      }
      out
    },
    error = function(cond) {
      writeLines(capture.output(cond), paste0(fname_output1, "_err"))
      return(cond)
    }
  )
  if (any(class(susie.res) == "error")) {
    return()
  }
  susie.res$summary$ori_nsnps <- length(dat1$snp)

  write.table(susie.res$summary,
    fname_output1,
    quote = F, row.names = F, col.names = T, sep = "\t"
  )
  # pdf(fname_output2, width = 8, height = 5)
  # sensitivity(out, "H3 > 0.5 & H4 < 0.5", row = 1, dataset1 = dat1, dataset2 = dat2)
  # dev.off()
  return(susie.res$summary)
}


run_colocPropTest <- function(res, dir_results) {
  RF1 <- res$names$risk_factors[1]
  RF2 <- res$names$risk_factors[2]

  fname_output <- file.path(dir_results, "colocPropTest", "colocPropTest.txt")
  ensure_dir(dirname(fname_output))
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
  res_propcoloc <- "insufficient"
  if (file.exists(fname_propcoloc)) {
    df_temp <- readRDS(fname_propcoloc)
    p_het <- ifelse(is.factor(df_temp$p_cond) & df_temp$p_cond == F, 1, df_temp$p_cond)
    p_slope <- ifelse(is.factor(df_temp$LM_cond) & df_temp$LM_cond == F, 1, df_temp$LM_cond)
    res_propcoloc <- ifelse(is.na(p_slope) | p_slope > 0.05, "insufficient",
      ifelse(p_het > 0.05, "coloc", "non_coloc")
    )
  }

  return(res_propcoloc)
}



get_susie_res <- function(dir_results, list_factors) {
  RF1 <- list_factors[1]
  RF2 <- list_factors[2]

  fname_susie <- file.path(dir_results, "susie", "susie.txt")
  res_susie <- "insufficient"
  if (file.exists(fname_susie)) {
    res_temp <- with(
      read.delim(fname_susie),
      data.frame(RF1, RF2,
        H0 = PP.H0.abf,
        H1 = PP.H1.abf,
        H2 = PP.H2.abf,
        H3 = PP.H3.abf,
        H4 = PP.H4.abf
      )
    )

    res_susie <- ifelse(any(res_temp$H4 >= 0.5), "coloc",
      ifelse(all(res_temp$H4 < 0.5) && any(res_temp$H3 >= 0.5), "non_coloc", "insufficient")
    )
  }
  return(res_susie)
}

get_coloc_res <- function(dir_results, list_factors) {
  RF1 <- list_factors[1]
  RF2 <- list_factors[2]
  fname_coloc <- file.path(dir_results, "coloc", "coloc.txt")
  res_coloc <- "insufficient"
  if (file.exists(fname_coloc)) {
    res_coloc <- with(
      as.data.frame(t(read.delim(fname_coloc))),
      ifelse(PP.H4.abf >= 0.5, "coloc", ifelse(PP.H3.abf >= 0.5, "non_coloc", "insufficient"))
    )
  }
  return(res_coloc)
}

get_colocPropTest_res <- function(dir_results, list_factors) {
  res_colocPropTest <- "insufficient"
  RF1 <- list_factors[1]
  RF2 <- list_factors[2]

  fname_colocPropTest <- file.path(dir_results, "colocPropTest", "colocPropTest.txt")
  if (file.exists(fname_colocPropTest) && file.info(fname_colocPropTest)$size > 1) {
    df1 <- read.delim(fname_colocPropTest)
    res_temp <- data.frame(RF1 = RF1, RF2 = RF2, min_p = min(df1$p), min_fdr = min(df1$fdr))
    res_colocPropTest <- ifelse(res_temp$min_fdr > 0.05, "coloc", "non_coloc")
  }
  return(res_colocPropTest)
}


get_all_results <- function(runID, names_risk_factor) {
  lst <- list.files(file.path("results", runID))
  merged <- NULL
  for (gene_of_interest in lst) {
    dir_results <- file.path("results", runID, gene_of_interest)
    res_propcoloc <- get_propcoloc_res(dir_results, names_risk_factor)
    res_susie <- get_susie_res(dir_results, names_risk_factor)
    res_coloc <- get_coloc_res(dir_results, names_risk_factor)
    res_colocPropTest <- get_colocPropTest_res(dir_results, names_risk_factor)

    merged <- rbind(merged, data.frame(
      gene = gene_of_interest,
      propcoloc = res_propcoloc,
      coloc = res_coloc,
      susie = res_susie,
      colocPropTest = res_colocPropTest
    ))
  }

  summary_coloc <- t(data.frame(prop.table(table(factor(merged$coloc, levels = c("coloc", "non_coloc", "insufficient"))))))[2, , drop = F]
  summary_susie <- t(data.frame(prop.table(table(factor(merged$susie, levels = c("coloc", "non_coloc", "insufficient"))))))[2, , drop = F]
  summary_propcoloc <- t(data.frame(prop.table(table(factor(merged$propcoloc, levels = c("coloc", "non_coloc", "insufficient"))))))[2, , drop = F]
  summary_colocPropTest <- t(data.frame(prop.table(table(factor(merged$colocPropTest, levels = c("coloc", "non_coloc", "insufficient"))))))[2, , drop = F]
  summary_combined <- Reduce("cbind", list(summary_coloc, summary_susie, summary_propcoloc, summary_colocPropTest))
  summary_combined <- data.frame(runID, summary_combined)

  methods <- c("coloc", "susie", "propcoloc", "colocPropTest")
  categories <- c("coloc", "non_coloc", "insufficient")
  cat_abbr <- c("C", "NC", "IS")

  n_methods <- length(methods)
  n_categories <- length(categories)
  row_names <- paste(rep(methods, each = n_categories), cat_abbr, sep = ".")
  col_names <- paste(rep(methods, each = n_categories), rep(cat_abbr, n_methods), sep = ".")
  result_table <- matrix(NA,
    nrow = n_methods * n_categories, ncol = n_methods * n_categories,
    dimnames = list(row_names, col_names)
  )

  for (i in 1:n_methods) {
    for (j in 1:n_methods) {
      method1 <- methods[i]
      method2 <- methods[j]

      for (k in 1:n_categories) {
        for (l in 1:n_categories) {
          cat1 <- categories[k]
          cat2 <- categories[l]
          match_rate <- mean(merged[[method1]] == cat1 & merged[[method2]] == cat2, na.rm = TRUE)
          row_idx <- (i - 1) * n_categories + k
          col_idx <- (j - 1) * n_categories + l
          result_table[row_idx, col_idx] <- round(match_rate, 2) # 소수점 3자리로 반올림
        }
      }
    }
  }
  write.table(result_table, paste0("results/", runID, ".txt"), quote = F, row.names = F, col.names = F, sep = "\t")
  return(summary_combined)
}

run_colocalization_analysis <- function(runID, gene_of_interest, WINDOW_SIZE, list_data, LD_type) {
  dir_output <- file.path("results", runID, gene_of_interest)
  res <- harmonize(runID, gene_of_interest, WINDOW_SIZE, list_data, LD_type, dir_output)

  if (length(LD_type) == 1) {
    run_coloc(res, dir_output)
    run_susie(res, dir_output)
    run_propcoloc(res, dir_output)
    run_colocPropTest(res, dir_output)
  } else if (length(LD_type) == 2) {
    LD_UKBB <- res[names(res) %in% c("EGA", "UKBB")][[1]]$LD
    LD_FinnGen <- res[["FinnGen"]]$LD

    run_coloc(res, dir_output)
    run_susie(res, dir_output)

    runID0 <- paste0(runID, "_UKBBLD")
    dir_output0 <- file.path("results", runID0, gene_of_interest)
    res[[1]]$LD <- res[[2]]$LD <- LD_UKBB
    run_propcoloc(res, dir_output0)
    run_colocPropTest(res, dir_output0)
    run_susie(res, dir_output0)

    runID1 <- paste0(runID, "_FinnGenLD")
    dir_output1 <- file.path("results", runID1, gene_of_interest)
    res[[1]]$LD <- res[[1]]$LD <- LD_FinnGen
    run_propcoloc(res, dir_output1)
    run_colocPropTest(res, dir_output1)
    run_susie(res, dir_output1)
  }
}
