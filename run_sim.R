source("load-hapdata.R")
library(randomFunctions)
library(magrittr)
library(coloc)

####################################################################
source("utils.R")
suppressMessages(library(prop.coloc))
suppressMessages(library(colocPropTest))
suppressMessages(library(easyGgplot2))
suppressMessages(library(gnomadR))
suppressMessages(library(reticulate))
suppressMessages(library(synapser))
suppressMessages(library(synapserutils))


#####################################################################

## simulations options
## 5 = share 1, weakest effect for one, strongest for other
## 4 = share 2, opposite effects
## 3 = share 1, equal effect (weakest) + indep each
## 2 = share 2, equal effects
## 1 = share 1, equal effect (strongest), + indep each
## 0 = share 0

for (ld in c("highld", "lowld")) {
  for (NCV in 3:4) {
    # if (NCV == 3) list_SPECIAL <- 0:3
    # if (NCV == 4) list_SPECIAL <- 0:5
    list_SPECIAL <- 0
    for (SPECIAL in list_SPECIAL) {
      for (NSIM in 1:1000) {
        set.seed((NSIM))
        N <- 1000
        NSNP <- 1000

        test <- setRefClass("test", fields = list(N = "numeric", NSIM = "numeric", NCV = "numeric", NSNP = "numeric", SPECIAL = "numeric", ld = "character", pop = "character"))
        args <- test(N = N, NSIM = NSIM, NCV = NCV, NSNP = NSNP, SPECIAL = SPECIAL, ld = ld, pop = "eur")

        print(args)

        hdata <- readRDS(paste0("~/scratch/", args$ld, ".RDS"))
        LD <- hdata$LD
        LD.alt <- hdata$LD.alt
        dfsnps <- hdata$snps
        h <- hdata$h
        makeD <- function(y, X, best = NULL) {
          if (is.null(best)) {
            m <- snp.rhs.estimates(y ~ 1, snp.data = X, family = "Gaussian")
          } else {
            df <- data.frame(y = y, best = as(X[, best], "numeric"), row.names = rownames(X))
            m <- snp.rhs.estimates(y ~ ., data = df, snp.data = X, family = "Gaussian")
          }
          nulls <- sapply(m, is.null)
          if (any(nulls)) {
            m <- m[!nulls]
          }
          b <- sapply(m, "[[", "beta")
          v <- sapply(m, "[[", "Var.beta")
          z <- b / sqrt(v)
          list(
            N = args$N,
            MAF = dfsnps$maf[!nulls],
            beta = b,
            varbeta = v,
            type = "quant",
            sdY = sd(y),
            snp = dfsnps$id[!nulls]
          )
        }

        simone <- function() {
          pr.var1 <- rnorm(1, 0, 0.15) %>% abs()
          pr.var2 <- rnorm(1, 0, 0.15) %>% abs()
          if (args$NCV == 4) { ## sample common CVs
            if (args$SPECIAL %in% c("2", "4")) {
              CV <- sample(which(dfsnps$maf > 0.1 & dfsnps$maf < 0.9), 2)
              nCV1 <- nCV2 <- dfsnps$id[CV]
              beta1 <- beta2 <- sort(sample(c(1:9) / 6, 2))
              if (args$SPECIAL == "4") {
                beta2 <- rev(beta2)
              }
            } else if (args$SPECIAL == "0") {
              CV <- sample(which(dfsnps$maf > 0.1 & dfsnps$maf < 0.9), 4)
              nCV1 <- dfsnps$id[CV[1:2]]
              nCV2 <- dfsnps$id[CV[3:4]]
              beta1 <- sort(sample(c(1:9) / 6, 2))
              beta2 <- beta1 # sort(sample(c(1:6)/6,args$NCV))
            } else if (args$SPECIAL %in% c("1", "3", "5")) {
              CV <- sample(which(dfsnps$maf > 0.1 & dfsnps$maf < 0.9), 3)
              nCV1 <- dfsnps$id[CV[1:2]] # A B
              nCV2 <- dfsnps$id[CV[c(1, 3)]] # A C
              b <- sort(sample(c(1:9) / 6, 2))
              if (args$SPECIAL == "1") {
                beta1 <- b[c(2, 1)]
                beta2 <- b[c(2, 1)]
              } else if (args$SPECIAL == "3") {
                beta1 <- b[c(1, 2)]
                beta2 <- b[c(1, 2)]
              } else if (args$SPECIAL == "5") {
                beta1 <- b[c(2, 1)]
                beta2 <- b[c(1, 2)]
              }
            } else {
              stop("special not programmed yet: ", args$SPECIAL)
            }
          } else if (args$NCV == 3) {
            if (args$SPECIAL == "0") {
              CV <- sample(which(dfsnps$maf > 0.1 & dfsnps$maf < 0.9), 3)
              nCV1 <- dfsnps$id[CV[1:2]] # A, B
              nCV2 <- dfsnps$id[CV[3]] # C
              beta1 <- sort(sample(c(1:9) / 6, 2))
              beta2 <- beta1[2] # sort(sample(c(1:6)/6,1))
            } else {
              CV <- sample(which(dfsnps$maf > 0.1 & dfsnps$maf < 0.9), 2)
              nCV1 <- dfsnps$id[CV] # A, B
              nCV2 <- nCV1[1] # A
              b <- sort(sample(c(1:9) / 6, 2))
              if (args$SPECIAL == "1") {
                beta1 <- b[c(2, 1)]
                beta2 <- b[c(2)] # share strong
              } else if (args$SPECIAL == "2") {
                beta1 <- b[c(1, 2)]
                beta2 <- b[1] # share weak
              } else if (args$SPECIAL == "3") {
                beta1 <- b[c(2, 1)]
                beta2 <- b[1] # share strong but with weak effect
              }
            }
          } else {
            stop("nCV not programmed yet: ", args$nCV)
          }

          ## make genotypes
          nr <- nrow(h)
          G1 <- h[sample(1:nr, args$N, replace = TRUE), ] + h[sample(1:nr, args$N, replace = TRUE), ]
          G2 <- h[sample(1:nr, args$N, replace = TRUE), ] + h[sample(1:nr, args$N, replace = TRUE), ]
          colnames(G1) <- colnames(G2) <- dfsnps$id
          rownames(G1) <- rownames(G2) <- paste0("I", 1:args$N)
          X1 <- new("SnpMatrix", G1 + 1)
          X2 <- new("SnpMatrix", G2 + 1)

          ## outcomes
          y1 <- rnorm(args$N) + tcrossprod(G1[, nCV1], t(beta1))
          y2 <- rnorm(args$N) + tcrossprod(G2[, nCV2], t(beta2))

          # hist(y,breaks=100)
          usedata <- function(A) max(abs(A$beta) / sqrt(A$varbeta)) > 4.89
          ## stepwise regressions
          A1 <- makeD(y1, X1) # all signals, first dataset
          A2 <- makeD(y2, X2) # all signals, first dataset
          return(list(A1 = A1, A2 = A2, nCV1 = nCV1, nCV2 = nCV2))
        }
        data <- simone()

        runID <- paste0("SPECIAL_", SPECIAL, "_N_", N, "_NCV_", NCV, "_", ld)
        dir_output1 <- file.path("results", runID, paste0("sim", NSIM))
        # dir_output1 <- file.path("test_NC_propcoloc_upgraded", runID, paste0("sim", NSIM))

        res <- list(A1 = data$A1, A2 = data$A2, names = list(risk_factors = c("A1", "A2")))
        res[[1]]$effect_AF <- data$A1$MAF
        res[[2]]$effect_AF <- data$A2$MAF
        res[[1]]$LD <- res[[2]]$LD <- LD
        res[[1]]$type <- res[[2]]$type <- "quant"
        res[[1]]$position <- res[[2]]$position <- 1:length(res[[1]]$snp)
        res[[1]]$effect <- res[[2]]$position <- 1:length(res[[1]]$snp)
        res[[1]]$effect <- res[[2]]$effect <- NA
        res[[1]]$other <- res[[2]]$other <- NA
        res[[1]]$se <- sqrt(res[[1]]$varbeta)
        res[[2]]$se <- sqrt(res[[2]]$varbeta)
        res[[1]]$pval <- 2 * (1 - pnorm(abs(res[[1]]$beta / res[[1]]$se)))
        res[[2]]$pval <- 2 * (1 - pnorm(abs(res[[2]]$beta / res[[2]]$se)))


        run_coloc(res, dir_output1) # coloc does not use LD
        run_susie(res, dir_output1) # SuSiE uses LD, per dataset
        run_propcoloc(res, dir_output1) # propcoloc uses one shared LD
        run_colocPropTest(res, dir_output1) # colocPropTest uses one shared LD
        run_sharepro(res, dir_output1) # colocPropTest uses one shared LD
      }
    }
  }
}


# Constants

methods <- c("coloc", "susie", "sharepro", "propcoloc", "colocPropTest")
categories <- c("coloc", "non_coloc", "insufficient")
cat_abbr <- c("C", "NC", "IS")

get_merged_combined <- function(PATH) {
  merged_combined <- NULL
  list_runIDs <- list.files(PATH)
  list_runIDs <- grep("SPECIAL_0_", list_runIDs, value = T)
  for (runID in list_runIDs) {
    merged <- NULL
    list_sim <- paste0("sim", 1:1000)
    for (simID in list_sim) {
      dir_results <- file.path(PATH, runID, simID)

      res_propcoloc <- get_propcoloc_res(dir_results)$coloc
      res_susie <- get_susie_res(dir_results)$coloc
      res_coloc <- get_coloc_res(dir_results)$coloc
      res_colocPropTest <- get_colocPropTest_res(dir_results)$coloc
      res_sharepro <- get_sharepro(dir_results)$coloc

      merged <- rbind(merged, data.frame(
        runID = runID,
        simID = simID,
        propcoloc = res_propcoloc,
        coloc = res_coloc,
        susie = res_susie,
        colocPropTest = res_colocPropTest,
        sharepro = res_sharepro
      ))
    }

    summary_list <- lapply(methods, function(method) {
      out <- round(prop.table(table(factor(merged[[method]], levels = c("coloc", "non_coloc", "insufficient")))), 3)
      t(as.data.frame(out))[2, , drop = FALSE]
    })

    summary_combined <- do.call(cbind, summary_list)

    summary_combined <- data.frame(runID, summary_combined)

    colnames(summary_combined) <- c("runID", paste0(rep(methods, rep(3, length(methods))), "-", cat_abbr))
    rownames(summary_combined) <- summary_combined$runID

    merged_combined <- rbind(merged_combined, summary_combined)
  }
  return(merged_combined)
}

merged_combined_sim = get_merged_combined("results")
write.table(merged_combined_sim, "merged_combined_propcoloc_upgraded.txt", quote = F, row.names = F, col.names = T, sep = "\t")

