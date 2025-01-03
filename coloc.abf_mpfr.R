coloc.abf_mpfr <- function (dataset1, dataset2, MAF = NULL, p1 = 1e-04, p2 = 1e-04, p12 = 1e-05) {
    if (!("MAF" %in% names(dataset1)) & !is.null(MAF)) 
        dataset1$MAF <- MAF
    if (!("MAF" %in% names(dataset2)) & !is.null(MAF)) 
        dataset2$MAF <- MAF
    check_dataset(d = dataset1, 1)
    check_dataset(d = dataset2, 2)
    df1 <- process.dataset_mpfr(d = dataset1, suffix = "df1")
    df2 <- process.dataset_mpfr(d = dataset2, suffix = "df2")
    p1 = adjust_prior(p1, length(df1$snp), "1")
    p2 = adjust_prior(p2, length(df2$snp), "2")

    if (!all(df1$snp == df2$snp))
        stop("dataset1 and dataset2 should contain the same snps in the same order")
    df2$snp = df2$position = NULL
    merged.df <- c(df1, df2)
    p12 = adjust_prior(p12, length(merged.df$snp), "12")
    merged.df$internal.sum.lABF <- with(merged.df, lABF.df1 + 
        lABF.df2)
    my.denom.log.abf <- logsum(merged.df$internal.sum.lABF)
    merged.df$SNP.PP.H4 <- exp(merged.df$internal.sum.lABF - 
        my.denom.log.abf)
    pp.abf <- combine.abf(merged.df$lABF.df1, merged.df$lABF.df2, 
        p1, p2, p12)
    common.snps <- nrow(merged.df)
    results <- c(nsnps = common.snps, pp.abf)
    output <- list(summary = results, results = merged.df, priors = c(p1 = p1, 
        p2 = p2, p12 = p12))
    class(output) <- c("coloc_abf", class(output))
    return(output)
}

process.dataset_mpfr <- function (d, suffix) 
{
    nd <- names(d)
    if ("beta" %in% nd && "varbeta" %in% nd) {
        if (d$type == "quant" && !("sdY" %in% nd)) 
            d$sdY <- sdY.est(d$varbeta, d$MAF, d$N)
        df <- approx.bf.estimates(z = d$beta/sqrt(d$varbeta), 
            V = d$varbeta, type = d$type, suffix = suffix, sdY = d$sdY)
        df$snp <- as.character(d$snp)
        if ("position" %in% nd) 
            df <- cbind(df, position = d$position)
        return(df)
    }
    if ("pvalues" %in% nd & "MAF" %in% nd & "N" %in% nd) {
        df <- list(pvalues = d$pvalues, MAF = d$MAF, N = d$N,
            snp = as.character(d$snp))
        snp.index <- which(names(df) == "snp")
        names(df)[-snp.index] <- paste(names(df)[-snp.index],
            suffix, sep = ".")
        abf <- approx.bf.p_mpfr(p = df$pvalues, f = df$MAF, type = d$type,
            N = df$N, s = d$s, suffix = suffix)
        df <- c(df, abf)
        if ("position" %in% nd) 
            df <- c(df, list(position = d$position))
        return(df)
    }
    stop("Must give, as a minimum, one of:\n(beta, varbeta, type, sdY)\n(beta, varbeta, type, MAF)\n(pvalues, MAF, N, type)")
}


# https://rdrr.io/cran/coloc/src/R/claudia.R#sym-approx.bf.p
approx.bf.p_mpfr <- function(p,f,type, N, s, suffix=NULL) {
    if(type=="quant") {
        sd.prior <- 0.15
        V <- Var.data(f, N)
    } else {
        sd.prior <- 0.2
        V <- Var.data.cc(f, N, s)
    }
    # FIX Here!!! ===
    z <- qnorm(0.5 * p, lower.tail = FALSE)
    # FIX Here!!! ===

    ## Shrinkage factor: ratio of the prior variance to the total variance
    r <- sd.prior^2 / (sd.prior^2 + V)
    ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
    lABF = 0.5 * (log(1-r) + (r * z^2))
    ret <- list(V = V,z = z,r = r,lABF = lABF)
    if(!is.null(suffix))
        names(ret) <- paste(names(ret), suffix, sep=".")
    return(ret)
}

Var.data <- function(f, N) {
    1 / (2 * N * f * (1 - f))
}


adjust_prior=function(p,nsnps,suffix="") {
    if(nsnps * p >= 1) { ## for very large regions
        warning(paste0("p",suffix," * nsnps >= 1, setting p",suffix,"=1/(nsnps + 1)"))
        1/(nsnps + 1)
    } else {
        p
    }
}

logsum <- function(x) {
    my.max <- max(x)                              ##take out the maximum value in log form
    my.res <- my.max + log(sum(exp(x - my.max ))) 
    return(my.res)
}

combine.abf <- function(l1, l2, p1, p2, p12, quiet=FALSE) {
    stopifnot(length(l1)==length(l2))
    lsum <- l1 + l2
    lH0.abf <- 0
    lH1.abf <- log(p1) + logsum(l1)
    lH2.abf <- log(p2) + logsum(l2)
    lH3.abf <- log(p1) + log(p2) + logdiff(logsum(l1) + logsum(l2), logsum(lsum))
    lH4.abf <- log(p12) + logsum(lsum)

    all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
    my.denom.log.abf <- logsum(all.abf)
    pp.abf <- exp(all.abf - my.denom.log.abf)
    names(pp.abf) <- paste("PP.H", (1:length(pp.abf)) - 1, ".abf", sep = "")
    if(!quiet) {
        print(signif(pp.abf,3))
        print(paste("PP abf for shared variant: ", signif(pp.abf["PP.H4.abf"],3)*100 , '%', sep=''))
    }
    return(pp.abf)
}

logdiff <- function(x,y) {
    my.max <- max(x,y)                              ##take out the maximum value in log form
    my.res <- my.max + log(exp(x - my.max ) - exp(y-my.max))
    return(my.res)
}
