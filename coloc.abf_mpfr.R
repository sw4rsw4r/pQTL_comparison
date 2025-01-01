coloc.abf_mpfr <- function (dataset1, dataset2, MAF = NULL, p1 = 1e-04, p2 = 1e-04, p12 = 1e-05) {
    if (!("MAF" %in% names(dataset1)) & !is.null(MAF)) 
        dataset1$MAF <- MAF
    if (!("MAF" %in% names(dataset2)) & !is.null(MAF)) 
        dataset2$MAF <- MAF
    check_dataset(d = dataset1, 1)
    check_dataset(d = dataset2, 2)
    df1 <- process.dataset_mpfr(d = dataset1, suffix = "df1")
    df2 <- process.dataset_mpfr(d = dataset2, suffix = "df2")
    p1 = adjust_prior(p1, nrow(df1), "1")
    p2 = adjust_prior(p2, nrow(df2), "2")
    merged.df <- merge(df1, df2)
    p12 = adjust_prior(p12, nrow(merged.df), "12")
    if (!nrow(merged.df)) 
        stop("dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified")
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
        df <- cbind(df, abf)
        if ("position" %in% nd) 
            df <- cbind(df, position = d$position)
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
    z <- qnorm(0.5 * p, lower.tail = FALSE)
    ## Shrinkage factor: ratio of the prior variance to the total variance
    r <- sd.prior^2 / (sd.prior^2 + V)
    ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
    lABF = 0.5 * (log(1-r) + (r * z^2))
    ret <- data.frame(V,z,r,lABF)
    if(!is.null(suffix))
        colnames(ret) <- paste(colnames(ret), suffix, sep=".")
    return(ret)
}

Var.data <- function(f, N) {
    1 / (2 * N * f * (1 - f))
}