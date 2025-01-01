coloc.abf_mpfr <- function (dataset1, dataset2, MAF = NULL, p1 = 1e-04, p2 = 1e-04, p12 = 1e-05) {
    if (!("MAF" %in% names(dataset1)) & !is.null(MAF)) 
        dataset1$MAF <- MAF
    if (!("MAF" %in% names(dataset2)) & !is.null(MAF)) 
        dataset2$MAF <- MAF
    check_dataset(d = dataset1, 1)
    check_dataset(d = dataset2, 2)
    df1 <- process.dataset(d = dataset1, suffix = "df1")
    df2 <- process.dataset(d = dataset2, suffix = "df2")
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