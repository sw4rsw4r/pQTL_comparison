obj <- readRDS("data/RData/UKBB/A1BG_window_5e+04_processed.RDS")
snp_list <- with(obj[["Olink"]], sub("^chr", "", gsub("_", "-", ID)))

batch_size <- 10
snp_batches <- split(snp_list, ceiling(seq_along(snp_list) / batch_size))

result_list <- list()
for (batch in snp_batches) {
    batch_result <- getVariantPopData(varids = batch, genomes = "GRCh38")
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
merged <- merge(obj[[1]], final_df, by.x = "position", by.y = "pos")

with(merged, plot(ifelse(effect_AF < .5, effect_AF, 1 - effect_AF), ifelse(AF < .5, AF, 1 - AF), main = "A1BG", xlab = "MAF", ylab = "gnomAD (NFE, Non-Finnish European)"))
