# Set the path for the sample files
sample_bgen="/path/to/ukb22828_c1_b0_v3_s487160.sample"
file_withdraw1="/path/to/withdraw98032_19.txt"
file_withdraw2="/path/to/withdraw98032_20241217.csv"
path_output="/your/data/directory/sample_group"

# Load BGEN sample file (skip the first header row)
df_bgen = read.table(sample_bgen, header = T)
df_bgen = df_bgen[-1,]

# Load withdrawn participant IDs
df_withdraw1 = read.table(file_withdraw1)
df_withdraw2 = read.table(file_withdraw2)

# Number of samples
#> dim(df_plink)
#[1] 488377      6
#> dim(df_bgen)
#[1] 487409      4

# Exclude withdrawn samples
all_ids <- setdiff(df_bgen[,1], c(df_withdraw1[,1], df_withdraw2[,1]))

# Randomly split remaining samples into two equal-sized groups
set.seed(1)
group1_ids <- sample(all_ids, length(all_ids) / 2)
group2_ids <- setdiff(all_ids, group1_ids)

# Save sample ID files in PLINK-compatible format
# Output format example:
# FID     IID
# 4559402 4559402
# 1294479 1294479
write.table(data.frame(FID = all_ids, IID = all_ids), file.path(path_output, "all_samples.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(data.frame(FID = group1_ids, IID = group1_ids), file.path(path_output, "group1.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(data.frame(FID = group2_ids, IID = group2_ids), file.path(path_output, "group2.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
