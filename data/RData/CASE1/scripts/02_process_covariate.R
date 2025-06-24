library(data.table)

# Load UK Biobank phenotype data
bd <- fread("/path/to/ukb672224.tab", header=TRUE, sep="\t")

# Create covariate dataframe including sex, age, genotyping array, assessment center,
# and the first 20 genetic principal components (PCs)
data_to_write = with(bd, data.frame(FID = f.eid,
				    IID = f.eid,
				    sex = f.31.0.0,
				    age = f.21022.0.0,
				    age_sq = f.21022.0.0**2,
				    assessment_center = f.54.0.0,
				    genotyping_array = f.22000.0.0,
				    pc1 = f.22009.0.1,
				    pc2 = f.22009.0.2,
				    pc3 = f.22009.0.3,
				    pc4 = f.22009.0.4,
				    pc5 = f.22009.0.5,
				    pc6 = f.22009.0.6,
				    pc7 = f.22009.0.7,
				    pc8 = f.22009.0.8,
				    pc9 = f.22009.0.9,
				    pc10 = f.22009.0.10,
				    pc11 = f.22009.0.11,
				    pc12 = f.22009.0.12,
				    pc13 = f.22009.0.13,
				    pc14 = f.22009.0.14,
				    pc15 = f.22009.0.15,
				    pc16 = f.22009.0.16,
				    pc17 = f.22009.0.17,
				    pc18 = f.22009.0.18,
				    pc19 = f.22009.0.19,
				    pc20 = f.22009.0.20))

write.table(data_to_write, "your/data/directory/covariate.txt", quote= F, row.names=F, col.names = T, sep = "\t")

