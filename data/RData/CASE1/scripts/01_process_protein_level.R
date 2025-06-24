library(data.table)

# Define output directory for processed protein level data and create it if needed
path_data = "/your/data/directory/protein_level/"
dir.create(path_data, recursive = T)


# Load raw proteomics data files (replace with actual paths or mount points as needed)
olink = fread("/path/to/olink_data.txt.gz")
olink_samples = fread("/path/to/ukb676343.csv")
olink_coding = fread("/path/to/coding143.tsv")

# Extract gene names and descriptions from the 'meaning' column
olink_coding[,gene:= gsub(";.*","",meaning)]
olink_coding[,description:= gsub(".*;","",meaning)]

for (genename in olink_coding$gene){
  # Subset data for the current gene's protein at baseline (ins_index == 0)
  data_subset = olink[protein_id == olink_coding[gene==genename,coding] & ins_index==0,]
  # Prepare data frame with sample IDs and protein levels
  data_to_write = with(data_subset, data.frame(FID = eid, IID = eid, protein_level = result))
  # Write gene-specific protein level data to a tab-delimited text file
  filename = paste0(path_data, "/", genename, ".txt")
  write.table(data_to_write, filename, quote= F, row.names=F, col.names = T, sep = "\t")

  # Example output format (PLINK phenotype format):
  # FID     IID     protein_level
  # 1000250 1000250 1.1816
  # 1000407 1000407 2.5131
  # 1000499 1000499 0.8153
  # 1000679 1000679 2.2766
}
