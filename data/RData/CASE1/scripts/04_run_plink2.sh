#!/bin/bash
#SBATCH --array=0-695%10
#SBATCH -J run_plink
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem=50G
#SBATCH --time 12:00:00
#SBATCH --mail-type TIME_LIMIT
#SBATCH -p cclake-himem
#SBATCH --output=logs/run_plink_%A_%a.out
#SBATCH --error=logs/run_plink_%A_%a.err

. /etc/profile.d/modules.sh # Leave this line (enables the module command)
module purge                # Removes all modules still loaded

module load plink/2.00-alpha
module load bgen/1.2

path_plink='/path/to/genotypes'
output='/path/to/results'
mkdir -p $output logs

# Get the list of protein files
FILE_GENE="/path/to/genes_to_run154.txt"
LIST_PHENO=($(cat $FILE_GENE | awk '{print $1}'))
# Select the protein file for the current task
GENE=${LIST_PHENO[$((SLURM_ARRAY_TASK_ID - 1))]}

# Extract CHR, START, END from the gene coordinates file
FILE_GENE_ANNO="/path/to/olink_protein_map_3k_v1_liftOverToHg19.tsv"
read CHR START END < <(grep "^${GENE}\b" ${FILE_GENE_ANNO} | awk 'BEGIN {FS="\t"} {print $2, $3, $4}')

# Set window size (50kb)
WINDOW_SIZE=50000

# Calculate locus lower and upper region with window size
locus_lower=$((START - WINDOW_SIZE))
locus_upper=$((END + WINDOW_SIZE))

# Ensure locus_lower is not lower than 0
if [ "$locus_lower" -lt 1 ]; then
    locus_lower=0
fi

# Skip if coordinates are missing
if [ -z "$CHR" ] || [ -z "$locus_lower" ] || [ -z "$locus_upper" ]; then
    echo "Gene ${GENE} not found in coordinate file. Skipping."
    exit 1
fi

path_pheno="/your/data/directory/protein_level/"
file_pheno=${path_pheno}/${GENE}.txt
path_groups="/your/data/directory/sample_group/"
file_cov="/your/data/directory/covariate.txt"
cols_cov="--covar-name sex,age,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10,pc11,pc12,pc13,pc14,pc15,pc16,pc17,pc18,pc19,pc20"

#bgenix -g ${path_bgen}/"ukb22828_c${CHR}_b0_v3.bgen" \
#        -incl-range ${CHR}":"${locus_lower}"-"${locus_upper} > ${output}/${GENE}.bgen

# Run PLINK separately for each group (e.g., group1 and group2)
groups=("group1" "group2")
for group in "${groups[@]}"; do

    plink2 \
       --bed ${path_plink}/ukb22418_c${CHR}_b0_v2.bed \
       --bim ${path_plink}/ukb_snp_chr${CHR}_v2.bim \
       --fam ${path_plink}/ukb22418_c${CHR}_b0_v2_s488131.fam \
       --chr ${CHR} \
       --from-bp ${locus_lower} \
       --to-bp ${locus_upper} \
       --pheno ${file_pheno} \
       --pheno-name "protein_level" \
       --covar ${file_cov} \
       ${cols_cov} \
       --maf 0.001 \
       --geno 0.1 \
       --hwe 1e-06 \
       --glm hide-covar \
       --freq \
       --keep ${path_groups}/${group}.txt \
       --threads 4 \
       --out ${output}/${GENE}_${group}
       #--bgen ${output}/${GENE}.bgen ref-first \
       #--sample ${path_bgen}/"ukb22828_c${CHR}_b0_v3_s487160.sample" \

done

echo "Finished processing: $GENE"
