# Check if BiocManager is installed and install it if necessary

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

# Load the biomaRt package
library(biomaRt)
# Connect to the ENSEMBL Variation Mart for human SNP data using the GRCh37 genome version
snpMart = useEnsembl(biomart = "snps", dataset = "hsapiens_snp",version = 'GRCh37')
attributes = listAttributes(snpMart) # Retrieve available attributes in the Mart
attributes[1:5,]  # View the first 5 attributes

# Load GWAS summary statistics
sumstat <- fread("C:/Users/rasyed2/OneDrive/Documents/LDscore_V2/DrinksPerWeek_V2.txt.gz")
head(sumstat); dim(sumstat) 

# SNP_ID   CHR    POS     A1     A2        Beta          SE  Pval
# <char> <int>  <int> <char> <char>       <num>       <num> <num>
#   1: rs3094315     1 752566      A      G 0.001977795 0.002592506 0.446
# 2: rs3131972     1 752721      G      A 0.002114778 0.002574026 0.411
# 3: rs3131969     1 754182      G      A 0.001618614 0.002779992 0.560
# 4: rs1048488     1 760912      T      C 0.001855046 0.002592506 0.474
# 5: rs3115850     1 761147      C      T 0.001846789 0.002598792 0.477
# 6: rs2286139     1 761732      T      C 0.001268789 0.002686811 0.637
# [1] 1200257       8

# Query SNP metadata from biomaRt using SNP IDs from the summary statistics
snp_sumstat = getBM(
  attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'allele','allele_1','minor_allele','minor_allele_freq'),
  filters = 'snp_filter',
  values = sumstat$SNP_ID,
  mart = snpMart
)

# Preview the queried SNP data
head(snp_sumstat); dim(snp_sumstat)
snp_sumstat[snp_sumstat$refsnp_id == "rs3094315",]

# Process SNP data for compatibility with GWAS data
dat <- as.data.frame(snp_sumstat)
dat$allele_split <- strsplit(as.character(dat$allele), "/") # Split allele column into individual alleles

# Extract the first two alleles (A1 and A2) from the split data
dat$A1 <- sapply(dat$allele_split, function(x) ifelse(length(x) > 0, x[1], NA))
dat$A2 <- sapply(dat$allele_split, function(x) ifelse(length(x) > 1, x[2], NA))

# Reformat and rename columns for clarity
dat <- dat[,c("refsnp_id","chr_name","chrom_start","A1","A2","minor_allele","minor_allele_freq")]
colnames(dat) <- c("SNP_ID","CHR","POS","A1","A2","MAF","MAF_freq")
head(dat); dim(dat)

# SNP_ID CHR     POS A1 A2 MAF MAF_freq
# 1 rs2817134   1 3098715  T  A   T 0.210463
# 2 rs2817138   1 3105276  T  A   T 0.261581
# 3 rs2817148   1 3138136  A  G   G 0.176717
# 4 rs2817158   1 3149145  G  A   A 0.335663
# 5 rs2817168   1 3013690  A  C   A 0.339457
# 6 rs2817174   1 3044181  T  C   T 0.359625
# [1] 1222627       7

# Save the processed SNP data to a compressed file
fwrite(dat,"snp_sumstat.txt.gz",sep = "\t")


#------