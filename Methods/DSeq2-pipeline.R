# Pipeline for DEG analysis of RNA-seq transcriptome data using DESeq2

# Version info: R 4.2.2, Biobase 2.64.0 , GEOquery 2.66.0, DESeq2 1.38.3, 
#               openxlsx 4.2.7.1, rstudioapi 0.17.1

# parameters to be set:
GSE_ID <- "GSE161245"

################################################################
#   Differential expression analysis with DESeq2

#load needed libraries
library(DESeq2)
#library(R.utils)
library(GEOquery)
library(rstudioapi)
library(openxlsx)

# load series and platform data from GEO

# setting the work directory
script_dir <- dirname(getActiveDocumentContext()$path)
setwd(script_dir)
getwd()

# downloads the GEO files and puts them into a folder in the work directory
getGEOSuppFiles(GSE_ID)
# inside the folder is a single file with all files zipped in .tar format so we
# unzip the files to a new folder
untar(paste0(GSE_ID, "/", GSE_ID, "_RAW.tar"), exdir= paste0(GSE_ID, "_files"))

# now there is a folder with files that are still gunzipped but we can just read
# the files in normally before manually decompressing them
files <- list.files(paste0(GSE_ID, "_files"), pattern = "\\.gz$", full.names = TRUE)  # List all .gz files
#print(files)  # View the list of files

# Read each .gz file containing the count data and put them in a list
counts_list <- lapply(files, function(file) {
  read.table(file, header = TRUE, sep = "\t", row.names = 1)  # Load each .gz file into a list
})

# Check dimensions and structure of the first count matrix to ensure it's as expected
dim(counts_list[[1]])  # Check dimensions of the first file
head(counts_list[[1]])  # Preview first few rows of the first count matrix

# looking at the first sample we see a different sample name that is not the GSM number
# therefore we exchange the colnames with the actual GSM numbers

# we can extract the GSM numbers from the file names:
GSMs <- substr(files, 17, 26)

# and then change the colnames to their GSM number
for(i in 1:length(GSMs)){
  colnames(counts_list[[i]]) <- GSMs[i]
}

# this should now show the GSM number as the colname
head(counts_list[[1]])  # Preview first few rows of the first count matrix



# Extract the gene names from the first count matrix
gene_names <- rownames(counts_list[[1]])

# Initialize an empty list to store the combined matrices (by sample)
counts_combined <- counts_list[[1]]

# Iterate over the remaining count matrices
for (i in 2:length(counts_list)) {
  # Check if gene names match
  if (!all(gene_names == rownames(counts_list[[i]]))) {
    warning(paste("Gene names do not match between file 1 and file", i))
    # Optionally, you can identify the differences:
    mismatched_genes <- setdiff(gene_names, rownames(counts_list[[i]]))
    cat("Mismatched genes between files:", mismatched_genes, "\n")
    break  # Exit the loop early if there's a mismatch
  }
  
  # Use cbind to add the current count matrix to the combined matrix
  counts_combined <- cbind(counts_combined, counts_list[[i]])
}

# Step 7: Check the Combined Count Matrix
dim(counts_combined)  # Check dimensions (genes x samples)
head(counts_combined)  # Preview first few rows

# Step 8: Verify if the Data is Raw Counts (integers expected)
# summary(counts_combined)  # Check for raw counts (should be integers)
is.integer(as.matrix(counts_combined))  # Should return TRUE if raw counts


# Load the GEO dataset
gse <- getGEO(GSE_ID, GSEMatrix = TRUE)  # Set GSEMatrix=TRUE to fetch data as a matrix

# extract the metadata from the data set
sample_metadata <- pData(gse[[1]])
# there is a lot of information here but we are only interested in the donor condition
head(sample_metadata)

# here we only get the col with the donors condition, overwrite the col name 
# and make sure that the conditions are valid names (no spaces etc.) 
metadata <- sample_metadata["donor:ch1"]
colnames(metadata) = "condition"
metadata$condition <- make.names(metadata$condition)

# now all preprocessing is finished and we have all we need for DESeq2
head(metadata)
head(counts_combined)


# Now we can create the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts_combined, colData = metadata, design = ~ condition)
dds 

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# this step should have filtered out about half of the genes (seen in the dims)
dds

# make sure the healthy controls are used as the references
dds$condition <- relevel(dds$condition, ref = "Healthy.controls")

# Now we can start the DEG analysis with DESeq2. This step may take some time
dds <- DESeq(dds)
dds

# now for the results:
res <- results(dds, alpha = 0.05)
summary(res)

# contrasts
resultsNames(dds)
summary(res[2])
summary(res[3])
summary(res[4])

# to look at specific results:
res.mild <- results(dds, contrast = c("condition", "Mild.asthma.patients",
                                      "Healthy.controls"), alpha = 0.05)
summary(res.mild)
head(res.mild[order(res.mild$padj), ])


res.moderate <- results(dds, contrast = c("condition", "Moderate.asthma.patients",
                                          "Healthy.controls"), alpha = 0.05)
summary(res.moderate)
head(res.moderate[order(res.moderate$padj), ])

res.severe <- results(dds, contrast = c("condition", "Severe.asthma.patients", "Healthy.controls"), alpha = 0.05)
summary(res.severe)
head(res.severe[order(res.severe$padj), ])

# save the DEG tables as an excel file
res.mild <- res.mild[which(res.mild$padj<=0.05), ]
write.xlsx(res.mild, paste0("../../Results/", GSE_ID, "/", GSE_ID, "_mild.xlsx"))
res.moderate <- res.moderate[which(res.moderate$padj<=0.05), ]
write.xlsx(res.moderate, paste0("../../Results/", GSE_ID, "/", GSE_ID, "_moderate.xlsx"))
res.severe <- res.severe[which(res.severe$padj<=0.05), ]
write.xlsx(res.severe, paste0("../../Results/", GSE_ID, "/", GSE_ID, "_severe.xlsx"))

