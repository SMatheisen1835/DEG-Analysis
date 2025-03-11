# This is a pipeline that downloads expression data from NCBIs GEO and writes
# the top most DEGs into an excel file

# Version info: R 4.4.1, Biobase 2.64.0 , GEOquery 2.72.0, limma 3.60.6, 
#               umap 0.2.10.0, openxlsx 4.2.7.1

# parameters to be set:
GSE_ID <- "GSE68801"


################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)
library(openxlsx)

# load series and platform data from GEO

# this line downloads the data from the database into the work environment
# the parameter GSEMatrix = TRUE means that the data undergoes some pre-processing
# so that it is ready to use
# AnnotGPL = TRUE adds annotations such as gene symbols, names and descriptions
# to the matrix
gset <- getGEO(GSE_ID, GSEMatrix =TRUE, AnnotGPL=TRUE)

# some of GEO's data sets have different platforms. The GLP IDs represent these.
# There are different IDs for example for a specific array chip or sequencing technique.
# If there is data from multiple platforms, this code makes sure that we only
# use the data corresponding to our GPL-ID
# Check if there are multiple platforms in the data set
if (length(gset) > 1) {
  platforms <- attr(gset, "names")
  print(platforms)
} else {
  print("Only one platform included")
}
# if the last code block printed out multiple GLP IDs choose one. Looking at 
# the GEO entry in the data base might help decide which GLP to use.
GLP_ID <- ""
if (length(gset) > 1) idx <- grep(GLP_ID, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# This ensures that the column names in gset's feature data are compatible with 
# subsequent analysis steps, such as matching column names in the topTable output 
# or using column names programmatically. Many R functions (e.g., subsetting or 
# plotting) require column names to follow specific syntax rules.
fvarLabels(gset) <- make.names(fvarLabels(gset))

# Now we need to select which samples belong to which group. Therefore here is
# code that will view a table similar to how it would look in GEO2R
# first we extract the phenotype data
view <- pData(gset)
# we select which columns we want to look at. To view which are available use:
# colnames(pData(gset))
view <- subset(view, select=c("title","geo_accession", "characteristics_ch1", "source_name_ch1", "age:ch1", "gender:ch1"))
View(view)
# if you have decided how many groups and which groups (like control and case) 
# you want to have input them here:
groups <- make.names("control", "case")
# to decide which sample belongs to which group, put a 0 for one and a 1 for 
# another group. The numbers should correspond to the groups you made earlier 
assign.groups <- "1111111100001111111111100001111111100001111111111111000011110"
# running these next 3 lines will add the groups to the view so you can check if
# everything is correct
sml <- strsplit(assign.groups, split="")[[1]]
view$groups <- sml
View(view)
# if there is something wrong, you can edit it with the window that pops up here:
edit(view)
sml <- view$groups
paste(view$groups, collapse = "")


# log2 transformation
# first we extract the data into an expression matrix ex
ex <- exprs(gset)
# now we collect some statistics about or data so we can decide how to handle it.
# we get the different quantiles telling us the min, max, etc. values
# na.rm=T ignores not available (missing) data 
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# next we will decide if we will need to make a log2 transformation. If LogC is 
# true we will transform the expression values. If LogC = FALSE it suggests that
# the data might already be normalized (e.g. already log-transformed)
# LogC is TRue if the 99% quantile is bigger than 100 OR
# LogC is TRUE if the distance between the min and max value is greater than 50 AND
# the 25% quantile is greater than 0
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
# if LogC is TRUE then we will now do the log2-transformation
# since log2 of 0  and of negative numbers is undefined, we will set all of those
# expressions to NaN (Not a Number). All other values undergo log2 transformation
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
# sml is the list of characters indicating which sample belongs to which group
# first we transform this list into a factor
gs <- factor(sml)
# make.names() makes sure that our names a viable
# for example it removes spaces and other problematic characters
groups <- make.names(c("control","case"))
# now we change the levels in gs from 0 and 1 to control and case
levels(gs) <- groups
# the group information is added the gset matrix by adding a new column gs 
gset$group <- gs
# now we create the design matrix that is needed for analysis by limma for example
# ~group + 0 is a formula that specifies the way the design matrix should be created
design <- model.matrix(~group + 0, gset)
# the column names of the design matrix should also have the same name as our 
# levels (which are case and control here)
colnames(design) <- levels(gs)

# here we remove any samples(columns) that contain missing values (NA or NaN)
gset <- gset[complete.cases(exprs(gset)), ]

# now we use the limma package to fit a linear model on the expression data (gset)
# using the design matrix
fit <- lmFit(gset, design)

# set up contrasts of interest and recalculate model coefficients
# in cts we specify for later which groups we want to compare. Here we only 
# have 2 groups but in case we had more than 2 we would have to specify which we
# want to use here  
cts <- paste(groups[1], groups[2], sep="-")
# here we create the contrast matrix that we need so we can specify which groups
# we want to compare in our differential expression analysis
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
# here we do another fit. This will help us compare the control with the case group
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
# here we apply empirical Bayes moderation. It stabilizes the variance estimates
# and prevents extreme t-statistics or p-values that might result from small 
# sample sizes or noisy data.
fit2 <- eBayes(fit2, 0.01)
# here we extract the top "number" of expressed genes sorted by the B value
tT <- topTable(fit2, adjust="fdr", sort.by="B", number = Inf)
# another possible query would be 
# topGenes <- topTable(fit2, adjust="BH", sort.by="P", number=10)
# where we sort by P and adjust by Benjamin-Hochberg


# here we choose which columns we want to keep
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GenBank.Accession","Gene.symbol","Gene.title"))
# here we write our table to a file
write.table(tT, file=paste0(GSE_ID,".tsv"), row.names=F, sep="\t")

write.xlsx(tT, paste0(GSE_ID,".xlsx"), rowNames = FALSE)
