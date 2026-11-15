setwd("~/github/private/coursera")
file <- "BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"

# cell #1
mrnaIDs <- read.table(file, 
                      header = F, fill = T, nrows = 1)
mrnaIDs <- mrnaIDs[, -1][, -1]  # remove the first solumn "Hybridization REF" which spans two columns because of the space between the two words

# cell #2
mrnaNorm <- read.delim(file, skip = 2, header = FALSE)
colnames(mrnaNorm) <- c("gene", unlist(mrnaIDs))

# cell #4
samp <- lapply(as.list(t(mrnaIDs)), function(t) substr(unlist(strsplit(t, "-"))[4], 1, 2))
# extract the sample type (tumor / normal), excluding the first column since it is a gene ID
sampleType <- as.data.frame(samp)
dim(sampleType)
# 1 1212

# cell #5
unique(sampleType)
# extracts how many unique objects there are 

tab <- table(unlist(sampleType))
tab
# count how many of each type
#   01   06   11 
#  1093   7  112 

# cell #6
sampClass <- lapply(samp, function(t) (if (t < 10) return("1") else return("0")))
mrnaClass <- as.data.frame(sampClass)
dim(mrnaClass)
table(unlist(sampClass))
#   0    1 
# 112 1100 

# cell #7
sampClassNum <- lapply(samp, function(t) (if (t < 10) return(1) else return(0)))
mrnaClassNum <- as.data.frame(sampClassNum) 

# cell #8
bssWssFast <- function (X, givenClassArr, numClass=2)
# between squares / within square feature selection
{
	classVec <- matrix(0, numClass, length(givenClassArr))
	for (k in 1:numClass) {
		temp <- rep(0, length(givenClassArr))
		temp[givenClassArr == (k - 1)] <- 1
		classVec[k, ] <- temp
	}
	classMeanArr <- rep(0, numClass)
	ratio <- rep(0, ncol(X))
	for (j in 1:ncol(X)) {
		overallMean <- sum(X[, j]) / length(X[, j])
		for (k in 1:numClass) {
			classMeanArr[k] <- 
				sum(classVec[k, ] * X[, j]) / sum(classVec[k, ])
		}
	  classMeanVec <- classMeanArr[givenClassArr + 1]
	  bss <- sum((classMeanVec - overallMean)^2)
	  wss <- sum((X[, j] - classMeanVec)^2)
	  ratio[j] <- bss/wss
	}
	    list(x = sort(ratio, decreasing = TRUE), ix = order(ratio, decreasing = TRUE))
}

# cell #9
bss <- bssWssFast(t(mrnaNorm[, -1]), t(mrnaClassNum), 2)
# returns a list with two elements: a list of bss/wss ratios in decreasing order bss[1] and a list of features bss[2]

bss$ix[1:50]
# show list of 50 ranked gene indexes

# cell #10
genes <- mrnaNorm[bss$ix[1:50],1]
genes
# show list of 50 tanked gene names
# [1] FHL1|2273        ADAMTS4|9507     ABAT|18          PAG1|55824       ARHGAP1|392     
# [6] TP53TG3B|729355  SCRN3|79634      TIE1|7075        SNRPA|6626       SCCPDH|51097    
#[11] PCYT1B|9468      DSE|29940        LGI4|163175      LIX1L|128077     LCE3D|84648     
#[16] POTEE|445582     GIMAP2|26157     CABYR|26256      PLGLA|285189     LRRC57|255252   
#[21] C21orf33|8209    CDC14A|8556      CD300LF|146722   CEP57|9702       MAGEC2|51438    
#[26] LGALS3|3958      CAT|847          LOC647979|647979 ABCA8|10351      PPIL3|53938     
#[31] EFCAB1|79645     CCDC45|90799     ABCA5|23461      CAV1|857         TP53I3|9540     
#[36] ADRB1|153        SOD1|6647        LDLRAP1|26119    LYPLA1|10434     KLHL10|317719   
#[41] HLA-DPB1|3115    SGCG|6445        MLL2|8085        FGD4|121512      CD33|945        
#[46] BMS1|9790        MARCKSL1|65108   BTNL8|79908      SS18L2|51188     ABCA7|10347   


# cell #11
# merge class values into mrnaNorm as the last row - after the genes and select only top 100 genes
mrnaSetClass <- rbind(mrnaNorm[bss$ix[1:100],-1], setNames(mrnaClassNum, names(mrnaNorm[,-1])))
# name the rows with gene symbols ("FIGF|2277" -> "FIGF"); the last row is the class
rownames(mrnaSetClass) <- c(sub("\\|.*", "", mrnaNorm[bss$ix[1:100], 1]), "class")
dim(mrnaSetClass)
# 101  1212       101 genes as rows and 1212 patients as columns
# to select by class, transpose the matrix
transpose <- as.data.frame(t(mrnaSetClass))
# select normals
normals <- subset(transpose, transpose$class == 0)
dim(normals)
#[1] 112 101
# select tumors
tumors <- subset(transpose, transpose$class == 1)
dim(tumors)
#[1] 1100 101
# combine them together
both <- rbind(normals, tumors[1:112,])
dim(both)
# select a color scheme - red for normal and blue for tumor
color.map <- function(class) { if (class==0) "#0000FF" else "#FF0000" } # red cancer, blue for normal (1)
# create heatmap (the "class" column is left out so it is not drawn as a gene)
groupColors <- unlist(lapply(both$class, color.map))

png("heatmap_top100.png", width = 1000, height = 1400, res = 300)
heatmap(as.matrix(t(both[, colnames(both) != "class"])), scale = "row", col = topo.colors(100),
        ColSideColors = groupColors, labCol = FALSE, cexRow = 0.5)
dev.off()

#=====================================
# Next script
# ====================================

library(randomForest)
library(class)

geneNames <- mrnaNorm[1] # extract the gene names from mrnaNorm as its first column
dim(geneNames)

gc()

# cell #6
trainSet <- mrnaData
testSet <- mrnaData
trainClasses <- unlist(mrnaClassNum[1,], use.names=FALSE)
testClasses <- unlist(mrnaClassNum[1,], use.names=FALSE)
knn.predic <- knn(trainSet, testSet, trainClasses, testClasses,k=1)
cbr.predic = as.vector(knn.predic)
table(cbr.predic, testClasses)
tab <- table(cbr.predic, t(testClasses))
error <- sum(tab) - sum(diag(tab))
accuracy <- round(100- (error * 100 / length(testClasses)))
print(paste("accuracy= ", as.character(accuracy), "%"), quote=FALSE)

# cell #8
# select features
dim(mrnaData)
# 1212 20531  matrix
dim(mrnaClass)
# 1 1212
dim(mrnaClassNum)
# 1 1212
dim(geneNames)
# 20531 genes
bss <- bssWssFast(mrnaData, t(mrnaClassNum), 2)

# cell #9
mrnaDataReduced <- mrnaData[,bss$ix[1:100]]
dim(mrnaDataReduced)
# 1212  100
trainSet <- mrnaDataReduced
testSet <- mrnaDataReduced
trainClasses <- unlist(mrnaClassNum[1,], use.names=FALSE)
# or as.numeric(mrnaClassNum[1,])
testClasses <- unlist(mrnaClassNum[1,], use.names=FALSE)

# cell #10
knn.predic <- knn(trainSet, testSet, trainClasses, testClasses,k=1) # knn form 'class' package
knn.predic = as.vector(knn.predic)  # change knn.predic to become a vector
table(knn.predic, testClasses)      # build the confusion matrix
tab <- table(knn.predic, t(testClasses))
error <- sum(tab) - sum(diag(tab))  # calculate acuracy
accuracy <- round(100- (error * 100 / length(testClasses)))
print(paste("accuracy= ", as.character(accuracy), "%"), quote=FALSE)   # display acuracy after formating it as a character string

# cell #11
trainSetClass <- as.data.frame(cbind(trainSet, t(mrnaClassNum[1,])))  # concatenate gene expressions and class data
testSetClass <- as.data.frame(cbind(testSet, t(mrnaClassNum[1,])))    # concatenate gene expressions and class data
colnames(trainSetClass)[101] <- "class"     # give a name to the class column
#trainSetClass$class <- as.numeric(trainSetClass$class) # for regression
trainSetClass$class <- as.factor(trainSetClass$class)  # for classification
class(trainSetClass$class)      # should be factor or categorical for classification
rf <- randomForest(class ~., trainSetClass,
                   ntree=100,
                   importance=T)      # build randomForest classifier
colnames(testSetClass)[101] <- "class"     # give a name to the class column
testSetClass$class <- as.factor(testSetClass$class)  # for classification
rf.predic <- predict(rf ,testSetClass)  # test the randomForest built model on the test set
rf.predic = as.vector(rf.predic)        # change rf.predic to become a vector
table(rf.predic, testClasses)           # build the confusion matrix
tab <- table(rf.predic, t(testClasses))
error <- sum(tab) - sum(diag(tab))      # calculate acuracy
accuracy <- round(100- (error * 100 / length(testClasses)))
print(paste("accuracy= ", as.character(accuracy), "%"), quote=FALSE)

# cell #12
nbRows <- nrow(mrnaDataReduced)
set.seed(33)       # seet random seed so that we always get same samples drawn - since they are random
trainRows <- sample(1:nbRows, .70*nbRows)
trainSet <- mrnaDataReduced[trainRows, ]
testSet <- mrnaDataReduced[-trainRows, ]
dim(trainSet)
dim(testSet)

# cell #13
trainClasses <- unlist(mrnaClassNum[1,trainRows], use.names=FALSE)
testClasses <- unlist(mrnaClassNum[1,-trainRows], use.names=FALSE)
knn.predic <- knn(trainSet, testSet, trainClasses, testClasses,k=1)
knn.predic = as.vector(knn.predic)
table(knn.predic, testClasses)
tab <- table(knn.predic, t(testClasses))
error <- sum(tab) - sum(diag(tab))
accuracy <- round(100- (error * 100 / length(testClasses)))
print(paste("accuracy= ", as.character(accuracy), "%"), quote=FALSE)

# cell #14
trainSetClass <- as.data.frame(cbind(trainSet, t(mrnaClassNum[1,trainRows])))
testSetClass <- as.data.frame(cbind(testSet, t(mrnaClassNum[1,-trainRows])))
colnames(trainSetClass)[101] <- "class"
trainSetClass$class <- as.factor(trainSetClass$class)  # for classification
class(trainSetClass$class)
# should be factor for classification
rf <- randomForest(class ~., trainSetClass,
                   ntree=100,
                   importance=T)
colnames(testSetClass)[101] <- "class"
testSetClass$class <- as.factor(testSetClass$class)  # for classification
rf.predic <- predict(rf ,testSetClass)
rf.predic = as.vector(rf.predic)
table(rf.predic, testClasses)
tab <- table(rf.predic, t(testClasses))
error <- sum(tab) - sum(diag(tab))
accuracy <- round(100- (error * 100 / length(testClasses)))
print(paste("accuracy= ", as.character(accuracy), "%"), quote=FALSE)

#=====================================
# Next script
# ====================================


# cell #1
file2 <-"BRCA.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg"
cnvLogs <- read.table(file2, 
                      header = T, fill = T)
dim(cnvLogs)
# 284458 6
head(cnvLogs)

# cell #2
summary(cnvLogs)
# mean = -0.1132 median = 0

# cell #3
nrow(subset(cnvLogs, Sample == "TCGA-A8-A06R-01A-11D-A011-01"))
# 968

# cell #4
nrow(subset(cnvLogs, Sample == "TCGA-A8-A06R-01A-11D-A011-01" & Segment_Mean > 0))
# 585

# cell #5
nrow(subset(cnvLogs, Sample == "TCGA-A8-A06R-01A-11D-A011-01" & Segment_Mean < 0))
# 383

# cell #6
mean(subset(cnvLogs, Sample == "TCGA-A8-A06R-01A-11D-A011-01")[["Segment_Mean"]])
# 0.201533574380165

#=====================================
# Next script
# ====================================
library(enrichplot)
library(cluster)
library(ReactomePA)


# cell #5
mrnaData = t(mrnaNorm[, -1])

# cell #7
# select features

bss <- bssWssFast(mrnaData, t(mrnaClassNum), 2)
mrnaDataReduced <- mrnaData[,bss$ix[1:100]]
dim(mrnaDataReduced)
# 1212  100
trainClasses <- unlist(mrnaClassNum[1,], use.names=FALSE)

# cell #8
set.seed(20)
kmeans.clusters <- kmeans(mrnaDataReduced, 2, nstart = 20)
table(kmeans.clusters$cluster)

# cell #9
table(kmeans.clusters$cluster, trainClasses)

# cell #10
png("clusplot.png", width = 2000, height = 2800, res = 300)
clusplot(mrnaDataReduced, kmeans.clusters$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)
dev.off()

# cell #11
pam.clusters <- pam(mrnaDataReduced, 2)
table(pam.clusters$clustering)

# cell #12
table(pam.clusters$clustering, trainClasses)

# cell #13
png("clusplot2.png", width = 2000, height = 2800, res = 300)
clusplot(mrnaDataReduced, pam.clusters$clustering, color=TRUE, shade=TRUE, 
         labels=2, lines=0)
dev.off()

# cell #14
genes_f <- geneNames[bss$ix[1:100],1]
genes_f <- as.character(genes_f)
hugoNames <- lapply(genes_f, function(t) substr(t, 1, regexpr("\\|", t) - 1))   # extract HUGO names
entrezNames <- lapply(genes_f, function(t) substr(t, regexpr("\\|", t) + 1, nchar(t)))  # extract Entrez names
paths = enrichPathway(unlist(entrezNames), pvalueCutoff=1)


# cell #15
paths <- pairwise_termsim(paths)
png("ema.png", width = 2000, height = 2800, res = 300)
emapplot(paths, showCategory = 30, layout = igraph::layout_with_kk, node_label_size = 3)
dev.off()
