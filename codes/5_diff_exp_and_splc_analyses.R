setwd ("/media/berheldefin/BERFIN/rna_seq/QoRTsPipelineWalkthrough - Kopya")
###differential expression analyses
##DESeq2
suppressPackageStartupMessages(library(DESeq2));
decoder.bySample <- read.table(
"inputData/annoFiles/decoder.bySample.txt",
header=T,stringsAsFactors=F);
directory <- "outputData/countTables/";
sampleFiles <- paste0(
decoder.bySample$sample.ID,
"/QC.geneCounts.formatted.for.DESeq.txt.gz"
);
sampleCondition <- decoder.bySample$group.ID;
sampleName <- decoder.bySample$sample.ID;
sampleTable <- data.frame(sampleName = sampleName,
fileName = sampleFiles,
condition = sampleCondition);
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
directory = directory,
design = ~ condition);
dds;
dds <- DESeq(dds);
res <- results(dds);
res
#save results
write.table(res, file="/media/berheldefin/BERFIN/rna_seq/QoRTsPipelineWalkthrough/outputData/analyses/DESeq/DESeq_results.tsv", sep="\t")

##edgeR
suppressPackageStartupMessages(library(edgeR));
decoder.bySample <- read.table(
"inputData/annoFiles/decoder.bySample.txt",
header=T,stringsAsFactors=F);
directory <- "outputData/countTables/";
files <- paste0(directory,
decoder.bySample$sample.ID,
"/QC.geneCounts.formatted.for.DESeq.txt.gz"
);
countData <- lapply(files, function(f){
ct <- read.table(f,header=F,stringsAsFactors=F)$V2;
ct <- ct[1:(length(ct)-5)]
});
countMatrix <- do.call(cbind.data.frame,countData);
colnames(countMatrix) <- decoder.bySample$sample.ID;
rownames(countMatrix) <- read.table(files[1],header=F,
stringsAsFactors=F)$V1[1:nrow(countMatrix)]
group <- factor(decoder.bySample$group.ID);
design <- model.matrix(~group)
y <- DGEList(counts = countMatrix, group = group);
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)
#save results ??????
lrt_matrix <- as.matrix(lrt)
lrt_df <- as.data.frame(lrt_matrix)
write.table(lrt_df, file="/media/berheldefin/BERFIN/rna_seq/QoRTsPipelineWalkthrough/outputData/analyses/edgeR/edgeR_results.txt", sep="\t")

###differential splicing analyses
##JunctionSeq
library("JunctionSeq");
#The sample decoder:
decoder <- read.table(
  "inputData/annoFiles/decoder.bySample.txt",
  header=T,stringsAsFactors=F);
#The count files:
countFiles <- paste0(
  "outputData/countTables/",
  decoder$sample.ID,
  "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz");
#Run the analysis:
jscs <- runJunctionSeqAnalyses(
  sample.files = countFiles,
  sample.names = decoder$sample.ID,
  condition = decoder$group.ID,
  flat.gff.file = "outputData/countTables/withNovel.forJunctionSeq.gff.gz",
  nCores = 1,
  verbose=TRUE,
  debug.mode = TRUE);
writeCompleteResults(jscs,
            outfile.prefix="outputData/analyses/JunctionSeq/",
            save.jscs = TRUE
);
writeCompleteResults(jscs,
outfile.prefix="outputData/analyses/JunctionSeq/",
save.jscs = TRUE
);

buildAllPlots(
  jscs=jscs,
  FDR.threshold = 0.01,
  use.plotting.device = "Cairo",
  outfile.prefix = "outputData/analyses/JunctionSeq/results/",
  variance.plot = TRUE,
  ma.plot = TRUE,
  rawCounts.plot=TRUE,
  verbose = TRUE);

##DEXSeq
suppressPackageStartupMessages(library(DEXSeq));
decoder <- read.table(
  "inputData/annoFiles/decoder.bySample.txt",
  header=T,stringsAsFactors=F);
directory <- "outputData/countTables/";
countFiles <- paste0(
  directory,
  decoder$sample.ID,
  "/QC.exonCounts.formatted.for.DEXSeq.txt.gz");
dexseq.gff <- "outputData/forDEXSeq.gff.gz";
sampleTable <- data.frame(
  row.names = decoder$sample.ID,
  condition = decoder$group.ID);
dxd <- DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData = sampleTable,
  design = ~sample + exon + condition:exon,
  flattenedfile=dexseq.gff);
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)
dxd <- testForDEU(dxd);
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition");
dxres <- DEXSeqResults(dxd);
dxres
#save results ????????