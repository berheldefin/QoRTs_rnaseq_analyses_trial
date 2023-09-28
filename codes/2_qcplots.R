#generating qc plots
setwd ("/media/berheldefin/BERFIN/rna_seq/QoRTsPipelineWalkthrough - Kopya")
library(QoRTs);
#Read in the QC data:
res <- read.qc.results.data("outputData/qortsData/",
decoder.files = "inputData/annoFiles/decoder.byUID.txt",
calc.DESeq2 = TRUE, calc.edgeR = TRUE)
makeMultiPlot.all(res,
outfile.dir = "outputData/qortsPlots/summaryPlots/",
plot.device.name = "CairoPNG");
makeMultiPlot.all(res,
outfile.dir = "outputData/qortsPlots/summaryPDFs/",
plot.device.name = "pdf");
makeMultiPlot.basic(res,
outfile.dir = "outputData/qortsPlots/basicPlots/",
plot.device.name = "CairoPNG",
separatePlots = TRUE);

#get size factors
get.size.factors(res, outfile = "outputData/sizeFactors.GEO.txt");
