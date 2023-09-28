cd "/media/berheldefin/BERFIN/rna_seq/QoRTsPipelineWalkthrough - Kopya"
#merge replicates
java -jar softwareRelease/QoRTs.jar \
mergeAllCounts \
outputData/qortsData/ \
inputData/annoFiles/decoder.byUID.txt \
outputData/countTables/

