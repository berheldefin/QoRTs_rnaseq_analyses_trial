#Generate ”flat” annotation files for JunctionSeq
java -Xmx1G -jar softwareRelease/QoRTs.jar \
makeFlatGff \
--stranded \
inputData/annoFiles/anno.gtf.gz \
outputData/forJunctionSeq.gff.gz

#for DEXSeq
java -Xmx1G -jar softwareRelease/QoRTs.jar \
makeFlatGff \
--stranded \
--DEXSeqFmt \
inputData/annoFiles/anno.gtf.gz \
outputData/forDEXSeq.gff.gz

#adding novel functions 
java -Xmx1G -jar softwareRelease/QoRTs.jar \
mergeNovelSplices \
--minCount 6 \
--stranded \
outputData/countTables/ \
outputData/sizeFactors.GEO.txt \
inputData/annoFiles/anno.gtf.gz \
outputData/countTables/