# set directory
cd "/media/berheldefin/BERFIN/rna_seq/QoRTsPipelineWalkthrough - Kopya"
#reset everything
sh exampleScripts/step0/step0.reset.sh

#data processing and qc
while read line
do
mkdir outputData/qortsData/$line/
java -Xmx1G -jar softwareRelease/QoRTs.jar \
QC \
--stranded \
--chromSizes inputData/annoFiles/chrom.sizes \
inputData/bamFiles/$line.bam \
inputData/annoFiles/anno.gtf.gz \
outputData/qortsData/$line/
done < "inputData/annoFiles/uniqueID.list.txt"