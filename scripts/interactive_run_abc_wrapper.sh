# pwd
cd /home/rad/users/gaurav/projects/abc

# Running the ABC Model

# Running the ABC model consists of the following steps:
# 1) Define candidate enhancer regions
# 2) Quantify enhancer activity
# 3) Compute ABC Scores

# Step 1. Define candidate elemets
# 'Candidate elements' are the set of putative enhancer elements for which ABC Scores will be computed. A typical way to define candidate elements is by calling peaks on a DNase-Seq or ATAC-Seq bam file. In this implementation we first call peaks using MACS2 and then process these peaks using makeCandidateRegions.py.
# makeCandidateRegions.py will take as input the narrowPeak file produced by MACS2 and then perform the following processing steps:
# 1.1) Count DNase-seq reads in each peak and retain the top N peaks with the most read counts
# 1.2) Resize each of these N peaks to be a fixed number of base pairs centered on the peak summit
# 1.3) Remove any blacklisted regions and include any whitelisted regions
# 1.4) Merge any overlapping regions
# We recommend using --nStrongestPeaks 150000 when making genome-wide peak calls.

# Step 2. Quantifying Enhancer Activity:
# run.neighborhoods.py will count DNase-seq (or ATAC-seq) and H3K27ac ChIP-seq reads in candidate enhancer regions. It also makes GeneList.txt, which counts reads in gene bodies and promoter regions.
# Replicate epigenetic experiments should be included as comma delimited list of files. Read counts in replicate experiments will be averaged when computing enhancer Activity.
# Main output files:
# - EnhancerList.txt: Candidate enhancer regions with Dnase-seq and H3K27ac ChIP-seq read counts
# - GeneList.txt: Dnase-seq and H3K27ac ChIP-seq read counts on gene bodies and gene promoter regions

# Step 3. Computing the ABC Score
# Compute ABC scores by combining Activity (as calculated by run.neighborhoods.py) and Hi-C.

# # The main output files are:
# - EnhancerPredictions.txt: all element-gene pairs with scores above the provided threshold. Only includes expressed genes and does not include promoter elements. This file defines the set of 'positive' predictions of the ABC model.
# - EnhancerPredictionsFull.txt: same as above but includes more columns. See https://docs.google.com/spreadsheets/d/1UfoVXoCxUpMNPfGypvIum1-RvS07928grsieiaPX67I/edit?usp=sharing for column definitions
# - EnhancerPredictions.bedpe: Same as above in .bedpe format. Can be loaded into IGV.
# - EnhancerPredictionsAllPutative.txt.gz: ABC scores for all element-gene pairs. Includes promoter elements and pairs with scores below the threshold. Only includes expressed genes. This file includes both the 'positive' and 'negative' predictions of the model. (use --make_all_putative to generate this file).
# - EnhancerPredictionsAllPutativeNonExpressedGenes.txt.gz: Same as above for non-expressed genes. This file is provided for completeness but we generally do not recommend using these predictions.
# The default threshold of 0.02 corresponds to approximately 70% recall and 60% precision [1].


# 1. Define candidate elemets
# narroPeakFile='/media/rad/HDD1/atacseq/christine/ckatac/results/bwa/mergedLibrary/macs/broadPeak/53646_PPT-1_S1_R1.mLb.clN_peaks.broadPeak'
# sortedNarroPeakFile='/media/rad/HDD1/atacseq/christine/ckatac/results/bwa/mergedLibrary/macs/broadPeak/53646_PPT-1_S1_R1.mLb.clN_peaks.broadPeak.sorted'
# atacBamFile='/media/rad/HDD1/atacseq/christine/ckatac/results/bwa/mergedLibrary/53646_PPT-1_S1_R1.mLb.clN.sorted.bam'
narroPeakFile='/media/rad/HDD1/atacseq/christine/ckatac/gjatac/peaks/53646_PPT-1_S1_R1_001_rmdup/macs2peaks/53646_PPT-1_S1_R1_001_rmdup_peaks.broadPeak'
sortedNarroPeakFile="${narroPeakFile}.sorted"
atacBamFile='/media/rad/HDD1/atacseq/christine/ckatac/gjatac/bams/trimmed/53646_PPT-1_S1_R1_001_rmdup.bam'
sampleName="53646_PPT-1_S1_R1"
outputDir="/media/rad/HDD1/abc/christine/pdacBatch1/${sampleName}"
chromSizesFile='/home/rad/users/gaurav/projects/abc/input/reference/mm10_ucsc.chromsizes'
blacklistFile='/home/rad/users/gaurav/projects/abc/input/reference/GRCm38-blacklist_ucsc.bed'
nStrongestPeaks=150000

# Sort narrowPeak file
bedtools sort -faidx ${chromSizesFile} -i ${narroPeakFile} > ${sortedNarroPeakFile}
python scripts/ABC-Enhancer-Gene-Prediction/src/makeCandidateRegions.py --narrowPeak ${sortedNarroPeakFile} --bam ${atacBamFile} --outDir ${outputDir} --chrom_sizes ${chromSizesFile} --regions_blacklist ${blacklistFile} --peakExtendFromSummit 250 --nStrongestPeaks ${nStrongestPeaks}  --ignoreSummits

# 2. Quantifying Enhancer Activity
candidateRegionsFile=${outputDir}/$(basename ${sortedNarroPeakFile}).candidateRegions.bed
genesBedFile="input/reference/mmu_GRCm38_gencode_vM24_genes_ucsc.bed"
outputNeighborDir="/media/rad/HDD1/abc/christine/pdacBatch1/${sampleName}/neighbors"
h3k27acBamFile="/media/rad/HDD1/nfchip/christine/pdacBatch1/gjchip/mapping/pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT-1_mmu_chipseq_se_R1_rmdup.bam"
python scripts/ABC-Enhancer-Gene-Prediction/src/run.neighborhoods.py --candidate_enhancer_regions ${candidateRegionsFile} --genes ${genesBedFile} --H3K27ac ${h3k27acBamFile} --DHS ${atacBamFile} --chrom_sizes ${chromSizesFile} --outdir ${outputNeighborDir}
