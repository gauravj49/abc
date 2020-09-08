#!/bin/bash

# USAGE: bash scripts/run_abc_wrapper.sh '/media/rad/HDD1/atacseq/christine/ckatac/gjatac/peaks/53646_PPT-1_S1_R1_001_rmdup/macs2peaks/53646_PPT-1_S1_R1_001_rmdup_peaks.broadPeak' '/media/rad/HDD1/atacseq/christine/ckatac/gjatac/bams/trimmed/53646_PPT-1_S1_R1_001_rmdup.bam' '/media/rad/HDD1/nfchip/christine/pdacBatch1/gjchip/mapping/pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT-1_mmu_chipseq_se_R1_rmdup.bam' '53646_PPT-1_S1_R1' '/media/rad/HDD1/abc/christine/pdacBatch1/test' 'mouse'

# Parameters for the script
atacPeakFile=${1:-""}
atacBamFile=${2:-""}
h3k27acBamFile=${3:-""}
sampleName=${4:-""}
outputDir=${5:-""}
species=${6:-""}
nStrongestPeaks=${7:-"150000"}
hicdir=${8:-"/home/rad/users/gaurav/projects/abc/input/hicBedpe"}
blacklistFile=${9:-"/home/rad/users/gaurav/projects/abc/input/reference/GRCm38-blacklist_ucsc.bed"}
chromSizesFile=${10:-"/home/rad/users/gaurav/projects/abc/input/reference/mm10_ucsc.chromsizes"}
genesBedFile=${11:-"/home/rad/users/gaurav/projects/abc/input/reference/mmu_GRCm38_gencode_vM24_genes_ucsc.bed"}
jobdir=${12:-"/home/rad/users/gaurav/projects/abc"}

echo
echo "################# COMMAND LINE ARGUMENTS ######################"
echo "
atacPeakFile=${1:-""}
atacBamFile=${2:-""}
h3k27acBamFile=${3:-""}
sampleName=${4:-""}
outputDir=${5:-""}
species=${6:-""}
nStrongestPeaks=${7:-"150000"}
hicdir=${8:-"/home/rad/users/gaurav/projects/abc/input/input/hicBedpe"}
blacklistFile=${9:-"/home/rad/users/gaurav/projects/abc/input/reference/GRCm38-blacklist_ucsc.bed"}
chromSizesFile=${10:-"/home/rad/users/gaurav/projects/abc/input/reference/mm10_ucsc.chromsizes"}
genesBedFile=${11:-"/home/rad/users/gaurav/projects/abc/input/reference/mmu_GRCm38_gencode_vM24_genes_ucsc.bed"}
jobdir=${12:-"/home/rad/users/gaurav/projects/abc"}
"
echo "##########################################################"
echo


echo "# Sort narrowPeak file"
# Sort narrowPeak file
sortedAtacPeakFile="${atacPeakFile}.sorted"
bedtools sort -faidx ${chromSizesFile} -i ${atacPeakFile} > ${sortedAtacPeakFile}
echo 

# 1. Define candidate elemets
echo "# 1. Define candidate elemets"
python scripts/ABC-Enhancer-Gene-Prediction/src/makeCandidateRegions.py --narrowPeak ${sortedAtacPeakFile} --bam ${atacBamFile} --outDir ${outputDir} --chrom_sizes ${chromSizesFile} --regions_blacklist ${blacklistFile} --peakExtendFromSummit 250 --nStrongestPeaks ${nStrongestPeaks}  --ignoreSummits
echo 

# 2. Quantifying Enhancer Activity
echo "# 2. Quantifying Enhancer Activity"
candidateRegionsFile=${outputDir}/$(basename ${sortedAtacPeakFile}).candidateRegions.bed
outputNeighborDir="${outputDir}/${sampleName}/neighbors"
python scripts/ABC-Enhancer-Gene-Prediction/src/run.neighborhoods.py --candidate_enhancer_regions ${candidateRegionsFile} --genes ${genesBedFile} --H3K27ac ${h3k27acBamFile} --DHS ${atacBamFile} --chrom_sizes ${chromSizesFile} --outdir ${outputNeighborDir}
echo 

# 3. Compute ABC Scores
echo "# 3. Compute ABC Scores"
enhancerList=${outputNeighborDir}/EnhancerList.txt
geneList=${outputNeighborDir}/GeneList.txt
outputPredictionsDir="${outputDir}/${sampleName}/predictions"
# python scripts/ABC-Enhancer-Gene-Prediction/src/predict.py --enhancers ${enhancerList} --genes ${geneList} --score_column powerlaw.Score --threshold .02 --hic_resolution 5000 --outdir ${outputPredictionsDir} --make_all_putative
python scripts/ABC-Enhancer-Gene-Prediction/src/predict.py --enhancers ${enhancerList} --genes ${geneList} --HiCdir ${hicdir} --hic_resolution 70000  --scale_hic_using_powerlaw --hic_type bedpe --threshold .02 --outdir ${outputPredictionsDir} --make_all_putative

echo ""

#################### DOCUMENTATION ########################
# # pwd
# cd /home/rad/users/gaurav/projects/abc

# Source: https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction

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

# Contact data in bedpe formats
# Contact data can also be loaded in a generic bedpe format. This supports promoter-capture Hi-C and mixed-resolution Hi-C matrices. Use --hic_type bedpe in predict.py to enable this. Note that if contact data is provided in bedpe format, some of the default Hi-C normalizations are not applied. (For example, the code will not adjust the diagonal of the contact matrix)

# The bedpe file should be a tab delimited file containing 8 columns (chr1,start1,end1,chr2,start2,end2,name,score) where score denotes the contact frequency. An example is given here: input_data/HiC/raw/chr22/chr22.bedpe.gz

# python src/predict.py \
# --enhancers example_chr22/ABC_output/Neighborhoods/EnhancerList.txt \
# --genes example_chr22/ABC_output/Neighborhoods/GeneList.txt \
# --HiCdir example_chr22/input_data/HiC/raw/ \
# --hic_type bedpe \
# --hic_resolution 5000 \
# --scale_hic_using_powerlaw \
# --threshold .02 \
# --cellType K562 \
# --outdir example_chr22/ABC_output/Predictions/ \
# --make_all_putative
# ABC model without experimental contact data
# If experimentally derived contact data is not available, one can run the ABC model using the powerlaw estimate only. In this case the --HiCdir argument should be excluded from predict.py and the --score_column powerlaw.Score argument should be included in predict.py. In this case the ABC.Score column of the predictions file will be set to NaN. The powerlaw.Score column of the output prediction files will be the relevant Score column to use.

# +------------------------------------+---------------------------------+-------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------+
# |                chr                 |              chr1               |                                           Chromosome of the enhancer and gene                                           |                                                               |
# +------------------------------------+---------------------------------+-------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------+
# | start                              | 1000000                         | Start coordinate of the enhancer element                                                                                |                                                               |
# | end                                | 1000500                         | End coordinate of the enhancer element                                                                                  |                                                               |
# | class                              | intergenic                      | Annotation of enhancer element as {promoter,genic,intergenic}                                                           |                                                               |
# | name                               | intergenic|chr1:1000000-1000500 | unique enhancer region identifier                                                                                       |                                                               |
# | distance                           | 5000                            | Distance in bp between enhancer and TSS of target gene                                                                  |                                                               |
# | isSelfPromoter                     | FALSE                           | Boolean denoting whether element is the promoter of the TargetGene                                                      | Note that EnhancersPredictions.txt does not contain promoters |
# | TargetGene                         | TPTEP1                          | Target Gene                                                                                                             |                                                               |
# | TargetGeneTSS                      | 1005250                         | Transcription Start Site of Target Gene (used to extract Hi-C Contact)                                                  |                                                               |
# | TargetGeneExpression               | 37.96                           | Target Gene Expression                                                                                                  |                                                               |
# | TargetGenePromoterActivityQuantile | 0.95                            | Quantile of Activity at Promoter of target gene                                                                         |                                                               |
# | hic_contact                        | 0.035                           | K-R normalized Hi-C contacts between element and TargetGene TSS                                                         |                                                               |
# | hic_contact_pl_scaled              | 0.035                           | hic_contact scaled by the difference in powerlaw fits between target cell type and reference cell type                  |                                                               |
# | hic_pseudocount                    | 0.001                           | pseudocount added to HiC Contact                                                                                        |                                                               |
# | hic_contact_pl_scaled_adj          | 0.0351                          | Powerlaw scaled KR Normalized HiC plus pseudocount. This is the Contact used in the ABC Score                           |                                                               |
# | activity_base                      | 10.7                            | Geometric mean DHS (or ATAC) and H3K27ac. This is the Activity used in the ABC Score                                    |                                                               |
# | ABC.Score.Numerator                | 0.375                           | The numator of the ABC Score. Activity x Contact without dividing by the scores for other elements near the Target Gene |                                                               |
# | ABC.Score                          | 0.03                            | ABC Score                                                                                                               |                                                               |
# +------------------------------------+---------------------------------+-------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------+
