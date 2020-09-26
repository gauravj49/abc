# pwd
cd /home/rad/users/gaurav/projects/abc

# Define the bams and the peaks files
h3k27acBAMs=(/media/rad/HDD1/nfchip/christine/pdacBatch1/gjchip/mapping/pdacBatch1_20200424121526_A0000498_127abcam-5320-PPT-1_mmu_chipseq_se_R1_rmdup.bam /media/rad/HDD1/nfchip/christine/pdacBatch1/gjchip/mapping/pdacBatch1_20200424121526_A0000498_128abcam-5320-LungMet-1_mmu_chipseq_se_R1_rmdup.bam /media/rad/HDD1/nfchip/christine/pdacBatch1/gjchip/mapping/pdacBatch1_20200424121526_A0000498_129abcam-5320-LivMet-1_mmu_chipseq_se_R1_rmdup.bam /media/rad/HDD1/nfchip/christine/pdacBatch1/gjchip/mapping/pdacBatch1_20200424121526_A0000498_130abcam-5320-LivMet-3_mmu_chipseq_se_R1_rmdup.bam /media/rad/HDD1/nfchip/christine/pdacBatch1/gjchip/mapping/pdacBatch1_20200424121526_A0000498_123abcam-53646-PPT-1_mmu_chipseq_se_R1_rmdup.bam /media/rad/HDD1/nfchip/christine/pdacBatch1/gjchip/mapping/pdacBatch1_20200424121526_A0000498_124abcam-53646-LivMet-1_mmu_chipseq_se_R1_rmdup.bam /media/rad/HDD1/nfchip/christine/pdacBatch1/gjchip/mapping/pdacBatch1_20200424121526_A0000498_125abcam-53646-LivMet-2_mmu_chipseq_se_R1_rmdup.bam /media/rad/HDD1/nfchip/christine/pdacBatch1/gjchip/mapping/pdacBatch1_20200424121526_A0000498_126abcam-53646-LivMet-3_mmu_chipseq_se_R1_rmdup.bam)
sampleNames=(5320_PPT-1 5320_LungMet-1 5320_LivMet-1 5320_LivMet-3 53646_PPT-1 53646_LivMet-1 53646_LivMet-2 53646_LivMet-3)

# Add chr to bam and broad peaks file
broadPeaksUcscDir="/media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/ucsc"; mkdir -p ${broadPeaksUcscDir}
atacBamsUcscDir="/media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/ucsc"; mkdir -p ${atacBamsUcscDir}
atacpeaks=(/media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/5320_PPT-1_005_R1.mLb.clN_peaks.broadPeak /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/5320_LungMet-1_005_R1.mLb.clN_peaks.broadPeak /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/5320_LivMet-1_005_R1.mLb.clN_peaks.broadPeak /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/5320_LivMet-3_005_R1.mLb.clN_peaks.broadPeak /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/53646_PPT-1_005_R1.mLb.clN_peaks.broadPeak /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/53646_LivMet-1_005_R1.mLb.clN_peaks.broadPeak /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/53646_LivMet-2_005_R1.mLb.clN_peaks.broadPeak /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/53646_LivMet-3_005_R1.mLb.clN_peaks.broadPeak)
atacbams=(/media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/5320_PPT-1_005_R1.mLb.clN.sorted.bam /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/5320_LungMet-1_005_R1.mLb.clN.sorted.bam /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/5320_LivMet-1_005_R1.mLb.clN.sorted.bam /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/5320_LivMet-3_005_R1.mLb.clN.sorted.bam /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/53646_PPT-1_005_R1.mLb.clN.sorted.bam /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/53646_LivMet-1_005_R1.mLb.clN.sorted.bam /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/53646_LivMet-2_005_R1.mLb.clN.sorted.bam /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/53646_LivMet-3_005_R1.mLb.clN.sorted.bam)
for ((i=0;i<${#atacpeaks[@]};++i)); do
  # Add chr to atac broadpeak and replacing chrMT to chrM
  echo "- Adding hr to atac broadpeak and replacing chrMT to chrM: ${sampleNames[i]}"
  sed 's/^/chr/' ${atacpeaks[i]}| sed 's/^chrMT/chrM/' > ${broadPeaksUcscDir}/$(basename ${atacpeaks[i]} .broadPeak)_ucsc.broadPeak
  
  # Add chr to atac bams
  echo "- Adding chr to atac bams and creating the index: ${sampleNames[i]}"
  # samtools view -h ${atacbams[i]}| sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' | samtools view -bS > ${atacBamsUcscDir}/$(basename ${atacbams[i]} .bam)_ucsc.bam
  samtools view -h ${atacbams[i]}| sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader -  ${atacbams[i]}> ${atacBamsUcscDir}/$(basename ${atacbams[i]} .bam)_ucsc.bam
  samtools index ${atacBamsUcscDir}/$(basename ${atacbams[i]} .bam)_ucsc.bam
done

# Get the expression data
# Initial data 272KC_Mouse_DGE_Matrix.txt and 272KC_SampleInfo.xlsx
# From this file, generated pdac1_expression_data_Raw_TPM.xlsx containing raw and tpm normalized data
# Extract two column genes <> expression in an output file for each cellLine from the pdac1_expression_data_TPM.txt
# Get the counts data:
# head /media/rad/HDD1/abc/christine/pdacBatch1/expression/pdac1_expression_data_TPM.txt | cut -f 1-6 | column -t
# GENE           53646-PPT-1        53646-LivMet-1     53646-LivMet-2     53646-LivMet-3     6075-PPT-1
# 0610005C13Rik  0.161751680114701  0                  0                  0.206100704514038  0
# 0610009B22Rik  40.4379200286753   44.6233516246271   43.2285692155048   52.1434782420517   44.6481431289655
# 0610009E02Rik  0.323503360229403  0.642062613303987  0.295075557785015  0.206100704514038  0.983439275968403
# ...
infile="/media/rad/HDD1/abc/christine/pdacBatch1/expression/pdac1_expression_data_TPM.txt"; dos2unix ${infile}
countsDIR="/media/rad/HDD1/abc/christine/pdacBatch1/expression/counts"; mkdir -p ${countsDIR}
ncols=`expr $(awk '{print NF; exit}' ${infile})`;
for ((i=2;i<=ncols;i++)); do
 bname=`head -1 ${infile}|cut -f${i}`
 ofname="${countsDIR}/${bname}_Counts.txt"
 echo "processing $bname..."
 cut -f1,${i} ${infile} --output-delimiter=$'\t' | grep -v ${bname} > ${ofname}
done

# ==> 5320-LivMet-1_Counts.txt <==
# 0610005C13Rik   0
# 0610009B22Rik   43.4488503042964

# ==> 5320-LivMet-3_Counts.txt <==
# 0610005C13Rik   0.187353173659718
# 0610009B22Rik   46.2762338939502

# ==> 5320-LungMet-1_Counts.txt <==
# 0610005C13Rik   0
# 0610009B22Rik   44.3845881044952

# Define the bams and the peaks files
sampleNames=(5320_PPT-1 5320_LungMet-1 5320_LivMet-1 5320_LivMet-3 53646_PPT-1 53646_LivMet-1 53646_LivMet-2 53646_LivMet-3)
peaksdir="/media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/ucsc"
atacbamdir="/media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/ucsc"
h3k27acbamdir="/media/rad/HDD1/nfchip/christine/pdacBatch1/gjchip/mapping"
countsdir="/media/rad/HDD1/abc/christine/pdacBatch1/expression/counts"
species="mouse"

for ((i=0;i<${#sampleNames[@]};++i)); do
  echo "########################################"
  echo "- Running for sample: ${sampleNames[i]}"
  atacpeak=$(find ${peaksdir} -name "*${sampleNames[i]}*_ucsc.broadPeak")
  atacbam=$(find ${atacbamdir} -name "*${sampleNames[i]}*.mLb.clN.sorted_ucsc.bam")
  h3k27acbam=$(find ${h3k27acbamdir} -name "*${sampleNames[i]}*_rmdup.bam" | grep -v input)
  countsfile=$(find ${countsdir} -name "*${sampleNames[i]}*")
  outputDir="/media/rad/HDD1/abc/christine/pdacBatch1/${sampleNames[i]}"; mkdir -p ${outputDir}
  echo "bash scripts/run_abc_wrapper.sh ${atacpeak} ${atacbam} ${h3k27acbam} ${countsfile} ${sampleNames[i]} ${outputDir} ${species}"
  bash scripts/run_abc_wrapper.sh ${atacpeak} ${atacbam} ${h3k27acbam} ${countsfile} ${sampleNames[i]} ${outputDir} ${species}
done



# atacpeaks=(/media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/ucsc/5320_PPT-1_005_R1.mLb.clN_peaks_ucsc.broadPeak /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/ucsc/5320_LungMet-1_005_R1.mLb.clN_peaks_ucsc.broadPeak /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/ucsc/5320_LivMet-1_005_R1.mLb.clN_peaks_ucsc.broadPeak /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/ucsc/5320_LivMet-3_005_R1.mLb.clN_peaks_ucsc.broadPeak /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/ucsc/53646_PPT-1_005_R1.mLb.clN_peaks_ucsc.broadPeak /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/ucsc/53646_LivMet-1_005_R1.mLb.clN_peaks_ucsc.broadPeak /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/ucsc/53646_LivMet-2_005_R1.mLb.clN_peaks_ucsc.broadPeak /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/ucsc/53646_LivMet-3_005_R1.mLb.clN_peaks_ucsc.broadPeak)
# atacbams=(/media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/ucsc/5320_PPT-1_005_R1.mLb.clN.sorted_ucsc.bam /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/ucsc/5320_LungMet-1_005_R1.mLb.clN.sorted_ucsc.bam /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/ucsc/5320_LivMet-1_005_R1.mLb.clN.sorted_ucsc.bam /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/ucsc/5320_LivMet-3_005_R1.mLb.clN.sorted_ucsc.bam /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/ucsc/53646_PPT-1_005_R1.mLb.clN.sorted_ucsc.bam /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/ucsc/53646_LivMet-1_005_R1.mLb.clN.sorted_ucsc.bam /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/ucsc/53646_LivMet-2_005_R1.mLb.clN.sorted_ucsc.bam /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/ucsc/53646_LivMet-3_005_R1.mLb.clN.sorted_ucsc.bam)
# countsfile=(/media/rad/HDD1/abc/christine/pdacBatch1/expression/counts/5320_LivMet-1_Counts.txt /media/rad/HDD1/abc/christine/pdacBatch1/expression/counts/5320_PPT-1_Counts.txt /media/rad/HDD1/abc/christine/pdacBatch1/expression/counts/53646_LivMet-3_Counts.txt /media/rad/HDD1/abc/christine/pdacBatch1/expression/counts/5320_LivMet-3_Counts.txt /media/rad/HDD1/abc/christine/pdacBatch1/expression/counts/53646_LivMet-1_Counts.txt /media/rad/HDD1/abc/christine/pdacBatch1/expression/counts/53646_PPT-1_Counts.txt /media/rad/HDD1/abc/christine/pdacBatch1/expression/counts/5320_LungMet-1_Counts.txt /media/rad/HDD1/abc/christine/pdacBatch1/expression/counts/53646_LivMet-2_Counts.txt)
# # Run ABC pipeline on the above samples
# for ((i=0;i<${#atacpeaks[@]};++i)); do
#   echo "########################################"
#   echo "- Running for sample: ${sampleNames[i]}"
#   outputDir="/media/rad/HDD1/abc/christine/pdacBatch1/${sampleNames[i]}"; mkdir -p ${outputDir}
#   echo "bash scripts/run_abc_wrapper.sh ${atacpeaks[i]} ${atacbams[i]} ${h3k27acBAMs[i]} ${sampleNames[i]} ${outputDir}"
#   # bash scripts/run_abc_wrapper.sh ${atacpeaks[i]} ${atacbams[i]} ${h3k27acBAMs[i]} ${sampleNames[i]} ${outputDir}
#   echo ""
# done

# # bash scripts/run_abc_wrapper.sh /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/macs/broadPeak/ucsc/5320_PPT-1_005_R1.mLb.clN_peaks_ucsc.broadPeak /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/results/bwa/mergedLibrary/ucsc/5320_PPT-1_005_R1.mLb.clN.sorted_ucsc.bam /media/rad/HDD1/nfchip/christine/pdacBatch1/gjchip/mapping/pdacBatch1_20200424121526_A0000498_127abcam-5320-PPT-1_mmu_chipseq_se_R1_rmdup.bam 5320_PPT-1 /media/rad/HDD1/abc/christine/pdacBatch1/5320_PPT-1
