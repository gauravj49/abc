# pwd

# Get input files

# 1) Get genes files
# Convert the tsv file to 6 column bed file
# Input file
# +-----+------------+---------+---------+---------+------------------+----------------------+--------------------+-------------+-------------------------+
# | chr |    start   |   end   | feature |  strand |      geneName    |      geneIDVersion   |       geneID       | geneVersion |       geneBiotype       |
# +-----+------------+---------+---------+---------+------------------+----------------------+--------------------+-------------+-------------------------+
# |   1 |    3073253 | 3074322 | gene    | +       |    4933401J01Rik | ENSMUSG00000102693.1 | ENSMUSG00000102693 |           1 |    TEC                  |
# |   1 |    3102016 | 3102125 | gene    | +       |    Gm26206       | ENSMUSG00000064842.1 | ENSMUSG00000064842 |           1 |    snRNA                |
# |   1 |    3205901 | 3671498 | gene    | -       |    Xkr4          | ENSMUSG00000051951.5 | ENSMUSG00000051951 |           5 |    protein_coding       |
# |   1 |    3252757 | 3253236 | gene    | +       |    Gm18956       | ENSMUSG00000102851.1 | ENSMUSG00000102851 |           1 |    processed_pseudogene |
# |   1 |    3365731 | 3368549 | gene    | -       |    Gm37180       | ENSMUSG00000103377.1 | ENSMUSG00000103377 |           1 |    TEC                  |
# +-----+------------+---------+---------+---------+------------------+----------------------+--------------------+-------------+-------------------------+

# output file
# === ============ ========== =============== === === 
# 1      3073253    3074322   4933401J01Rik   0   +  
# 1      3102016    3102125   Gm26206         0   +  
# 1      3205901    3671498   Xkr4            0   -  
# 1      3252757    3253236   Gm18956         0   +  
# 1      3365731    3368549   Gm37180         0   -  
# === ============ ========== =============== === === 
input_file="/media/rad/SSD1/Genome_References/mmu/GRCm38/gencode/Lookups/M24/mmu_GRCm38_gencode_vM24_genes.tsv"
outbed_file="/home/rad/users/gaurav/projects/abc/input/reference/mmu_GRCm38_gencode_vM24_genes.bed"
python -c "import sys; import pandas as pd; import datatable as dt; input_file=sys.argv[1]; outbed_file=sys.argv[2]; genesDT = dt.fread(input_file, sep='\t', header=True, nthreads=16); genesDF = genesDT.to_pandas(); genesDF.insert (4, 'score', 0); newCols=['chr','start','end','geneName', 'score', 'strand']; genesDF = genesDF[newCols]; genesDF.drop_duplicates('geneName', inplace=True); genesDF.to_csv(outbed_file, header=False, index=False, sep='\t', float_format='%.0f'); " ${input_file} ${outbed_file}

# Get ucsc genes list
sed 's/^/chr/' input/reference/mmu_GRCm38_gencode_vM24_genes.bed > input/reference/mmu_GRCm38_gencode_vM24_genes_ucsc.bed