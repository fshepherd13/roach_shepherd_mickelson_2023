import pandas as pd
import sys

rna_ratios = sys.argv[1]
splicing = sys.argv[2]
expression = sys.argv[3]
sample_id = sys.argv[4]
out_file = sys.argv[5]


#Read in ratios of mRNA:total mRNA for PB2, PB1, PA, HA, NP, NA, M AND S
ratios = pd.read_table(rna_ratios, index_col=0, keep_default_na=False, na_values="N/A")

#Read in splicing information for M1/M2 and NS1/NS2
splice_counts = pd.read_table(splicing)

#Read in cuffdiff expression levels, extracting only the fpkm values for the sample ID being analyzed
expression_levels = pd.read_table(expression, index_col=4, keep_default_na=False)[sample_id+"_FPKM"]

#Calculate total fpkm using all genes in cuffdiff results by summing all fpkm values
total_FPKM = expression_levels.sum()

#Get rid of any non-flu genes, as quantification from here only requires IAV segment expression
expression_levels = expression_levels["HA":"vPB2"].to_frame(name='FPKM')


############First, perform calculations necessary for M transcript quantification############
#Define M2 <-depth of spliced M2 mRNA
#Takes difference in read depth at right and left splice junctions (to get the depth of spliced M2 transcript), and averages them.
M2 = (1/2)*((splice_counts['Depth_total_left'].values[0] - splice_counts['Depth_unspliced_left'].values[0]) + (splice_counts['Depth_total_right'].values[0] - splice_counts['Depth_unspliced_right'].values[0]))

#Define mM <- read counts of total mRNA of M gene. Takes average read depth at splice junctions where spliced+unspliced transcripts are present, multiplied by ratio of mRNA to total positive RNA for M
mrna_total_rna_ratio_m = float(ratios.loc['M', 'mRNA:total_pos_RNA'])

mM = (1/2)*(splice_counts['Depth_total_left'].values[0] + splice_counts['Depth_total_right'].values[0]) * mrna_total_rna_ratio_m

#define "f_M" <- fraction of spliced M2 read counts to total M mRNA read counts:
if mM != 0:
    f_M = M2 / mM
else:
    f_M = 0

#Calculate l_M <- factor to adjust for length difference between spliced and unspliced transcripts
if splice_counts['Depth_unspliced_left'].values[0] !=0 and splice_counts['Depth_total_left'].values[0] != 0:
    l_M = 1027/((splice_counts['Depth_unspliced_left'].values[0]*1027 / splice_counts['Depth_total_left'].values[0]) + (M2 * 338 / splice_counts['Depth_total_left'].values[0]))
else:
    l_M = 0
    
############Perform same calculations for NS############
#Define cNEP <- read counts of spliced NEP
#Takes difference in read depth at right and left splice junctions (to get the depth of spliced transcript), and averages them.
cNEP = (1/2)*((splice_counts['Depth_total_left'].values[1] - splice_counts['Depth_unspliced_left'].values[1]) + (splice_counts['Depth_total_right'].values[1] - splice_counts['Depth_unspliced_right'].values[1]))

#Define c(mNS) <- read counts of total mRNA of NS gene. Takes average read depth at splice junctions where spliced+unspliced transcripts are present, multiplied by ratio of mrna to total positive rna for NS
mrna_total_rna_ratio_ns = float(ratios.loc['NS', 'mRNA:total_pos_RNA'])

cmNS = (1/2)*(splice_counts['Depth_total_left'].values[1] + splice_counts['Depth_total_right'].values[1]) * mrna_total_rna_ratio_ns

#define "f_NS" <- fraction of spliced NEP read counts to total NS mRNA read counts:
if cmNS != 0:
    f_NS = cNEP / cmNS
else:
    f_NS = 0

#Calculate l_NS <- factor to adjust for length difference between spliced and unspliced transcripts
if splice_counts['Depth_unspliced_left'].values[0] != 0 and splice_counts['Depth_total_left'].values[0] != 0:
    l_NS = 890/((splice_counts['Depth_unspliced_left'].values[0]*890 / splice_counts['Depth_total_left'].values[0]) + (M2 * 418 / splice_counts['Depth_total_left'].values[0]))
else:
    l_NS = 0

#Pull FPKM values for positive sense transcripts
HA_fpkm = expression_levels.loc['HA', 'FPKM']
NA_fpkm = expression_levels.loc['NA', 'FPKM']
M_fpkm = expression_levels.loc['M', 'FPKM']
NP_fpkm = expression_levels.loc['NP', 'FPKM']
PA_fpkm = expression_levels.loc['PA', 'FPKM']
NS_fpkm = expression_levels.loc['NS', 'FPKM']
PB1_fpkm = expression_levels.loc["PB1", 'FPKM']
PB2_fpkm = expression_levels.loc["PB2", 'FPKM']

#Pull FPKM values for negative sense transcripts
vHA_fpkm = expression_levels.loc['vHA', 'FPKM']
vNA_fpkm = expression_levels.loc['vNA', 'FPKM']
vM_fpkm = expression_levels.loc['vM', 'FPKM']
vNP_fpkm = expression_levels.loc['vNP', 'FPKM']
vPA_fpkm = expression_levels.loc['vPA', 'FPKM']
vNS_fpkm = expression_levels.loc['vNS', 'FPKM']
vPB1_fpkm = expression_levels.loc['vPB1', 'FPKM']
vPB2_fpkm = expression_levels.loc['vPB2', 'FPKM']

#Create dataframe to hold all results
#First, insert the FPKM values from the Cufflinks results and mrna:(mrna+crna) ratios
results = pd.DataFrame(
    [[HA_fpkm, vHA_fpkm, ratios.loc['HA', 'mRNA:total_pos_RNA']],
    [NA_fpkm, vNA_fpkm, ratios.loc['NA', 'mRNA:total_pos_RNA']],
    [M_fpkm, vM_fpkm, ratios.loc['M', 'mRNA:total_pos_RNA']],
    [NP_fpkm, vNP_fpkm, ratios.loc['NP', 'mRNA:total_pos_RNA']],
    [PA_fpkm, vPA_fpkm, ratios.loc['PA', 'mRNA:total_pos_RNA']],
    [NS_fpkm, vNS_fpkm, ratios.loc['NS', 'mRNA:total_pos_RNA']],
    [PB1_fpkm, vPB1_fpkm, ratios.loc['PB1', 'mRNA:total_pos_RNA']],
    [PB2_fpkm, vPB2_fpkm, ratios.loc['PB2', 'mRNA:total_pos_RNA']],
    [None, None, None], 
    [None, None, None],
    [None, None, None],
    [None, None, None]],
    columns=['FPKM_positive_sense', 'FPKM_negative_sense', 'mrna:(mrna,crna)_ratio'],
    index = ['HA', 'NA', 'M', 'NP', 'PA', 'NS', 'PB1', 'PB2', 'NS1', 'NEP', 'M1', 'M2'])

#Transform FPKM to TPM for positive and negative sense RNA transcripts
results['TPM_positive_sense'] = (results['FPKM_positive_sense'] / total_FPKM)*(1000000)
results['TPM_negative_sense'] = (results['FPKM_negative_sense'] / total_FPKM)*(1000000)

#Create column for vRNA TPM counts (same as the direct cufflinks output)
results['vrna_TPM'] = results['TPM_negative_sense']

#Calculate crna TPM
results['crna_TPM'] = results['TPM_positive_sense'] * (1-results['mrna:(mrna,crna)_ratio'])

#Define function for calculating mRNA read counts by multiplying total positive RNA by mrna:total positive rna ratio
def mrna_calc(gene):
    mrna = results.loc[str(gene), 'TPM_positive_sense'] * results.loc[str(gene), 'mrna:(mrna,crna)_ratio']
    return mrna

#initialize empty dictionary to hold mrna TPM calculations. Run the mrna_calc function on all genes except M and NS
d = dict()
for x in ['HA', 'NA', 'NP', 'PA', 'PB1', 'PB2']:
    d[x] = mrna_calc(x)

#M and NS don't exist as mrna, so add these entries to dictionary with "NA"
d['M'] = None
d['NS'] = None


#Calculate TPM of mrna for unspliced M1 and NS1: TPM(all positive sense) * ratio of mRNA:positive sense RNA * (1 - ratio of spliced reads:mRNA) * (length correction factor "L")
d['NS1'] = results.loc['NS', 'TPM_positive_sense'] * (1 - f_NS) * results.loc['NS', 'mrna:(mrna,crna)_ratio'] * l_NS

d['M1'] = results.loc['M', 'TPM_positive_sense'] * (1-f_M) * results.loc['M', 'mrna:(mrna,crna)_ratio'] * l_M

#Calculate mrna for spliced M2 and NEP:
d['NEP'] = results.loc['NS', 'TPM_positive_sense'] * f_NS * results.loc['NS', 'mrna:(mrna,crna)_ratio'] * l_NS

d['M2'] = results.loc['M', 'TPM_positive_sense'] * (f_M) * results.loc['M', 'mrna:(mrna,crna)_ratio'] * l_M

#Map the full dictionary back to the results dataframe
results['mrna_TPM'] = results.index.map(d)

results.to_csv(out_file)