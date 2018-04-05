######################################################
#
# DATASET accompagnying publication of
# Salignon, Richard, Fulcrand, Duplus-Bottin and Yvert: 
# Genomics of cellular proliferation in periodic environmental fluctuations
#
######################################################
#
# In this document: N and S corresponds to the two steady conditions. N is the No stress condition. S is the stressful condition. 
# "Stressful" is either Salt 0.2 M (sup_file_Salt.rda) or methionine 0 mM (sup_file_Met.rda) depending on the file/experiment.
#
# The dataset is provided as R objects. You need to install R on your computer, 
# please refer to www.r-project.org to do so, and for documentation on how to use R.
#
# To access the data, open an R console and, after replacing ‘path’ by the path to the directory containing the file on your computer, type in the console:
#
load("path/sup_file_Salt.rda")
#
# this will load data from the salt experiment
#
# or type:
load("path/sup_file_Met.rda")
# to load data from the methionine experiment
#
# You will then have in your R environment the following objects:



#########
pop.inf 
#########
This is a table of annotation describing populations. Data matrices are organized with each column being one population. The rows of this table correspond to the columns in the data matrices, with the same ordering.
Columns in this table are : 
period: N and S correspond to steady environments, 6h-42h correspond to periodic environments with the indicated period.
day: Days 0 to 3 of sampling
replicate: biological replicate of the same condition
wells: coordinate on the plate
primer.barcode: index used for multiplexing
size.factor: normalization factor produced by DESeq2


#########
mut.inf
#########
This is a table of annotation describing mutants. Data matrices are organized with each row being one mutant of the yeast library. Each row of this table corresponds to one row of the data matrices, with the same ordering.
Columns in this table are :
orf.unique: Some open reading frames are present twice in the library (different batches). They are labeled as ORF and ORF.1 in this column
orf: Open reading frame (not necessarily unique)
batch: Batch origin of the mutant, according to Smith et al. Quantitative phenotyping via deep barcode sequencing. Genome Res. 19, 1836–1842 (2009).
gene: Gene name
gene.names: Alternative names of the gene
mutant.barcode: Uptag barcode of the mutant, according to Smith et al. Quantitative phenotyping via deep barcode sequencing. Genome Res. 19, 1836–1842 (2009).
acc.nb: accession number from Euroscarf


#############################
tfit, count.raw, count.norm: 
#############################
These are the data matrices. Mutants are in rows and conditions in columns. 
count.raw: raw count matrix. 
count.norm: count matrix normalized using the vst function from DESEQ2. 
tfit: fitness values. The values described in the paper are the ones for which pop.inf$day equals 3.

##############
tfit.summary: 
##############
Median values of fitness (w) and wdev (wd) across replicates.
Each row corresponds to one mutant.
Columns correspond to:
orf: Open reading frame
gene: Gene name
w.N: Median fitness in N condition
w.S: Median fitness in S condition
w.rat: Ratio of steady controls (w.N/w.S).
w.exp.1: Expected fitness for periods 6h, 12h, 18h, 24h, under the hypothesis of homogeneity. 
w.exp.2: Expected fitness for the period 42h, under the hypothesis of homogeneity.
w.6h to w.42h: Median fitness in periodic environments with the corresponding period.
wd.6h to wd.42h: Deviation from homogeneity. This is the ratio w.observed/w.expected: w.6h/w.exp.1, w.12h/w.exp.1, etc…

######
glms: 
######
Coefficients of the general linear models applied to each mutant.
Each row corresponds to one model (i.e. one mutant, one period).
Columns are:
period: Fluctuation period of the environment
orf: Open reading frame
gene: Gene name
of.N, of.S and of.NS: Fixed offset of the corresponding condition (median of counts at day 0).
est.N, est.S and est.NS: Estimated (fitted) value of the coefficient. In the paper, beta.i.1, beta.i.2, beta.i.3 correspond to est.N, est.S and est.NS, respectively.
se.N, se.S and se.NS: Standard errors of the coefficients
pr.N, pr.S and pr.NS: p-value for the significance of the coefficient (deviation from zero)
qv.N, qv.S and qv.NS: q-value for the significance of the coefficient (deviation from zero)

