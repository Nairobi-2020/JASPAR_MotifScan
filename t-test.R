####################################################################################################################
####################################################################################################################
# Perform t-test to identify significant motifs.
# Author: Haiying Kong
# Last Modified: 14 July 2020
####################################################################################################################
####################################################################################################################
setwd('/home/kong/Documents/')
options(stringsAsFactors=FALSE)
rm(list=ls())

####################################################################################################################
# Set paramenters.
# dir.name = '/home/kong/Documents/'
thresh = 8

####################################################################################################################
####################################################################################################################
# Read in classes of the sequences.
Class = read.table('Data/all.data.processed.txt', header=TRUE, sep='\t')[ ,2]
Class = abs(Class-2)
idx.chip = which(Class==1)
idx.control = which(Class==0)

# Read in scores for JASPAR motifs.
scores = read.table('Lock/JASPAR/JASPAR_Scores.txt', header=FALSE, sep='\t')

# Filter scores with threshold.
scores = apply(scores, 2, function(x) {x[x<=thresh]=0; x})
flag = apply(scores[1:1000, ], 2, function(x) length(which(x>0)))
idx = which(flag>300)
scores = scores[ ,idx]

# Get dictionary.
dict = read.table('Lock/JASPAR/JASPAR_Header.txt', header=FALSE, sep='\t')
names(dict) = c('Motif_ID', 'Protein')
dict = dict[idx, ]

# Perform t-test.
p.values = apply(scores, 2, function(x) t.test(x[idx.chip], x[idx.control])$p.value)
FDR = p.adjust(p.values, 'fdr', n=length(p.values))

# Get final summary table for FDR and save.
apple = cbind(dict, FDR)
apple = apple[order(apple$FDR), ]
write.table(apple, 'Results/JASPAR_FDR.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')


####################################################################################################################
####################################################################################################################
