####################################################################################################################
####################################################################################################################
# Build regression model to predict transcription factor binding on DNA sequences.
# Author: Haiying Kong
# Last Modified: 14 July 2020
####################################################################################################################
####################################################################################################################
setwd('/home/kong/Documents/Personal_Add/Apple/Swiss/Hall_Heim/Project')
options(stringsAsFactors=FALSE)
rm(list=ls())

library(aod)
library(pROC)

####################################################################################################################
# Set paramenters.
# dir.name = '/home/kong/Documents/Personal_Add/Apple/Swiss/Hall_Heim/Project/'
thresh = 8

####################################################################################################################
####################################################################################################################
# Read in classes of the sequences.
Class = read.table('Data/all.data.processed.txt', header=TRUE, sep='\t')[ ,2]
Class = abs(Class-2)

# Read in scores for JASPAR motifs.
scores = read.table('Lock/JASPAR/JASPAR_Scores.txt', header=FALSE, sep='\t')

# Filter noise with threshold.
scores = apply(scores, 2, function(x) {x[x<=thresh]=0; x})

flag = apply(scores[1:1000, ], 2, function(x) length(which(x>0)))
idx = which(flag>300)
scores = scores[ ,idx]
aster = as.data.frame(cbind(Class, scores))

dict = read.table('Lock/JASPAR/JASPAR_Header.txt', header=FALSE, sep='\t')
names(dict) = c('Motif_ID', 'Protein')
dict = dict[idx, ]

names(aster)[-1] = c(dict$Motif_ID)

####################################################################################################################
####################################################################################################################
# Cross validation.
####################################################################################################################
# Shuffle the index of the data.
N = nrow(aster)
idx_all = sample(1:N, N, replace=FALSE)
block_size = N / 10

####################################################################################################################
# Set a vector to save AUC.
AUC = c()

# Logistic regression formula.
# predictors = paste(dict$Motif_ID, collapse=' + ')
# form = formula(paste('Class', predictors, sep=' ~ '))

# cross validation.
for (i_block in 1:10)   {
  # Get training and testing data index.
  idx = ((i_block-1)*block_size+1) : (i_block*block_size)
  idx_train = idx_all[-idx]
  idx_test = idx_all[idx]
  
  # Train logistic regression model.
  model = glm(Class~., family='binomial', data=aster[idx_train, ])

  # Test with test data.
  pred = predict(model, newdata=aster[idx_test,-1], type='response')

  # Get AUC.
  roc.obj = roc(aster$Class[idx_test], pred)
  AUC = c(AUC, auc(roc.obj))
  }

print(mean(AUC))

####################################################################################################################
# Logistic regression on full data and save the results.
model = glm(Class~., family='binomial', data=aster)
apple = summary(model)$coef[-1, c(1,4)]
apple = as.data.frame(cbind(row.names(apple), apple))
names(apple) = c('Motif_ID', 'EffectSize', 'p_value')
apple$EffectSize = as.numeric(as.character(apple$EffectSize))
apple$p_value = as.numeric(as.character(apple$p_value))
row.names(apple) = NULL

apple = merge(dict, apple, by='Motif_ID')
apple = apple[order(apple$p_value), ]

write.table(apple, 'Results/JASPAR_Logistic_pvalues.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')


####################################################################################################################
####################################################################################################################
