####################################################################################################################
####################################################################################################################
# Get scores of JASPAR motifs for all sequences.
# Author: Haiying Kong
# Last Modified: 13 July 2020
####################################################################################################################
####################################################################################################################
from Bio import motifs
import pymysql
pymysql.install_as_MySQLdb()
from Bio.motifs.jaspar.db import JASPAR5
from Bio.Seq import Seq
import pickle
import numpy as np

####################################################################################################################
# Set parameters.
threshold = 8.0

####################################################################################################################
# Read in sequence data.
pkl_file = open('/home/kong/Documents/Personal_Add/Apple/Swiss/Hall_Heim/Project/Data/data.pkl', 'rb')
data = pickle.load(pkl_file)
pkl_file.close()

seqs = data['Seq_Mat']

####################################################################################################################
# Get all JASPAR motifs.
with open("/home/kong/Database/JASPAR/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt") as f:
  pfms = motifs.parse(f, 'jaspar')

names = []
pfm_ids = []
apple = []
k = 0

for pfm in pfms:
  names.append(pfm.name)
  pfm_ids.append(pfm.matrix_id)
  motif = pfm.counts
  motif = np.array(list(pfm.counts.values()))
  motif = motif/motif.sum(axis=0)
  motif = motif.transpose()
  scores = []
  for seq in seqs:
    aster = []
    for i in range(seq.shape[0]-motif.shape[0]+1):
      aster.append(sum(sum(seq[i:(i+motif.shape[0]), :] * motif)))
    score = max(aster)
    scores.append(score)
  apple.append(scores)
  k = k+1
  print(k)

apple = np.array(apple)
apple = np.transpose(apple)

####################################################################################################################
# Save header as a file.
header = np.transpose(np.array([pfm_ids, names]))
np.savetxt('/home/kong/Documents/Personal_Add/Apple/Swiss/Hall_Heim/Project/Lock/JASPAR/JASPAR_Header.txt', header, fmt='%s', delimiter='\t')

# Save the scores along with protein names and motif IDs.
np.savetxt('/home/kong/Documents/Personal_Add/Apple/Swiss/Hall_Heim/Project/Lock/JASPAR/JASPAR_Scores.txt', apple, fmt='%f', delimiter='\t')

####################################################################################################################
####################################################################################################################

