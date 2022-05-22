'''
find BP in all introns in chromosome 
genome:  input data/hg19.fa
introns: input data/ann_gencode_v19.bed

chrom - parameter ('chr1', ..., 'chrY')

output: predicted BPs for all chromosomes/chrom.csv
'''

from Bio import SeqIO
from tqdm import tqdm
from collections import Counter
import seaborn as sns
import pandas as pd
import numpy as np
import itertools
import matplotlib.pyplot as plt

from sklearn.utils import shuffle
from sklearn.model_selection import train_test_split
from sklearn.utils.random import sample_without_replacement
from sklearn.model_selection import cross_val_score

from sklearn.ensemble import RandomForestClassifier
import sys

#Functions
def count_kmers(read, k):
    '''
    count k-mers distribution in read
    output: dictionary {k-mer:number of k-mer in this read}
    '''
    read = str(read).upper()
    counts = {}
    num_kmers = len(read) - k + 1
    for i in range(num_kmers):
        kmer = read[i:i+k]
        if 'N' in kmer: 
            continue
        if kmer not in counts:
            counts[kmer] = 0
        counts[kmer] += 1
    return counts

def kmers_pos(chrom, pos, n, k):
    '''
    count k-mers dist in (pos-n, pos+n)
    output: sorted by k-mer names disctionary {k-mer: number of k-mer in (pos-n, pos+n)}
    '''
    genom_chr = hg19[chrom]

    kmers = {}
    bp_round = genom_chr[pos-n:pos+n]
    cnt = count_kmers(bp_round.seq, k)

    kmers = Counter(cnt)
    kmers = dict(sorted(kmers.items()))
    
    return kmers

# Find BP in chromosome
def find_BP(chrom):

    X = np.load('Xy for training model (n=70, k=5)/X_70_5'+chrom+'.npy')
    X = pd.DataFrame(X, columns = ['pos', *[''.join(i) for i in itertools.product('ATGC', repeat = 5)], 'dist3'])
    y = np.load('Xy for training model (n=70, k=5)/y_70_5'+chrom+'.npy')

    model = RandomForestClassifier()
    model.fit(X, y)
    
    introns_chrom = introns[introns['chr'] == chrom]
    BP = pd.DataFrame(columns = [*introns_chrom.columns, 'BP'])
    
    print('Find BP for all introns', )
    for i in tqdm(range(len(introns_chrom))):
        near_to_3 = []
        dist_to_3 = {}
        start, stop, strand = introns_chrom.iloc[i]['start'], introns_chrom.iloc[i]['stop'], introns_chrom.iloc[i]['strand']
        for j in range(start, stop+1):
            if strand == '+': #start -> stop is 5' -> 3' 
                dist_to_3[j] = stop - j
                if stop - j <= 50:
                    near_to_3.append(j)
            else:             #start -> stop is 3' -> 5' 
                dist_to_3[j] = j - start
                if j - start <= 50:
                    near_to_3.append(j)

        X_bp = pd.DataFrame(columns = [''.join(i) for i in itertools.product('ATGC', repeat = k)])
        for pos in near_to_3:
            X_bp = X_bp.append(kmers_pos(chrom, pos, n, k), ignore_index=True)
        X_bp.insert(loc=0, column='pos', value = near_to_3)
        X_bp = X_bp.fillna(0)
        X_bp['dist3'] = [dist_to_3[x] for x in X_bp['pos']]
        p_pred = model.predict_proba(X_bp)
        bp_cur = X_bp['pos'].iloc[np.where(p_pred[:, 1] == max(p_pred[:, 1]))[0][0]]
        row_cur = [*list(introns_chrom.iloc[i]), bp_cur]
        BP = BP.append({BP.columns[i]: row_cur[i] for i in range(6)}, ignore_index=True).drop_duplicates()
    BP.to_csv('predicted BPs for all chromosomes/' + chrom + '.csv', index=False)
    return

#parameters
n, k = 70, 5
chrom = sys.argv[1]

print('Read genome from input data/hg19.fa')
#reference genome
hg19 = SeqIO.to_dict(SeqIO.parse("input data/hg19.fa", "fasta"))

print('Read introns positions from input data/ann_gencode_v19.bed')
#introns positions
introns = pd.read_csv('input data/ann_gencode_v19.bed', sep = '\t', header = None, names = ['chr', 'start', 'stop', 'ID', 'score', 'strand'])
introns = introns[introns['chr'].isin(hg19.keys())]
introns.drop('ID', inplace=True, axis=1)
introns = introns.drop_duplicates()

# Functions for model creation
find_BP(chrom)
print('Results are here: predicted BPs for all chromosomes/'+chrom+'.csv')