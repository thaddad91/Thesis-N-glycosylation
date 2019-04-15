#!/usr/bin/env python

# Author: Thierry Haddad
# Description: Perform PCA on kmer set

import pandas as pd
from sklearn.decomposition import PCA

df = pd.read_csv("kmer_list.txt",sep='\t',header=0, index_col=0)
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(df)
principalDf = pd.DataFrame(data = principalComponents,
                           columns = ['principal component 1', 'principal component 2'])
finalDf = pd.concat([principalDf, df[['glycosite']]], axis = 1)

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1)
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
targets = ["pos", "neg"]
colors = ['g', 'r']
for target, color in zip(targets,colors):
    indicesToKeep = finalDf['glycosite'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)
ax.legend(targets)
ax.grid()