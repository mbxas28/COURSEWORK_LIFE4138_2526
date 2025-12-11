#!/usr/bin/env python
# coding: utf-8

# In[7]:


#importing required packages to run the analysis

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# In[8]:


#reading in the file
A_vs_E = pd.read_csv(r'D:/UoN/modules/4138 - Coding/coursework/assignment work/Python script - pycharm/A_vs_E.deseq2.results.tsv', sep='\,')
A_vs_E.head()


# In[9]:


#reading in second file
A_vs_F = pd.read_csv('D:/UoN/modules/4138 - Coding/coursework/assignment work/Python script - pycharm/A_vs_F.deseq2.results.tsv', sep='\,')
A_vs_F.head()


# In[20]:


#setting a threshold value to use and compare to find significant genes later
pval_threshold = 0.05
log2FC_threshold = 1


# In[21]:


#making a variable with the significant values of genes depending on two separate column values for A vs E and A vs F datasets
significant_genes_E = A_vs_E[(A_vs_E['pvalue'] < pval_threshold) & (A_vs_E['log2FoldChange'].abs() >= log2FC_threshold)]

significant_genes_F = A_vs_F[(A_vs_F['pvalue'] < pval_threshold) & (A_vs_F['log2FoldChange'].abs() >= log2FC_threshold)]


# In[12]:


#counting the number of up and downregulated genes in each dataset
#for A vs E

upregulated_genes_E = significant_genes_E[significant_genes_E['log2FoldChange'] > 0].shape[0]
print(f"Number of upregulated genes are: {upregulated_genes_E}")

downregulated_genes_E = significant_genes_E[significant_genes_E['log2FoldChange'] < 0].shape[0]
print(f"Number of downregulated genes are: {downregulated_genes_E}")


# In[13]:


#for A vs F

upregulated_genes_F = significant_genes_F[significant_genes_F['log2FoldChange'] > 0].shape[0]
print(f"Number of upregulated genes are: {upregulated_genes_F}")

downregulated_genes_F = significant_genes_F[significant_genes_F['log2FoldChange'] < 0].shape[0]
print(f"Number of downregulated genes are: {downregulated_genes_F}")


# In[14]:


#summary statistics for pvalue and log fold changes
#for A_vs_E

A_vs_E[["pvalue", "log2FoldChange"]].describe()


# In[15]:


#for A_vs_F

A_vs_F[["pvalue", "log2FoldChange"]].describe()


# In[87]:


#volcano plots for both comparisons
#dropping all values that have NA before plotting 
A_vs_E['NegLog10PV'] = -np.log10(A_vs_E['pvalue'].dropna()) 

#for A vs E
plt.figure(figsize=(4,6), dpi = 100)
plt.scatter(A_vs_E["log2FoldChange"], A_vs_E["NegLog10PV"], color = 'skyblue', alpha = 0.7)
plt.xlabel("Log2 Fold Change", fontsize = 8)
plt.ylabel("-Log10 P-value", fontsize = 8)
plt.title("Volcano plot of significance and magnitude of changes in gene expression", fontsize = 10)
plt.axhline(y = 1.3, linestyle = '--', color = 'green') #pvalue threshold converted to -log10
plt.axvline(x = -1, linestyle = '--', color = 'red') #log2FoldChange threshold (+ve and -ve)
plt.axvline(x = 1, linestyle = '--', color = 'red')
plt.show


# In[78]:


#for A vs F

A_vs_F['NegLog10PV'] = -np.log10(A_vs_F['pvalue'].dropna())

plt.figure(figsize=(4,6), dpi = 100)
plt.scatter(A_vs_F["log2FoldChange"], A_vs_F["NegLog10PV"], color = 'darkblue' , alpha = 0.6)
plt.xlabel("Log2 Fold Change", fontsize = 8)
plt.ylabel("-Log10 P-value", fontsize = 8)
plt.title("Volcano plot of significance and magnitude of changes in gene expression", fontsize = 10)
plt.axhline(y = 1.3, linestyle = '--', color = 'green') #pvalue threshold converted to -log10
plt.axvline(x = -1, linestyle = '--', color = 'red') #log2FoldChange threshold (+ve and -ve)
plt.axvline(x = 1, linestyle = '--', color = 'red')
plt.show


# In[79]:


#MA plot
#for A vs E
#converting baseMean to log10baseMean to use in MA plot
A_vs_E["log10baseMean"] = np.log10(A_vs_E["baseMean"])


#declaring significant, up and downregulated genes
sig = (A_vs_E["pvalue"] < pval_threshold) & (A_vs_E["log2FoldChange"].abs() >= log2FC_threshold)
up = sig & (A_vs_E["log2FoldChange"] > 0)
down = sig & (A_vs_E["log2FoldChange"] < 0)
plt.figure(figsize=(6,4))

#non-significant genes selected from the columns specified 
plt.scatter(
    A_vs_E.loc[~sig, "log10baseMean"],
    A_vs_E.loc[~sig, "log2FoldChange"],
    color = 'gold',
    s=10, alpha=0.7)

#upregulated 
plt.scatter(
    A_vs_E.loc[up, "log10baseMean"],
    A_vs_E.loc[up, "log2FoldChange"],
    color = 'maroon',
    s=10, alpha=0.7)

#downregulated 
plt.scatter(
    A_vs_E.loc[down, "log10baseMean"],
    A_vs_E.loc[down, "log2FoldChange"],
    color = 'lightcoral',
    s=10, alpha=0.7)

#plotting the MA plot
plt.axhline(0, color="black", linewidth=1)
plt.xlabel("log10(baseMean)")
plt.ylabel("log2FoldChange")
plt.title("MA Plot: A vs E")
plt.show()


# In[80]:


#MA plot
#for A vs F
#converting baseMean to log10baseMean to use in MA plot
A_vs_F["log10baseMean"] = np.log10(A_vs_F["baseMean"])


#declaring significant, up and downregulated genes
sigF = (A_vs_F["pvalue"] < pval_threshold) & (A_vs_F["log2FoldChange"].abs() >= log2FC_threshold)
upF = sig & (A_vs_F["log2FoldChange"] > 0)
downF = sig & (A_vs_F["log2FoldChange"] < 0)
plt.figure(figsize=(6,4))

#non-significant genes selected from the columns specified 
plt.scatter(
    A_vs_F.loc[~sigF, "log10baseMean"],
    A_vs_F.loc[~sigF, "log2FoldChange"],
    color = 'turquoise',
    s=10, alpha=0.7)

#upregulated 
plt.scatter(
    A_vs_F.loc[upF, "log10baseMean"],
    A_vs_F.loc[upF, "log2FoldChange"],
    color = 'darkslategrey',
    s=10, alpha=0.7)

#downregulated 
plt.scatter(
    A_vs_F.loc[downF, "log10baseMean"],
    A_vs_F.loc[downF, "log2FoldChange"],
    color = 'deepskyblue',
    s=10, alpha=0.7)

#plotting the MA plot
plt.axhline(0, color="black", linewidth=1)
plt.xlabel("log10(baseMean)")
plt.ylabel("log2FoldChange")
plt.title("MA Plot: A vs F")
plt.show()


# In[24]:


#histogram of pvalues
#for A vs E

A_vs_E['pvalue'].plot.hist(bins = 20, color = 'lightblue', edgecolor = 'black')
plt.title('Distribution of Statistical Significance')
plt.xlabel('p-value')
plt.ylabel('frequency')
plt.show()


# In[25]:


#histogram of pvalues
#for A vs F

A_vs_F['pvalue'].plot.hist(bins = 20, color = 'darkblue', edgecolor = 'black')
plt.title('Distribution of Statistical Significance')
plt.xlabel('p-value')
plt.ylabel('frequency')
plt.show()


# In[89]:


#making subsets and organising data to use in plotting the heatmap
#creating a subset DataFrame with gene IDs and significant log2FC
A_vs_E_significant = A_vs_E[
    (A_vs_E['log2FoldChange'].notna()) &
    (A_vs_E['log2FoldChange'].abs() >= log2FC_threshold)]

#keeping only gene_id and log2FoldChange columns
A_vs_E_significant = A_vs_E_significant[['"gene_id', 'log2FoldChange']]

#sorting by absolute log2FC
A_vs_E_significant = A_vs_E_significant.reindex(
    A_vs_E_significant['log2FoldChange'].abs().sort_values(ascending=False).index)

print(A_vs_E_significant.head())


# In[90]:


#creating a subset DataFrame with gene IDs and significant log2FC
A_vs_F_significant = A_vs_F[
    (A_vs_F['log2FoldChange'].notna()) &
    (A_vs_F['log2FoldChange'].abs() >= log2FC_threshold)]

#keeping only gene_id and log2FC columns
A_vs_F_significant = A_vs_F_significant[['"gene_id', 'log2FoldChange']]

#sorting by absolute log2FC
A_vs_F_significant = A_vs_F_significant.reindex(
    A_vs_F_significant['log2FoldChange'].abs().sort_values(ascending=False).index)

print(A_vs_F_significant.head())


# In[100]:


#merging the significant gene sets using gene_id
merged_sets = pd.merge(
    A_vs_E_significant,
    A_vs_F_significant,
    on='"gene_id',
    how='outer')

#renaming columns for clarity
merged_sets = merged_sets.rename(columns={
    'log2FoldChange_x': 'A_vs_E',
    'log2FoldChange_y': 'A_vs_F'})

#preparing data for the heatmap
heatmap_data = merged_sets.set_index('"gene_id')[['A_vs_E', 'A_vs_F']]

sns.heatmap(heatmap_data, cmap = 'viridis')

