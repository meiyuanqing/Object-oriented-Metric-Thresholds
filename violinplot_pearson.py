#!/usr/bin/env python
# encoding:utf-8
"""
Author : Yuanqing Mei
Date : 2022/1/5
Time: 23:18
File: violinplot_pearson.py
HomePage : https://github.com/meiyuanqing
Email : dg1533019@smail.nju.edu.cn

draw violin plot of pearson for each OO metric.

"""
import matplotlib.pylab as plt
import numpy as np
import seaborn as sns
import pandas as pd

# display all columns and rows, and set the item of row of dataframe
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', 5000)

working_Directory = "F:\\MTmeta\\pearsonMeta\\"


# inverse Fisher Transformation
def inverse_Fisher_Z(fisher_Z):
    rp = (np.exp(2 * fisher_Z) - 1) / (np.exp(2 * fisher_Z) + 1)
    return rp

# df = sns.load_dataset('iris')
df = pd.read_csv(working_Directory + 'Pearson_effects.csv')
print(df.head())
cohesion = ['LCOM1', 'LCOM2', 'LCOM3', 'LCOM4', 'Co', 'NewCo', 'LCOM5', 'NewLCOM5', 'TCC', 'LCC', 'ICH', 'OCC', 'PCC',
            'DCd', 'DCi', 'CAMC', 'NHD', 'SNHD']
coupling = ['ACAIC', 'ACMIC', 'AMMIC', 'DMMEC', 'OCAEC', 'OCAIC', 'OCMEC', 'OCMIC', 'OMMEC', 'OMMIC', 'DCAEC',
            'DCMEC', 'CBI', 'CBO', 'DAC', 'ICP', 'IHICP', 'MPC', 'NIHICP', 'RFC']
inheritance = ["AID", "CLD", "DIT", "DP", "DPA", "DPD", "NMA", "NMI", "NMO", "NOA", "NOC", "NOD", "NOP", "SIX",
               "SP", "SPA", "SPD"]
size = ["NA", "NAIMP", "NM", "NMIMP", "NumPara", "SLOC", "stms"]

cohesion_1 = ['LCOM1', 'LCOM2', 'LCOM3', 'LCOM4', 'LCOM5', 'NewLCOM5']
cohesion_2 = ['Co', 'NewCo', 'TCC', 'LCC', 'ICH', 'OCC']
cohesion_3 = ['PCC', 'DCd', 'DCi', 'CAMC', 'NHD', 'SNHD']

df_cohesion_1 = df[df['metric'].isin(cohesion_1)].reset_index()
df_cohesion_2 = df[df['metric'].isin(cohesion_2)].reset_index()
df_cohesion_3 = df[df['metric'].isin(cohesion_3)].reset_index()

coupling_1 = ['ACAIC', 'ACMIC', 'AMMIC', 'DMMEC', 'OCAEC', 'OCAIC']
coupling_2 = ['OCMEC', 'OCMIC', 'OMMEC', 'OMMIC', 'DCAEC', 'DCMEC', 'CBI']
coupling_3 = ['CBO', 'DAC', 'ICP', 'IHICP', 'MPC', 'NIHICP', 'RFC']

df_coupling_1 = df[df['metric'].isin(coupling_1)].reset_index()
df_coupling_2 = df[df['metric'].isin(coupling_2)].reset_index()
df_coupling_3 = df[df['metric'].isin(coupling_3)].reset_index()

inheritance_1 = ["AID", "CLD", "DIT", "DP", "DPA", "DPD"]
inheritance_2 = ["NMA", "NMI", "NMO", "NOA", "NOC", "NOD"]
inheritance_3 = ["NOP", "SIX", "SP", "SPA", "SPD"]

df_inheritance_1 = df[df['metric'].isin(inheritance_1)].reset_index()
df_inheritance_2 = df[df['metric'].isin(inheritance_2)].reset_index()
df_inheritance_3 = df[df['metric'].isin(inheritance_3)].reset_index()

df_size = df[df['metric'].isin(size)].reset_index()

# p1 = sns.violinplot(y=inverse_Fisher_Z(df[df['metric'] == 'AID']['Fisher_Z']))
# p1 = sns.violinplot(y=inverse_Fisher_Z(df[df['metric'] == 'AID']['Fisher_Z']))
# p1 = sns.violinplot(x=df['metric'], y=inverse_Fisher_Z(df['Fisher_Z']))
# p1 = sns.violinplot(x=df_cohesion_1['metric'], y=inverse_Fisher_Z(df_cohesion_1['Fisher_Z']))
# p2 = sns.violinplot(x=df_cohesion_2['metric'], y=inverse_Fisher_Z(df_cohesion_2['Fisher_Z']))
# p3 = sns.violinplot(x=df_cohesion_3['metric'], y=inverse_Fisher_Z(df_cohesion_3['Fisher_Z']))
# plt.show()
# plt.savefig(working_Directory + 'cohesion_1.png', bbox_inches='tight')
# plt.savefig(working_Directory + 'cohesion_2.png', bbox_inches='tight')
# plt.savefig(working_Directory + 'cohesion_3.png', bbox_inches='tight')

# p1 = sns.violinplot(x=df_coupling_1['metric'], y=inverse_Fisher_Z(df_coupling_1['Fisher_Z']))
# p2 = sns.violinplot(x=df_coupling_2['metric'], y=inverse_Fisher_Z(df_coupling_2['Fisher_Z']))
# p3 = sns.violinplot(x=df_coupling_3['metric'], y=inverse_Fisher_Z(df_coupling_3['Fisher_Z']))

# plt.savefig(working_Directory + 'coupling_1.png', bbox_inches='tight')
# plt.savefig(working_Directory + 'coupling_2.png', bbox_inches='tight')
# plt.savefig(working_Directory + 'coupling_3.png', bbox_inches='tight')

# p1 = sns.violinplot(x=df_inheritance_1['metric'], y=inverse_Fisher_Z(df_inheritance_1['Fisher_Z']))
# p2 = sns.violinplot(x=df_inheritance_2['metric'], y=inverse_Fisher_Z(df_inheritance_2['Fisher_Z']))
# p3 = sns.violinplot(x=df_inheritance_3['metric'], y=inverse_Fisher_Z(df_inheritance_3['Fisher_Z']))

p4 = sns.violinplot(x=df_size['metric'], y=inverse_Fisher_Z(df_size['Fisher_Z']))

# plt.savefig(working_Directory + 'inheritance_1.png', bbox_inches='tight')
# plt.savefig(working_Directory + 'inheritance_2.png', bbox_inches='tight')
# plt.savefig(working_Directory + 'inheritance_3.png', bbox_inches='tight')

plt.savefig(working_Directory + 'size_3.png', bbox_inches='tight')
