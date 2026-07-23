	
import sys
import os
import pandas as pd
# from matplotlib import pyplot as plt

f = open('radrathe_QD.dat', 'r')

df = pd.read_csv(f, sep='\s+', names=["ni",'li','si','nf','lf','sf', "Energy",'AValue'])

df['LLi'] = df.loc[:, 'li']
df['LLf'] = df.loc[:, 'lf']
df['Vali'] = 1
df['Valf'] = 1

# df.insert(1, "LLi", 1)
print (df)

df = df[["ni",'li','LLi','si','Vali','nf','lf','LLf','sf','Vali', "Energy",'AValue']]
print (df)
df.to_csv('radrathe.dat',sep=" ", header = False, index = False)  