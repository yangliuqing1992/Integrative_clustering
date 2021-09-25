# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 11:33:09 2021

@author: yangl
"""

from math import *
import pandas as pd
import numpy as np
from scipy.spatial import distance

bulk=pd.read_csv("X.csv",index_col=0)

p=bulk.shape[0]
inv_dis=np.zeros([p,p])
for i in range(p):
    for j in range((i+1),p):
        dis=distance.euclidean(bulk.iloc[i,:], bulk.iloc[j,:])
        if dis!=0:
            inv_dis[i][j]=1/dis
        else:
            inv_dis[i][j]=1
inv_dis=inv_dis+inv_dis.transpose()
inv_dis=(inv_dis-inv_dis.min())/(inv_dis.max()-inv_dis.min())
np.savetxt("Xaff.csv",inv_dis,delimiter=",")