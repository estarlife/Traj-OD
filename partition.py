import pickle
import os
import re
import sys
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math
import scipy.io as sio
import scipy.sparse as sp
from scipy.sparse import csc_matrix
import csv
import matplotlib.pyplot as plt
import datetime
import time
import utils
import os

def loadData(filename):
    pkl_file = open(filename, 'rb')
    data = pickle.load(pkl_file)

    pkl_file.close()
    #loader = np.load(filename)
    return data

def traj_partition(traj):
    all_endpoint=[]
    endpoint=[]

    
    for i in range(len(traj)):
        tmp=traj[i][0][0]    
        tmp_index=0
        endpoint=[]
        for j in range(len(traj[i])):
            
#           print(j, len(traj[i]))
            if j==len(traj[i])-1:
                endpoint.append([tmp_index,j])
            elif  traj[i][j][0]-tmp>10000:
                endpoint.append([tmp_index,j])
                tmp_index=j+1
                tmp=traj[i][j+1][0]
        all_endpoint.append(endpoint)
    return all_endpoint

def subtrajectory(traj, endpoint): #endpoint in form [traj_index, []]
    traj_index, [st, ed]=endpoint
    subtraj=traj[traj_index][st:ed+1]
    #list slice：[start，stop)，include start but not stop
    return subtraj

def subtrajectoryfromst2ed(traj, endpoint):
    traj_index, [st, ed]=endpoint
    subtraj=[traj[traj_index][st],traj[traj_index][ed]]
    return subtraj

