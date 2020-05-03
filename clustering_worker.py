import numpy as np
import math
import scipy.io as sio
import scipy.sparse as sp
from scipy.sparse import csc_matrix
import csv
import matplotlib.pyplot as plt
import datetime
import time
import os
from scipy import interpolate
import scipy.spatial.distance as distance
from collections import OrderedDict
import copy
import sys
import scipy.spatial.distance as DIST
import scipy.cluster.hierarchy as HAC
import random
from sklearn import metrics
import operator
import utils
import vector_computation
from math import acos,pi
from math import sqrt


# vertical distance dissimilarity(project from t2 to t1)
def trajectoryDissimilarityVertical(t1, t2):
		# t1,t2 in  X, Y form 
	vec1=[t2[0][0]-t1[0][0],t2[0][1]-t1[0][1]] #t1.start to t2.start
	vec2=[t2[1][0]-t1[0][0],t2[1][1]-t1[0][1]] #t1.start to t2.end
	vec3=[t1[1][0]-t1[0][0],t1[1][1]-t1[0][1]] #t1.start to t1.end
	if vec3==[0,0]:
		return max(vector_computation.magnitude(vec1),vector_computation.magnitude(vec2))
		
	unit=vector_computation.normalized(vec3)

	l12=vector_computation.magnitude(vec1)**2-vector_computation.dot(vec1,unit)**2
	l22=vector_computation.magnitude(vec2)**2-vector_computation.dot(vec2,unit)**2
	if abs(l12)<0.000001: # prevent l12 l22 from being negative
		l12=0
	if abs(l22)<0.000001:
		l22=0

	if l12+l22==0:
		dissimilarity=0
	else:
		dissimilarity = (l12+l22)/(sqrt(l12)+sqrt(l22))
	
	return dissimilarity

# parallel distance dissimilarity(project from t2 to t1)
def trajectoryDissimilarityParallel(t1, t2):
		# t1,t2 in  X, Y form
	vec1=[t2[0][0]-t1[0][0],t2[0][1]-t1[0][1]] #t1.start to t2.start
	vec2=[t2[1][0]-t1[0][0],t2[1][1]-t1[0][1]] #t1.start to t2.end
	vec3=[t2[0][0]-t1[1][0],t2[0][1]-t1[1][1]] #t1.end to t2.start
	vec4=[t2[1][0]-t1[1][0],t2[1][1]-t1[1][1]] #t1.end to t2.end
	vec5=[t1[1][0]-t1[0][0],t1[1][1]-t1[0][1]] #t1.start to t1.end
	unit=vector_computation.normalized(vec5)
	dissimilarity=min(abs(vector_computation.dot(vec1,unit)), \
                      abs(vector_computation.dot(vec2,unit)),abs(vector_computation.dot(vec3,unit)),abs(vector_computation.dot(vec4,unit)))
	return dissimilarity

# angular distance dissimilarity(project from t2 to t1)
def trajectoryDissimilarityAngular(t1, t2):
		# t1,t2 in  X, Y form
	vec1=[t1[1][0]-t1[0][0],t1[1][1]-t1[0][1]] #t1.start to t1.end
	vec2=[t2[1][0]-t2[0][0],t2[1][1]-t2[0][1]] #t2.start to t2.end
	if vec1==[0,0] or vec2==[0,0]:
		return 0
	unit=vector_computation.normalized(vec1)
	tmp=vector_computation.magnitude(vec2)**2-vector_computation.dot(vec2,unit)**2
	if abs(tmp)<0.000001:
		tmp=0
	dissimilarity=sqrt(tmp)
	return dissimilarity

def distanceBetweenTwoPoint(p1, p2):


	return DIST.euclidean( \
		[p1[0], p1[1]], \
		[p2[0], p2[1]])

def getTrajectoryLength(t):
	"""
	t: trajectory in X, Y coordinates
	return the sum of length of the trajectory
	"""
	distance = 0.0
	cur_pos = t[0]
	for i in range(1, len(t)):
		next_pos = t[i]
		distance += distanceBetweenTwoPoint(cur_pos, next_pos)
		cur_pos = next_pos # forward current position

	return distance

def TrajectorytoVec(t):
	"""
	t: trajectory in X,Y coordinates
	return the (dx,dy) as displacement of the trajectory t
	"""
	return np.asarray([ \
		t[len(t) - 1][utils.data_dict_x_y_coordinate["x"]] - t[0][utils.data_dict_x_y_coordinate["x"]], \
		t[len(t) - 1][utils.data_dict_x_y_coordinate["y"]] - t[0][utils.data_dict_x_y_coordinate["y"]]
		])

def formClassTrajectoriesDict(cluster_label, data):
	"""
	cluster_label: length n
	data: length n
	returns: a dictionary of [cluster_label: [trajectories]]
	"""
	assert len(cluster_label) == len(data), "data and cluster_label length should be the same"
	class_trajectories_dict = {} # a dictionary of class label to its cluster
	for i in range(0, len(cluster_label)):
		class_label = cluster_label[i]
		if(not class_label in class_trajectories_dict):
			class_trajectories_dict[class_label] = []
		class_trajectories_dict[class_label].append(data[i])
	return class_trajectories_dict

	






	