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

def trajectoryDissimilarityL2(t1, t2):
		# t1,t2 in  X, Y form
	i = 0
	j = 0
	dissimilarity = 0.0
	while(i < len(t1) and j < len(t2)):
		dissimilarity += DIST.euclidean([t1[i][utils.data_dict_x_y_coordinate["x"]] ,t1[i][utils.data_dict_x_y_coordinate["y"]]], \
			[t2[j][utils.data_dict_x_y_coordinate["x"]], t2[j][utils.data_dict_x_y_coordinate["y"]]])
		i += 1
		j += 1
	# only one of the following loops will be entered
	while(i < len(t1)):
		dissimilarity += DIST.euclidean([t1[i][utils.data_dict_x_y_coordinate["x"]], t1[i][utils.data_dict_x_y_coordinate["y"]]], \
		[t2[j - 1][utils.data_dict_x_y_coordinate["x"]], t2[j - 1][utils.data_dict_x_y_coordinate["y"]]]) # j -1 to get the last point in t2
		i += 1

	while(j < len(t2)):
		dissimilarity += DIST.euclidean([t1[i - 1][utils.data_dict_x_y_coordinate["x"]], t1[i - 1][utils.data_dict_x_y_coordinate["y"]]], \
			[t2[j][utils.data_dict_x_y_coordinate["x"]], t2[j][utils.data_dict_x_y_coordinate["y"]]])
		j += 1
	return dissimilarity

# vertical distance dissimilarity(project from t2 to t1)
def trajectoryDissimilarityVertical(t1, t2):
		# t1,t2 in  X, Y form
	vec1=[t2[0][0]-t1[0][0],t2[0][1]-t1[0][1]] #t1.start to t2.start
	vec2=[t2[1][0]-t1[0][0],t2[1][1]-t1[0][1]] #t1.start to t2.end
	vec3=[t1[1][0]-t1[0][0],t1[1][1]-t1[0][1]] #t1.start to t1.end
	unit=normalized(vec3)

	l12=magnitude(vec1)**2-dot(vec1,unit)**2
	l22=magnitude(vec2)**2-dot(vec2,unit)**2
	dissimilarity = (l1+l2)/(sqrt(l12)+sqrt(l22))
	
	return dissimilarity

# parallel distance dissimilarity(project from t2 to t1)
def trajectoryDissimilarityParallel(t1, t2):
		# t1,t2 in  X, Y form
	vec1=[t2[0][0]-t1[0][0],t2[0][1]-t1[0][1]] #t1.start to t2.start
	vec2=[t2[1][0]-t1[0][0],t2[1][1]-t1[0][1]] #t1.start to t2.end
	vec3=[t2[0][0]-t2[0][0],t2[0][1]-t2[0][1]] #t1.end to t2.start
	vec4=[t2[1][0]-t1[0][0],t2[1][1]-t1[0][1]] #t1.end to t2.end
	
	vec5=[t1[1][0]-t1[0][0],t1[1][1]-t1[0][1]] #t2.start to t2.end
	unit=normalized(vec5)

	dissimilarity=min(abs(dot(vec1,unit)),abs(dot(vec2,unit)),abs(dot(vec3,unit)),abs(dot(vec4,unit)))
	return dissimilarity

# angular distance dissimilarity(project from t2 to t1)
def trajectoryDissimilarityAngular(t1, t2):
		# t1,t2 in  X, Y form
	vec1=[t1[1][0]-t1[0][0],t1[1][1]-t1[0][1]] #t1.start to t1.end
	vec2=[t2[1][0]-t2[0][0],t2[1][1]-t2[0][1]] #t2.start to t2.end
	unit=normalized(vec1)

	dissimilarity=sqrt(magnitude(vec2)**2-dot(vec2,unit)**2)
	return dissimilarity

def distanceBetweenTwoPoint(p1, p2):


	return DIST.euclidean( \
		[p1[utils.data_dict_x_y_coordinate["x"]], p1[utils.data_dict_x_y_coordinate["y"]]], \
		[p2[utils.data_dict_x_y_coordinate["x"]], p2[utils.data_dict_x_y_coordinate["y"]]])

def getTrajectoryLength(t):
	"""
	t: trajectory in X, Y coordinates
	return the sum of length of the trajectory
	"""
	distance = 0.0
	cur_pos = t[0]
	for i in range(1, len(t)):
		next_pos = t[i]
		distance += distanceBetweenTwoTrajectoryPoint(cur_pos, next_pos)
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

	






	