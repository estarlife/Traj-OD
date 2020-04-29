import os
import numpy as np
import math
import scipy.io as sio
import scipy.sparse as sp
from scipy.sparse import csc_matrix
import csv
import matplotlib.pyplot as plt
import datetime
import time
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
from collections import namedtuple

"""
CONSTANTS
"""

dataDict = {
"tid":0,
"ts":1,
"latitude":2,
"longitude":3
}

data_dict_x_y_coordinate = {
"tid":0,
"ts":1,
"y":2,
"x":3
}


ClusterCentroidTuple = namedtuple('ClusterCentroidTuple', ['cluster', 'centroid'])

# from Lat Lon to X Y coordinates in Km
def LatLonToXY (lat1,lon1,lat2, lon2): 
	# fix origin for display
	dx = (lon2-lon1)*40000*math.cos((lat1+lat2)*math.pi/360)/360
	dy = (lat1-lat2)*40000/360
	return dx, dy

# from X Y coordinates in Km to Lat Lon 
def XYToLatLonGivenOrigin(lat1, lon1, x, y):
	lat2 = lat1 - y*360/40000
	lon2 = lon1 + x/(40000*math.cos((lat1+lat2)*math.pi/360)/360)
	return lat2, lon2

def isErrorTrajectory(trajectory, center_lat_sg, center_lon_sg):

	if(len(trajectory) <= 1):
		return True

	for i in range(0, len(trajectory)):
		lat = trajectory[i][dataDict["latitude"]]
		lon = trajectory[i][dataDict["longitude"]]
		dx, dy = LatLonToXY (lat, lon, center_lat_sg, center_lon_sg)
		if(np.linalg.norm([dx, dy], 2) > MAX_DISTANCE_FROM_SG):
			return True
	return False

def removeErrorTrajectoryFromList(trajectories, center_lat_sg = 1.2, center_lon_sg = 103.8):

	i = 0
	while(i < len(trajectories)):
		if(isErrorTrajectory(trajectories[i], center_lat_sg, center_lon_sg)):
			if(isinstance(trajectories, list)): # if list, call the list's delete method
				trajectories.pop(i)
			elif(isinstance(trajectories, np.ndarray)): # if numpy.ndarray, call its delete method
				trajectories = np.delete(trajectories, i, 0)
		else:
			i += 1
	return trajectories

#TRAJECTORY IN LIST OR NP.NDARRAY 

def writeDataToCSV(data, path, file_name):
	"""
	path: without trailing '/'
	file_name: string of name of file, without .csv suffix
	"""
	# print path + "/" +npz_file_name[:npz_file_name.find(".")]+ ".csv"
	# raise ValueError
	dataDict = {
	"ts":0,
	"latitude":1,
	"longitude":2
	}

	# datetime.datetime.fromtimestamp(currentTS).strftime('%Y-%m-%dT%H:%M:%SZ')

	with open(path +"/"+ file_name+ ".csv", 'w') as csvfile:
		fieldnames = [
		'latitude', \
		'longitude', \
		'ts']
		writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
		
		writer.writeheader()
	
		for i in range (0, data.shape[0]):
			writer.writerow({'navigation_status': data[i][dataDict['navigation_status']], 
				'rate_of_turn':data[i][dataDict['rate_of_turn']], 
				'speed_over_ground':data[i][dataDict['speed_over_ground']], 
				'latitude':data[i][dataDict['latitude']], 
				'longitude':data[i][dataDict['longitude']], 
				'course_over_ground':data[i][dataDict['course_over_ground']], 
				'true_heading':data[i][dataDict['true_heading']], 
				'ts':data[i][dataDict['ts']],
				'ts_string':datetime.datetime.fromtimestamp(data[i][dataDict['ts']]).strftime('%Y-%m-%dT%H:%M:%SZ')
				})

	return