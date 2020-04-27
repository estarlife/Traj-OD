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
import writeToCSV
import copy
import sys
import scipy.spatial.distance as DIST
import scipy.cluster.hierarchy as HAC
import random
from sklearn import metrics
import operator
import utils
import clustering_worker


class Point(object):
	def __init__(self,_x,_y):
		self.x = _x
		self.y = _y


def extractTrajectoriesUntilOD(data, originTS, originLatitude, originLongtitude, endTS, endLatitude, endLongtitude, show = True, save = False, clean = False, fname = "", path = "plots"):
	"""
	returns: OD_trajectories: in x,y coordinate;
			 OD_trajectories_lat_lon: in lat, lon coordinate;
	"""
	
	maxSpeed = 0
	for i in range(0, data.shape[0]):
		speed_over_ground = data[i][utils.dataDict["speed_over_ground"]]
		if(speed_over_ground > maxSpeed and speed_over_ground != 102.3): #1023 indicates speed not available
			maxSpeed = speed_over_ground
	
	OD_trajectories = [] # origin destination endpoints trajectory
	i = 0
	while(i< data.shape[0]):
		cur_pos = data[i]
		if(utils.nearOrigin( \
			originLatitude, \
			originLongtitude, \
			cur_pos[utils.dataDict["latitude"]], \
			cur_pos[utils.dataDict["longitude"]], \
			thresh = 0.0) and \
			cur_pos[utils.dataDict["ts"]] == originTS): # must be exact point

			this_OD_trajectory = []
			this_OD_trajectory.append(cur_pos)
			i += 1
			while(i < data.shape[0] and \
				(not utils.nearOrigin( \
					endLatitude, \
					endLongtitude, \
					data[i][utils.dataDict["latitude"]], \
					data[i][utils.dataDict["longitude"]], \
					thresh = 0.0))):
				this_OD_trajectory.append(data[i])
				i += 1
			if(i < data.shape[0]):
				this_OD_trajectory.append(data[i]) # append the destination endpoint
			this_OD_trajectory = np.asarray(this_OD_trajectory) # make it to be an np 2D array

			""" box/radius approach in cleaning of points around origin"""
			j = 1
			while(j < this_OD_trajectory.shape[0] and \
				utils.nearOrigin( \
					originLatitude, \
					originLongtitude, \
					this_OD_trajectory[j][utils.dataDict["latitude"]], \
					this_OD_trajectory[j][utils.dataDict["longitude"]], \
					thresh = utils.NEIGHBOURHOOD_ORIGIN)):
				j += 1
			this_OD_trajectory_around_origin = this_OD_trajectory[0:j]

			"""Take the box mean, treat timestamp as averaged as well"""
			this_OD_trajectory_mean_origin = boxMeanTrajectoryPoints(this_OD_trajectory_around_origin, originLatitude, originLongtitude)

			OD_trajectories.append(np.insert(this_OD_trajectory[j:],0,this_OD_trajectory_mean_origin, axis = 0))
			break  # only one trajectory per pair OD, since OD might be duplicated
		i += 1

	OD_trajectories = np.array(OD_trajectories)
	OD_trajectories_lat_lon = copy.deepcopy(OD_trajectories)
	for i in range(0, len(OD_trajectories)):
		for j in range(0, len(OD_trajectories[i])):
			x, y = utils.LatLonToXY(originLatitude, originLongtitude, OD_trajectories[i][j][utils.dataDict["latitude"]], OD_trajectories[i][j][utils.dataDict["longitude"]])
			OD_trajectories[i][j][utils.data_dict_x_y_coordinate["y"]] = y
			OD_trajectories[i][j][utils.data_dict_x_y_coordinate["x"]] = x
		# plotting purpose
		plt.scatter(OD_trajectories[i][0:len(OD_trajectories[i]),utils.data_dict_x_y_coordinate["x"]], \
			OD_trajectories[i][0:len(OD_trajectories[i]),utils.data_dict_x_y_coordinate["y"]])
	if(not plt.gca().yaxis_inverted()):
		plt.gca().invert_yaxis()
	if(save):
		plt.savefig("./{path}/{fname}.png".format(path = path, fname = fname))
	if(show):
		plt.show()
	if(clean):
		plt.clf()

	return OD_trajectories, OD_trajectories_lat_lon

def getDistance(point1, point2):
	dx, dy = utils.LatLonToXY(point1[utils.dataDict["latitude"]], point1[utils.dataDict["longitude"]], point2[utils.dataDict["latitude"]], point2[utils.dataDict["longitude"]])
	return (np.linalg.norm([dx,dy],2))

def alreadyInEndpoints(endpoints, target):
	for i in range(0, len(endpoints)):
		if(getDistance(endpoints[i], target) < utils.NEIGHBOURHOOD_ENDPOINT):
			return True
	return False

def extractEndPoints(data):
	"""
	Note: if the trajectory is discontinued because out of detection range, add that last point before out of range, and the new point in range as end point as well
	TODO: further cleaning of data is needed to extract better end points, eg. 8514019.csv end point 1,2 are actually of the same place but 3 is added due to error point
	"""
	endpoints = []

	if (len(data) > 0):
		endpoints.append(data[0]) # assume first point is an endpoint
	i = 0
	while(i< data.shape[0]):
		start_point = data[i]
		start_index = i
		
		"""Find the next_point that marks the departure from endpoint"""
		while(i+1<data.shape[0]):
			next_point = data[i+1]
			"""
			If inter point distance > thresh and is not error signal (speed is indeed> 0)
			Or
			inter point time difference > thesh
			"""
			if((getDistance(start_point, next_point) > utils.NEIGHBOURHOOD_ENDPOINT \
				and next_point[utils.dataDict["speed_over_ground"]] > 0) or \
				(next_point[utils.dataDict["ts"]] - start_point[utils.dataDict["ts"]] > utils.BOUNDARY_TIME_DIFFERENCE) \
				and i == start_index # immediate point after start point
				):
				break;
			i += 1

		next_point = data[i] # back track to get the last data point that is still near start_point
		if(i - start_index > 0 and next_point[utils.dataDict["ts"]] - start_point[utils.dataDict["ts"]] > utils.STAYTIME_THRESH):
			if(len(endpoints) == 0 or (not (endpoints[len(endpoints) - 1] == start_point).all())): # if not just appended
				endpoints.append(start_point)
		elif((i+1) != data.shape[0]): # check boundary case
			"""TODO: is there a boundary informaiton on the area that AIS can detect?"""
			next_point_outside_neighbour = data[i+1]
			if(next_point_outside_neighbour[utils.dataDict["ts"]] - start_point[utils.dataDict["ts"]] > utils.BOUNDARY_TIME_DIFFERENCE and \
				(next_point_outside_neighbour[utils.dataDict["speed_over_ground"]] != 0 or \
				start_point[utils.dataDict["speed_over_ground"]] != 0)): # if start of new trajectory at a new position after some time, boundary case, (one of the speed should not be zero)

				if(len(endpoints) == 0 or (not (endpoints[len(endpoints) - 1] == next_point).all())): # if not just appended
					endpoints.append(next_point)
				if(len(endpoints) == 0 or (not (endpoints[len(endpoints) - 1] == next_point_outside_neighbour).all())): # if not just appended
					endpoints.append(next_point_outside_neighbour)

		elif((i+1) == data.shape[0]):
			if(len(endpoints) == 0 or (not (endpoints[len(endpoints) - 1] == next_point).all())): # if not just appended
				endpoints.append(next_point) # last point in the .csv record, should be an end point

		i += 1
	
	return endpoints

def convertListOfTrajectoriesToLatLon(originLatitude, originLongtitude, listOfTrajectories):
	for i in range(0, len(listOfTrajectories)):
		for j in range(0, len(listOfTrajectories[i])):
			lat, lon = utils.XYToLatLonGivenOrigin(originLatitude, originLongtitude, listOfTrajectories[i][j][utils.data_dict_x_y_coordinate["x"]], listOfTrajectories[i][j][utils.data_dict_x_y_coordinate["y"]])
			listOfTrajectories[i][j][utils.dataDict["latitude"]] = lat
			listOfTrajectories[i][j][utils.dataDict["longitude"]] = lon
	return listOfTrajectories

def convertListOfTrajectoriesToXY(originLatitude, originLongtitude, listOfTrajectories):
	for i in range(0, len(listOfTrajectories)):
		for j in range(0, len(listOfTrajectories[i])):
			x, y = utils.LatLonToXY(originLatitude, originLongtitude, listOfTrajectories[i][j][utils.dataDict["latitude"]], listOfTrajectories[i][j][utils.dataDict["longitude"]])
			listOfTrajectories[i][j][utils.data_dict_x_y_coordinate["y"]] = y
			listOfTrajectories[i][j][utils.data_dict_x_y_coordinate["x"]] = x
	return listOfTrajectories


def endPointMatchTrajectoryCentroid(endpoint, centroid, reference_lat, reference_lon):
	assert (len(centroid) > 0), "cluster centroid must be non empty"
	x, y = utils.LatLonToXY(reference_lat,reference_lon,endpoint[utils.dataDict["latitude"]], endpoint[utils.dataDict["longitude"]])
	centroid_start_x = centroid[0][utils.data_dict_x_y_coordinate["x"]]
	centroid_start_y = centroid[0][utils.data_dict_x_y_coordinate["y"]]
	if (np.linalg.norm([x - centroid_start_x, y - centroid_start_y], 2) < 20 * utils.NEIGHBOURHOOD_ENDPOINT):
		return True
	else:
		return False


def endPointsToRepresentativeTrajectoryMapping(endpoints, trajectories, cluster_label, reference_lat, reference_lon):
	"""
	trajectories: in XY coordinate by reference_lat, reference_lon
	endpoints: in lat, lon
	cluster_label: array of cluster label w.r.t array of trajectories, starting with cluster index 1
	"""
	endpoints_cluster_dict = {}
	class_trajectories_dict = clustering_worker.formClassTrajectoriesDict(cluster_label = cluster_label, data = trajectories)
	cluster_centroids_dict = {} # [cluster:centroid] dictionary
	for class_label, trajectories in class_trajectories_dict.iteritems():
		cluster_centroids_dict[class_label] =clustering_worker.getMeanTrajecotoryWithinClass(trajectories)

	for endpoint in endpoints:
		if (not "{lat}_{lon}".format(lat = endpoint[utils.dataDict["latitude"]], \
				lon = endpoint[utils.dataDict["longitude"]]) in endpoints_cluster_dict):
			endpoints_cluster_dict["{lat}_{lon}".format(lat = endpoint[utils.dataDict["latitude"]], \
				lon = endpoint[utils.dataDict["longitude"]])] = []
		for cluster, centroid in cluster_centroids_dict.iteritems():
			if (endPointMatchTrajectoryCentroid(endpoint, centroid, reference_lat, reference_lon)):
				endpoints_cluster_dict["{lat}_{lon}".format(lat = endpoint[utils.dataDict["latitude"]], \
				lon = endpoint[utils.dataDict["longitude"]])].append(utils.ClusterCentroidTuple(cluster = cluster - 1, centroid = centroid)) # offset by 1

	return endpoints_cluster_dict


def lookForEndPoints(endpoints, endpoint_str):
	for endpoint in endpoints:
		if ("{lat}_{lon}".format(lat = endpoint[utils.dataDict["latitude"]], \
				lon = endpoint[utils.dataDict["longitude"]]) == endpoint_str):
			return endpoint
	return None

	
