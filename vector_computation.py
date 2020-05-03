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
from math import acos,pi
from math import sqrt
from decimal import Decimal,getcontext


# 向量的大小
def magnitude(vec):
    return sqrt(vec[0] ** 2 + vec[1] ** 2)

# 标量乘法
def timesscalar(vec, m):
    return [m*vec[0] , m*vec[1]]

# 单位向量
def normalized(vec):
    return timesscalar(vec, 1.0/magnitude(vec))

# 两个向量的点积
def dot(vec,v):
    return vec[0]*v[0]+vec[1]*v[1]

# 两个向量之间的角度
def angle_with(vec1, vec2):
    u1 = normalized(vec1)
    u2 = normalized(vec2)
    dots = dot(u1,u2)
    if abs(abs(dots) - 1) < 1e-10:
        if dots < 0:
            dots = -1
        else:
            dots = 1
    angle_in_radians = acos(dots)
    return angle_in_radians

# vec1在vec2上的投影
def projection(vec1, vec2):
    if vec2==0:
        return vec1
    unit = normalized(vec2)
    weight = dot(vec1,unit)
    return timesscalar(unit, weight)