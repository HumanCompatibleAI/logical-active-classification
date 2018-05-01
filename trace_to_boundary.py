import numpy as np 
import pandas as pd
import funcy as fn
import multidim_threshold as mdt 
import operator
from functools import reduce
import stl
from stl.load import from_pandas
from scipy.spatial.distance import euclidean
from fastdtw import fastdtw
import matplotlib.pyplot as plt

def dtw_dist(x, y):
    x = from_pandas(x)
    y = from_pandas(y)
    return fastdtw(x['Y'].sample(0.1), y['Y'].sample(0.1), dist=euclidean)[0]


rectangle = mdt.to_rec([(0, 1), (0, 19.799)])

import funcy
@fn.autocurry
def phi(x, params):
    h, tau = params
    x= x['Y'].slice(tau, None)
    return all(map(lambda y:y[1] <= h, x))

def compute_boundary(trace, eps=0.1):
    refinements = mdt.volume_guided_refinement([rectangle], phi(trace))
    return list(fn.pluck(1, fn.first(fn.dropwhile(lambda x :-min(fn.pluck(0, x))> eps, refinements))))

def boundary(trace):
	np_boundary = np.array(trace)
	df = pd.DataFrame(np_boundary, columns=['X', 'Y']).set_index('X')
	t = from_pandas(df)
	bound = compute_boundary(t, 0.01)
	top = sorted([r.top for r in bound])
	t_copy = list(top)
	bottom = sorted([r.bot for r in bound])
	b_copy = list(bottom)
	size = len(top)
	for i in range(size):
		top[i] = t_copy[size - i - 1]
		top[i][0], top[i][1] = top[i][1], top[i][0]
		bottom[i]  = b_copy[size - i - 1]
		bottom[i][0], bottom[i][1] = bottom[i][1], bottom[i][0] 
	boundary = [[0, 0] for i in range(size)]
	for i in range(size):
		boundary[i][0] = bottom[i][0]
		boundary[i][1] = top[i][1]
	boundary.append([top[size-1][0], bottom[size-1][1]])
	return boundary
	#print(top)
	#print(bottom)
	#print(boundary)
	#x_axis = [boundary[i][0] for i in range(len(boundary))]
	#y_axis = [boundary[i][1] for i in range(len(boundary))]
	#plt.plot(x_axis, y_axis)
	#plt.show()





