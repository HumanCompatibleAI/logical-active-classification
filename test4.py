import z3
import boundary_to_trace_v4 as bt
num = 100

# xs = ['x%d' % i for i in range(num)]
# us = ['u%d' % i for i in range(num)]

# for i in range(num):
#     xs[i] = z3.Real('x%d' % i)
#     us[i] = z3.Real('u%d' % i)

xs = [z3.Real('x%d' % i) for i in range(num)]
us = [z3.Real('u%d' % i) for i in range(num)]

A = 1
B = 1
#A = [[ Real('a_%d_%d' % (i, j)) for j in range(dim)] for i in range(dim)]
#B = [[ Real('b_%d_%d' % (i, j)) for j in range(dim)] for i in range(dim)]

def make_phi(z, w):
    formula = True
    for i in range(num - 1):
    	formula = z3.And(formula, z[i+1] == A * z[i] + B * w[i], z3.And(z[i+1] - z[i] <= 0.05, z[i] - z[i+1] <= 0.05), z[i+1] != z[i], z[i] >= 0, z[i+1] >= 0)
    	#formula = z3.And(formula, z[i+1] == A * z[i] + B * w[i], z3.Or(z[i+1] - z[i] > 0.05, z[i+1] - z[i] < -0.05), z[i] >= 0)
        #formula = z3.And(formula, z[i+1] == A * z[i] + B * w[i], z[i+1] != z[i], z[i] >= 0, z3.Or(w[i] > 0.8, w[i] < -0.8))
    return formula
phi = make_phi(xs, us)

tp = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
hp1 = [0.725,0.725,0.725,0.725,0.725,0.61,0.53,0.38,0.15,0.12,0.12,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09]
hp2 = [0.8, 0.8, 0.72, 0.65, 0.6, 0.52, 0.37, 0.25, 0.18, 0.12, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
hp3 = [0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.74, 0.72, 0.71]

boundary = [[tp[i],hp1[i]] for i in range(20)]
boundary2 = [[tp[i],hp2[i]] for i in range(20)]
boundary3 = [[tp[i],hp3[i]] for i in range(20)]
eps = 0.005
end_time = 20

#print(boundary)
#print(bt.pre_process(boundary, eps, num))
#b1, b2 = bt.pre_process(boundary, eps, num)
#psi = bt.make_psi(b1, b2, values, end_time)
#print(psi)

#a, b = bt.pre_process(boundary, eps)
#bt.make_psi(a, b, xs, end_time, num)
trace = bt.trace(boundary, eps, num, xs, end_time, phi)

#boundary = tb.boundary(trace)


#print(m)
#x = []
#y = []
#z = []
#for i in range(20):
#	x[i] = float(m[xs[i]].as_decimal(10))
#	y[i] = float(m[ys[i]].as_decimal(10))
#	z[i] = x[i] - y [i]
#print(z)



