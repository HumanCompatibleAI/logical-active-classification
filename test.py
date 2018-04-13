import z3
import boundary_to_trace as bt
num = 20

xs = ['x%d' % i for i in range(num)]
us = ['u%d' % i for i in range(num)]

ys = ['y%d' % i for i in range(num)]
vs = ['v%d' % i for i in range(num)]
for i in range(num):
    xs[i] = z3.Real('x%d' % i)
    us[i] = z3.Real('u%d' % i)
    ys[i] = z3.Real('y%d' % i)
    vs[i] = z3.Real('v%d' % i)
A = 1
B = 1
#A = [[ Real('a_%d_%d' % (i, j)) for j in range(dim)] for i in range(dim)]
#B = [[ Real('b_%d_%d' % (i, j)) for j in range(dim)] for i in range(dim)]

def make_phi(z, w):
    formula = True
    for i in range(num - 1):
        formula = z3.And(formula, z[i+1] == A * z[i] + B * w[i], z[i+1] != z[i])
    return formula
phi = z3.And(make_phi(xs, us), make_phi(ys, vs))

tp = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
hp = [75,75,75,75,75,70,65,60,55,50,45,40,35,30,25,15,10,10,10,10]

boundary = [[tp[i],hp[i]-0.5] for i in range(20)]
eps = 0.5
num = 20
end_time = 20
values = ['z%d' % i for i in range(num)]

for i in range(num):
	values[i] = xs[i] - ys[i]
#print(boundary)
#print(bt.pre_process(boundary, eps, num))
#b1, b2 = bt.pre_process(boundary, eps, num)
#psi = bt.make_psi(b1, b2, values, end_time)
#print(psi)
bt.trace(boundary, eps, num, values, end_time, phi)
#bt.trace(boundary, eps, num, values, end_time, phi)


