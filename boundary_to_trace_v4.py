import z3
import math
import numpy as np
import matplotlib.pyplot as plt
import funcy as fn


import funcy
@fn.autocurry
def phi(x, params):
    h, tau = params
    t_plus.slice(tau, None)
    return all(map(lambda y:y[1] <= h, x))

def pre_process(boundary, eps): #boundary is a 2-d array: [[x_1, y_1],[x_2,y_2], ...]. Normally, x is tau and y is h.
    size = len(boundary)
    b_plus = [[0, 0] for i in range(len(boundary)-1)]
    b_minus = [[0, 0] for i in range(len(boundary)-1)]
    for i in range(size - 1):
        dx = boundary[i][0] - boundary[i+1][0]
        dy = boundary[i][1] - boundary[i+1][1]
        norm = math.sqrt(dx**2 + dy**2)
        dx = dx / norm * eps
        dy = dy / norm * eps
        b_plus[i][0] = boundary[i][0] + dy
        b_plus[i][1] = boundary[i][1] - dx
        b_minus[i][0] = boundary[i][0] - dy
        b_minus[i][1] = boundary[i][1] + dx

    return b_plus, b_minus

def non_perp_process(boundary, eps):
    size = len(boundary)
    b_plus = [[0, 0] for i in range(len(boundary))]
    b_minus = [[0, 0] for i in range(len(boundary))]
    for i in range(size - 1):
        b_plus[i][0] = boundary[i][0]
        b_plus[i][1] = boundary[i][1] + eps
        b_minus[i][0] = boundary[i][0] 
        b_minus[i][1] = boundary[i][1] - eps

    return b_plus, b_minus
     
def make_psi(b_plus, b_minus, vs, end_time, num):
    psi_plus = True
    psi_minus = True
    t_plus = [int(b_plus[i][0]) for i in range(len(b_plus))]
    h_plus = [b_plus[i][1] for i in range(len(b_plus))]
    t_minus = [int(b_minus[i][0]) for i in range(len(b_minus))]
    h_minus = [b_minus[i][1] for i in range(len(b_minus))]


    for i in range(len(t_plus)):
        temp = True
        for j in range(int(t_plus[i] * 1.0 / end_time * num), num):
            temp = z3.And(temp, vs[j] < h_plus[i])
        psi_plus = z3.And(psi_plus, temp)

    for i in range(len(t_minus)):
        temp = True
        for j in range(int(t_minus[i] * 1.0 / end_time * num), num):
            temp = z3.And(temp, vs[j] < h_minus[i])

        psi_minus = z3.And(psi_minus, z3.Not(temp))
    psi = z3.And(psi_plus, psi_minus)
    return psi

def output_trace(phi, psi, vs, num):
    S = z3.Solver()
    S.add(z3.And(phi, psi))
    S.check()
    m = S.model()
    y = [i for i in range(num)]
    x = [i * 1.0 / 5 for i in y]
    for i in range(num):
        y[i] = m[vs[i]].numerator_as_long() * 1.0 / m[vs[i]].denominator_as_long()
    trace = [[x[i], y[i]] for i in range(len(x))]
    return trace
    #print(y)
    #plt.plot(x, y)
    #plt.show()
    

def trace(boundary, eps, num, vs, end_time, phi):
    b_plus, b_minus = pre_process(boundary, eps)
    psi = make_psi(b_plus, b_minus, vs, end_time, num)
    return output_trace(phi, psi, vs, num)


