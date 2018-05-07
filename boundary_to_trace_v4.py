import z3
import math
import numpy as np
import funcy as fn
import utils

import funcy
@fn.autocurry
def phi(x, params):
    h, tau = params
    t_plus.slice(tau, None)
    return all(map(lambda y:y[1] <= h, x))

def pre_process(boundary, eps, num, end_time): #boundary is a 2-d array: [[x_1, y_1],[x_2,y_2], ...]. Normally, x is tau and y is h. 
    size = len(boundary)
    if (size != 0 and boundary[0][0] == 0 and boundary[size-1][0] == 0):
        for i in range(size):
            boundary[i][0] = boundary[i][0] + (end_time / num * 1.0)
        #print(boundary)
    b_plus = [[0, 0] for i in range(len(boundary)-1)]
    b_minus = [[0, 0] for i in range(len(boundary)-1)]
    for i in range(size - 1):
        b_plus[i] = utils.move_middle_out([boundary[i], boundary[i+1]], eps)
        b_minus[i] = utils.move_middle_out([boundary[i], boundary[i+1]],
                                           (-1.0)*eps)
        # dx = boundary[i][0] - boundary[i+1][0]
        # dy = boundary[i][1] - boundary[i+1][1]
        # norm = math.sqrt(dx**2 + dy**2)
        # dx = dx / norm * eps
        # dy = dy / norm * eps
        # b_plus[i][0] = boundary[i][0] + dy
        # b_plus[i][1] = boundary[i][1] - dx
        # b_minus[i][0] = boundary[i][0] - dy
        # b_minus[i][1] = boundary[i][1] + dx
    return b_plus, b_minus

# print(pre_process([[0,1], [1,0]], 0.70710678, 1000, 1000))

def non_perp_process(boundary, eps):
    # deprecated
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
    t_plus = [int(point[0] * num / end_time * 1.0) for point in b_plus]
    # t_plus = [(point[0] * num / end_time * 1.0) for point in b_plus]
    h_plus = [b_plus[i][1] for i in range(len(b_plus))]
    t_minus = [int(b_minus[i][0] * num / end_time * 1.0) for i in range(len(b_minus))]
    # t_minus = [(b_minus[i][0] * num / end_time * 1.0) for i in range(len(b_minus))]
    h_minus = [b_minus[i][1] for i in range(len(b_minus))]
    print("t_plus", t_plus)
    # print(h_plus)
    print("t_minus", t_minus)



    for i in range(len(t_plus)):
        temp = True
        # print(int(t_plus[i]))
        for j in range(t_plus[i], num):
            temp = z3.And(temp, vs[j] < h_plus[i])
        # for j in range(num):
        #     if j > t_plus[i]:
        #         temp = z3.And(temp, vs[j] < h_plus[i])
        psi_plus = z3.And(psi_plus, temp)

    for i in range(len(t_minus)):
        temp = True
        for j in range(t_minus[i], num):
            temp = z3.And(temp, vs[j] < h_minus[i])
        # for j in range(num):
        #     if j > t_minus[i]:
        #         temp = z3.And(temp, vs[j] < h_minus[i])
        psi_minus = z3.And(psi_minus, z3.Not(temp))
    psi = z3.And(psi_plus, psi_minus)
    return psi

def output_trace(phi, psi, vs, num, end_time):
    S = z3.Solver()
    S.add(z3.And(phi, psi))
    S.check()
    m = S.model()
    y = [i for i in range(num)]
    x = [i / (num / end_time * 1.0) for i in y]
    for i in range(num):
        y[i] = m[vs[i]].numerator_as_long() * 1.0 / m[vs[i]].denominator_as_long()
    trace = [[x[i], y[i]] for i in range(len(x))]
    return trace
    #print(y)
    #plt.plot(x, y)
    #plt.show()
    

def trace(boundary, eps, num, vs, end_time, phi):
    # print("Inside trace function of boundary_to_trace")
    b_plus, b_minus = pre_process(boundary, eps, num, end_time)
    print("b_plus", b_plus)
    print("b_minus", b_minus)
    psi = make_psi(b_plus, b_minus, vs, end_time, num)
    return output_trace(phi, psi, vs, num, end_time)


