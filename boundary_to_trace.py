import z3

test = [[0, 100], [2, 80], [3, 60]]
def pre_process(boundary, eps, num): #boundary is a 2-d array: [[x_1, y_1],[x_2,y_2], ...]. Normally, x is tau and y is h.
    size = len(boundary)
    b_plus = [[0, 0] for i in range(len(boundary))]
    b_minus = [[0, 0] for i in range(len(boundary))]
    for i in range(size):
        b_plus[i][0] = boundary[i][0]
        b_plus[i][1] = boundary[i][1] + eps
        b_minus[i][0] = boundary[i][0]
        b_minus[i][1] = boundary[i][1] - eps
    return b_plus, b_minus

def make_psi(b_plus, b_minus, vs, end_time):
    psi_plus = True
    psi_minus = True
    t_plus = [b_plus[i][0] for i in range(len(b_plus))]
    h_plus = [b_plus[i][1] for i in range(len(b_plus))]
    t_minus = [b_minus[i][0] for i in range(len(b_minus))]
    h_minus = [b_minus[i][1] for i in range(len(b_minus))]

    for i in range(len(t_plus)):
        temp = True
        for j in range(t_plus[i], end_time):
            temp = z3.And(temp, vs[j] < h_plus[i])
        psi_plus = z3.And(psi_plus, temp)

    for i in range(len(t_minus)):
        temp = True
        for j in range(t_minus[i], end_time):
            temp = z3.And(temp, vs[j] < h_minus[i])
        psi_minus = z3.And(psi_minus, z3.Not(temp))
    psi = z3.And(psi_plus, psi_minus)
    return psi

def output_trace(phi, psi):
    S = z3.Solver()
    S.add(z3.And(phi, psi))
    print(S.check())
    print(S.model())

def trace(boundary, eps, num, vs, end_time, phi):
    b_plus, b_minus = pre_process(boundary, eps, num)
    psi = make_psi(b_plus, b_minus, vs, end_time)
    output_trace(phi, psi)


