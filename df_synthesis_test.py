import z3

traj_length = 20

xs = z3.IntVector('xs', traj_length)
ys = z3.IntVector('ys', traj_length)
t  = z3.Int('t')
# should probably limit t to be within traj_length, but meh

right_start = z3.And(xs[0] == 0, ys[0] == 0)

no_jump_vert = z3.ForAll(t, z3.Or([xs[t+1] == xs[t], xs[t+1] == xs[t] + 1,
                                   xs[t+1] == xs[t] - 1]))

no_jump_horiz = z3.ForAll(t, z3.Or([ys[t+1] == ys[t], ys[t+1] == ys[t] + 1,
                                    ys[t+1] == ys[t] - 1]))

def ever_below(time, y_val):
    return z3.Exists(t, z3.And([t >= 0, t <= time, ys[t] < y_val]))

def eventually_up_right(time, x_val, y_val):
    return z3.ForAll(t, z3.Or(t < time, z3.And(xs[t] >= -x_val,
                                               ys[t] >= -y_val)))

def synthesise_me_a_trajectory(time1, time2, y_val1, x_val2, y_val2):
    S = z3.Solver()
    S.add(z3.And(right_start, no_jump_vert, no_jump_horiz))
    S.add(ever_below(time1, y_val1))
    S.add(eventually_up_right(time2, x_val2, y_val2))
    if S.check() == z3.sat:
        model = S.model()
        x_vals = model.eval(xs)
        y_vals = model.eval(ys)
        print(x_vals, y_vals)
    else:
        print("No trajectory like this exists!")

if __name__ == '__main__':
    synthesise_me_a_trajectory(5,15,-2,6,6)
