import z3

traj_length = 20

xs = z3.IntVector('xs', traj_length)
ys = z3.IntVector('ys', traj_length)

right_start = z3.And(xs[0] == 0, ys[0] == 0)

isnt_jumping_vert = []
for t in range(traj_length - 1):
    isnt_jumping_vert.append(z3.Or([xs[t+1] == xs[t], xs[t+1] == xs[t] + 1,
                                    xs[t+1] == xs[t] - 1]))
no_jump_vert = z3.And(isnt_jumping_vert)

isnt_jumping_horiz = []
for t in range(traj_length - 1):
    isnt_jumping_horiz.append(z3.Or([ys[t+1] == ys[t], ys[t+1] == ys[t] + 1,
                                     ys[t+1] == ys[t] - 1]))
no_jump_horiz = z3.And(isnt_jumping_horiz)

def ever_below(time, y_val):
    # means that before time, y is at some point less than or equal to y_val
    is_below = []
    for t in range(time):
        is_below.append(ys[t] <= y_val)
    return z3.Or(is_below)

def eventually_up_right(time, x_val, y_val):
    # means that after time, x is always greater than or equal to -x_val,
    # and y is always greater than or equal to -y_val
    is_up_right = []
    for t in range(time, traj_length):
        is_up_right.append(z3.And(xs[t] >= -x_val, ys[t] >= -y_val))
    return z3.And(is_up_right)

def synthesise_me_a_trajectory(time1, y_val1, time2, x_val2, y_val2):
    S = z3.Solver()
    S.add(right_start)
    S.add(z3.And(no_jump_vert, no_jump_horiz))
    S.add(ever_below(time1, y_val1))
    S.add(eventually_up_right(time2, x_val2, y_val2))
    if S.check() == z3.sat:
        model = S.model()
        x_vals = [model.eval(xs[t]) for t in range(traj_length)]
        y_vals = [model.eval(ys[t]) for t in range(traj_length)]
        print("x coordinates:", x_vals)
        print("y coordinates:", y_vals)
    else:
        print("No trajectory can be synthesised with parameters",
              (time1, y_val1, time2, x_val2, y_val2))

def synthesise_trajectory_within_boxes(true_params, false_params):
    # each argument is a list of tuples (time1, y_val1, time2, x_val2, y_val2)
    # the ones in true_params are true, the ones in false_params are false.

    assert isinstance(true_params, list)
    assert isinstance(false_params, list)
    assert len(true_params) == len(false_params)
    S = z3.Solver()
    S.add(right_start)
    S.add(z3.And(no_jump_vert, no_jump_horiz))
    for tup in true_params:
        # somehow assert that tup is a tuple of integers
        assert(len(tup) == 5)
        S.add(z3.And(ever_below(tup[0], tup[1]),
                     eventually_up_right(tup[2], tup[3], tup[4])))
    for tup in false_params:
        # somehow assert that tup is a tuple of integers
        assert(len(tup) == 5)
        S.add(z3.Not(z3.And(ever_below(tup[0], tup[1]),
                            eventually_up_right(tup[2], tup[3], tup[4]))))
    if S.check() == z3.sat:
        model = S.model()
        x_vals = [model.eval(xs[t]) for t in range(traj_length)]
        y_vals = [model.eval(ys[t]) for t in range(traj_length)]
        print("x coordinates:", x_vals)
        print("y coordinates:", y_vals)
    else:
        print("No trajectory can be synthesised")
        print("List of true parameters:", true_params)
        print("List of false parameters:", false_params)

        
if __name__ == '__main__':
    print("Synthesised trajectory for single set of parameters:")
    synthesise_me_a_trajectory(4,-3,15,-6,-6)
    true_params = [(4,-3,15,-4,-4), (8,-4,18,-7,-7)]
    false_params = [(6,-5,16,-8,-8), (15,-7,14,-10,-10)]
    print("Synthesised trajectory for lists of true and false parameters:")
    synthesise_trajectory_within_boxes(true_params, false_params)
