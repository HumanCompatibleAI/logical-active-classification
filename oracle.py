import boundary_to_trace_v4 as bt
import trace_to_boundary as tb
import z3
import utils

def label(boundary, should_invert=False, param_boundaries=[[0.0, 20.0],
                                                           [0,1.0]],
          epsilon=0.01, num_points=200, lipschitz_param=0.05):
    """
    Takes a boundary, and returns its proper label, which is True or False.

    Correct implementation will depend on context, and in the extreme case will
    require computation done by the human.
    """
    # TODO: add type assertions
    if should_invert:
        my_boundary = utils.invert(boundary, param_boundaries)
    else:
        my_boundary = boundary
    utils.plot(my_boundary)
    print("in the label function")
    
    # print("in the label function, should_invert is", should_invert)
    # print("boundary that label is synthesising", my_boundary)
    # print(num_points)
    
    xs = [z3.Real('x%d' % i) for i in range(num_points)]
    us = [z3.Real('u%d' % i) for i in range(num_points)]

    def make_phi(xs, us):
        # this basically makes a formula encoding dynamics of a system xs
        # controlled by us
        formula = True
        for i in range(num_points - 1):
            formula = z3.And(formula, xs[i+1] == xs[i] + us[i],
                             # xs[i+1] - xs[i] <= lipschitz_param,
                             # xs[i] - xs[i+1] <= lipschitz_param,
                             xs[i] >= 0, xs[i+1] >= 0, xs[i] <= 1, xs[i+1] <= 1)
        return formula

    trace = bt.trace(my_boundary, epsilon, num_points, xs,
                     param_boundaries[0][1], make_phi(xs, us))

    # print("trace that is being labelled is", trace)
    # utils.plot(trace)
    class_val = True
    
    for point in trace:
        if point[0] >= 10:
            class_val = class_val and (point[1] <= (-0.08)*point[0] + 1.8)

    return class_val

def get_positive_example(param_boundaries):
    """Samples a boundary that should be classified positively.

    This will be specific to one particular hypothesis.
    """
    print("We are now inside get_positive_example")

    pos_trace = []
    for i in range(201):
        x_val = ((param_boundaries[0][1] - param_boundaries[0][0])*(i/200.0)
                 + param_boundaries[0][0])
        m = (-0.8)/19.8
        b = 0.9 - 0.1*m
        y_val = m*x_val + b
        pos_trace.append([x_val, y_val])
    
    boundary = tb.boundary(pos_trace)
    nice_boundary = [point for point in boundary
                     if (point[0] >= param_boundaries[0][0]
                         and point[0] <= param_boundaries[0][1]
                         and point[1] >= param_boundaries[1][0]
                         and point[1] <= param_boundaries[1][1])]

    # true_label = True
    # for point in pos_trace:
    #     if point[0] >= 10:
    #         true_label = true_label and (point[1] <= (-0.08)*point[0] + 1.8)

    # print("positive trace we generate is", pos_trace)
            
    # print("true classification of your trace is", true_label)

    # print("boundary we generate is", nice_boundary)
    
    assert label(nice_boundary), "the positive example you generated isn't actually labelled positive!"
    return nice_boundary

# my_bound = get_positive_example([[0, 20], [0,1]])
# print(my_bound)
