import z3
import boundary_to_trace_v2
import numpy as np
import matplotlib.pyplot as plt

def label(boundary):
    """
    Takes a boundary, and returns its proper label, which is True or False.

    Correct implementation will depend on context, and in the extreme case will
    require computation done by the human.
    """
    pass

def get_positive_example():
    """Somehow samples a boundary that should be classified positively"""
    pass

def midpoint(endpoints):
    """Return the midpoint of a list of two endpoints"""
    assert isinstance(endpoints, list) or isinstance(endpoints, tuple), "endpoints of wrong type in function midpoint"
    assert len(endpoints) == 2, "you didn't give exactly two endpoints in the midpoint function"
    points = np.array(endpoints)
    return 0.5*(points[0] + points[1])

def endpoints_to_boundary(list_endpoints, tolerance):
    """
    Return a list of points 'tolerance' apart that lie between the endpoints

    Take a list of endpoints, representing a piecewise linear boundary, and
    return the 'actual boundary', i.e. a list of points, each of which is at 
    distance 'tolerance' (ideally a small number) from the next one.
    """
    assert isinstance(list_endpoints, list), "first argument of endpoints_to_boundary should be a list"
    assert isinstance(tolerance, float), "second argument of endpoints_to_boundary should be a float"
    assert tolerance > 0, "you gave a non-positive tolerance in endpoints_to_boundary"
    boundary = []
    for i in range(len(list_endpoints) - 1):
        start = list_endpoints[i]
        end = list_endpoints[i+1]
        total_dist = np.sqrt((end[0] - start[0])**2 + (end[1] - start[1])**2)
        n = int(np.floor(total_dist / tolerance))
        for j in range(n):
            frac = j * tolerance / total_dist
            waypoint = [frac*start[0] + (1-frac)*end[0],
                        frac*start[1] + (1-frac)*end[1]]
            boundary.append(waypoint)
    return boundary

# # test of endpoints_to_boundary
# my_list_endpoints = [[0,1], [0.7,0.5], [1,0]]
# my_tolerance = 0.01
# my_boundary = endpoints_to_boundary(my_list_endpoints, my_tolerance)
# # print(my_boundary)
# # print(zip(*my_boundary))
# plt.scatter(*zip(*my_boundary))
# plt.show()

def move_middle_out(endpoints, distance):
    """
    take a line segment, return one with the middle moved out perpendicularly

    this function takes a straight line segment, with endpoints given as an
    argument, and moves the middle out perpendicularly some distance. if the
    distance is positive, the middle should be pushed out into the top-right
    (i.e. the positive direction), and if distance is negative, the middle 
    should be pulled into the bottom-left (i.e. the negative direction).
    the return value is the final location of the midpoint.

    Argument types:
    endpoints -- a tuple of points in parameter space
    distance -- a float that can be positive or negative
    """
    mpoint = midpoint(endpoints)
    length = np.sqrt((endpoints[0][0] - endpoints[1][0])**2
                     + (endpoints[0][1] - endpoints[1][1])**2)
    my_ends = np.array(endpoints)
    direction = np.array([-1,1]) * ((my_ends[1] - my_ends[0])[::-1])
    return (mpoint + (distance/length)*direction).tolist()

def find_endpoint_bounds(positive_example, tolerance, param_boundaries):
    """
    Find upper+lower bounds for where the endpoints of the boundary can be

    Returns a tuple of 4 coordinates: (left_upper_bound, right_upper_bound,
    left_lower_bound, right_lower_bound)
    positive_example should be some boundary.
    tolerance should be a float.
    param_boundaries should be a list of lists giving lower and upper bounds for
    the possible parameter values.
    """
    # TODO: write type asserts
    example_ends = [positive_example[0], positive_example[-1]]
    top_end = [positive_example[0][0], param_boundaries[1][1]]
    test_top_boundary = endpoints_to_boundary([top_end, example_ends[1]])
    if label(test_top_boundary):
        something
    else:
        blah
    pass

def find_endpoint_bound(tolerance, vary_end, pos_ends, neg_ends):
    assert isinstance(tolerance, float), "first argument of find_endpoint_bound should be the tolerance to which you want to find your endpoint bound and generate boundaries"
    assert isinstance(vary_end, int), "vary_end should be an int in find_endpoint_bound"
    assert 0 <= vary_end and vary_end <= 1, "vary_end should be 0 or 1 in find_endpoint_bound"
    assert pos_ends[1 - vary_end] == neg_ends[1 - vary_end], "end you're not varying in find_endpoint_bound should be fixed between examples"
    # write more assertions later
    if pos_ends[vary_end][0] != neg_ends[vary_end][0]:
        vary_index = 0
    else:
        vary_index = 1
    distance = abs(pos_ends[vary_end][vary_index]
                   - neg_ends[vary_end][vary_index])
    num_iters = np.ceil((-1)*np.log2(tolerance / distance))
    for i in range(num_iters):
        test_ends = pos_ends
        test_val = 0.5*(pos_ends[vary_end][vary_index]
                        + neg_ends[vary_end][vary_index])
        test_ends[vary_end][vary_index] = test_val
        test_boundary = endpoints_to_boundary(test_ends, tolerance)
        if label(test_boundary):
            pos_ends[vary_end][vary_index] = test_val
        else:
            neg_ends[vary_end][vary_index] = test_val
    return pos_ends

def maximally_extend_segment(endpoints, index, tolerance_a, tolerance_b,
                             is_positive):
    pass

def interleave(small_list, big_list):
    """Interleave two lists of different length. First arg should be shorter."""
    assert isinstance(small_list, list), "First argument of interleave should be a list, is not."
    assert isinstance(big_list, list), "Second argument of interleave should be a list, is not."
    assert len(small_list) <= len(big_list), "Second argument of interleave should be longer list than first argument."
    new_list = []
    for i in range(len(small_list)):
        new_list.append(big_list[i])
        new_list.append(small_list[i])
    for j in range(len(big_list) - len(small_list)):
        new_list.append(big_list[j + len(small_list)])
    return new_list

def find_positive_set(iterations, tolerance_a, tolerance_b, param_boundaries):
    """
    Find the set of boundaries that should be classified positively.

    First, find the upper and lower edge bounds. Then, move the midpoints of 
    positively-classified boundaries out until they are negatively classified.
    Recurse for the midpoints of the line segments created in this process for
    some number of iterations.

    Argument types:
    iterations -- integer
    tolerance_a -- float. Should be positive.
    tolerance_b -- float. Should be positive.
    param_boundaries -- list of form [[lower_bound_x, upper_bound_x], 
                                      [lower_bound_y, upper_bound_y]],
                        describing the dimensions of parameter space

    note: should really separate out tolerance for how finely we're generating
    boundaries vs tolerance for how careful we are about where the boundary of 
    the convex set is. at the moment we'll just call them tolerance_a and 
    tolerance_b, need better names
    """
    assert isinstance(iterations, int) or isinstance(iterations, long), "iterations isn't an integer in find_positive_set"
    assert iterations > 0, "you can't have a non-positive number of iterations in find_positive_set"
    assert isinstance(tolerance_a, float), "tolerance_a should be a float in find_positive_set"
    assert tolerance_a > 0, "tolerance_a should be positive in find_positive_set"
    assert isinstance(tolerance_b, float), "tolerance_b should be a float in find_positive_set"
    assert tolerance_b > 0, "tolerance_b should be positive in find_positive_set"
    
    positive_example = get_positive_example()
    endpoint_bounds = find_endpoint_bounds(positive_example, tolerance_b,
                                           param_boundaries)

    upper_bound = [endpoint_bounds[0], endpoint_bounds[1]]
    lower_bound = [endpoint_bounds[2], endpoint_bounds[3]]

    for i in range(iterations):
        assert len(upper_bound) == len(lower_bound), "somehow upper and lower bounds became a different length in the loop of find_positive_set"
        upper_extensions = []
        lower_extensions = []
        for j in range(len(upper_bound) - 1):
            new_upper = maximally_extend_segment(upper_bound, j, tolerance_a,
                                                 tolerance_b, True)
            new_lower = maximally_extend_segment(lower_bound, j, tolerance_a,
                                                 tolerance_b, False)
            upper_extensions.append(new_upper)
            lower_extensions.append(new_lower)
        assert len(upper_extensions) == len(lower_extensions), "upper_extensions and lower_extensions should end up having the same length after the end of the inner loop of find_positive_set"
        assert len(upper_extensions) == len(upper_bound) - 1, "upper_extensions should have length 1 less than upper_bound after the inner loop of find_positive_set"
        upper_bound = interleave(upper_extensions, upper_bound)
        lower_bound = interleave(lower_extensions, lower_bound)
    return (upper_bound, lower_bound)

def classify_trace(bounds, trace):
    # how to tell if a boundary is between two bounds?
    # for every point on the boundary, find the two points on the upper boundary
    # with x-coordinates that surround you, take the weighted average of their
    # y-coordinates, see if you're below that. same deal with lower boundary.
    pass
