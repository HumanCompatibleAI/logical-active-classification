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
    endpoints -- a list of points in parameter space
    distance -- a float that can be positive or negative
    """
    assert isinstance(distance, float), "second argument of move_middle_out should be float"
    mpoint = midpoint(endpoints)
    length = np.sqrt((endpoints[0][0] - endpoints[1][0])**2
                     + (endpoints[0][1] - endpoints[1][1])**2)
    my_ends = np.array(endpoints)
    direction = np.array([-1,1]) * ((my_ends[1] - my_ends[0])[::-1])
    return (mpoint + (distance/length)*direction).tolist()

def find_upper_endpoint_bounds(positive_example, tolerance_a, tolerance_b,
                               param_boundaries):
    """
    Find upper bounds for where the endpoints of the boundary can be

    Returns a tuple of 4 coordinates: (left_upper_bound, right_upper_bound,
    left_lower_bound, right_lower_bound)
    positive_example should be some boundary.
    tolerance_a and tolerance_b should be floats.
    param_boundaries should be a list of lists giving lower and upper bounds for
    the possible parameter values. Specifically, it should be a list of the form
    [[lower_bound_x, upper_bound_x], [lower_bound_y, upper_bound_y]]
    """
    assert isinstance(positive_example, list), "first argument of find_endpoint_bounds should be a boundary"
    assert isinstance(tolerance_a, float), "second argument of find_endpoint_bounds should be a float"
    assert isinstance(tolerance_b, float), "third argument of find_endpoint_bounds should be a float"
    assert tolerance_a > 0, "second argument of find_endpoint_bounds should be positive"
    assert tolerance_b > 0, "third argument of find_endpoint_bounds should be positive"
    assert isinstance(param_boundaries, list) and len(param_boundaries) == 2 and isinstance(param_boundaries[0], list) and len(param_boundaries[0]) == 2 and isinstance(param_boundaries[1], list) and len(param_boundaries[1]) == 2, "fourth argument of find_endpoint_bounds should be list of starting and ending points in parameter space"
    have_found_top_right_end = False
    example_ends = [positive_example[0], positive_example[-1]]
    pos_ends = example_ends
    # first, move top left end of positive example to maximum y value, and check
    # if that makes it negative.
    top_left_end = [pos_ends[0][0], param_boundaries[1][1]]
    test_ends = [top_left_end, pos_ends[1]]
    test_top_boundary = endpoints_to_boundary(test_ends, tolerance_a)
    if (not label(test_top_boundary)):
        # if result is negative, then find the limit for the top left end
        # between our positive example and the result
        top_left_limit = find_endpoint_bound(tolerance_a, tolerance_b, 0,
                                             pos_ends, test_ends)
    else:
        # if not, move it as far right as the right end of the boundary and
        # check if it's negative now.
        pos_ends = test_ends
        top_left_end = [pos_ends[1][0], param_boundaries[1][1]]
        test_ends = [top_left_end, pos_ends[1]]
        test_top_boundary = endpoints_to_boundary(test_ends, tolerance_a)
        if (not label(test_top_boundary)):
            # if this new result is negative, find the limit for the top left
            # end between the old positive result and the new negative results
            top_left_limit = find_endpoint_bound(tolerance_a, tolerance_b, 0,
                                                 pos_ends, test_ends)
        else:
            # if the new result still isn't negative, find the limit of where
            # the right end can go
            # first, move the right end as far as possible to the right
            pos_ends = test_ends
            top_right_end = [param_boundaries[0][1], pos_ends[1][1]]
            test_ends = [pos_ends[0], top_right_end]
            test_top_boundary = endpoints_to_boundary(test_ends, tolerance_a)
            if (not label(test_top_boundary)):
                # if this makes the boundary negative, find the limit for the
                # top right end
                top_right_limit = find_endpoint_bound(tolerance_a, tolerance_b,
                                                      1, pos_ends, test_ends)
                have_found_top_right_end = True
                # then move the top left end that far right.
                pos_ends[1] = top_right_limit
                top_left_end = [top_right_limit[0], pos_ends[0][1]]
                test_ends = [top_left_end, pos_ends[1]]
                test_top_boundary = endpoints_to_boundary(test_ends,
                                                          tolerance_a)
                if (not label(test_top_boundary)):
                    # if that makes it negative, find the limit of where the top
                    # end can go
                    top_left_limit = find_endpoint_bound(tolerance_a,
                                                         tolerance_b, 0,
                                                         pos_ends, test_ends)
                else:
                    # if that doesn't make it go negative, then that's just the
                    # boundary
                    top_left_limit = top_left_end
            else:
                # if moving the right end as far as possible to the right didn't
                # make it a negative example, try moving it as far up as
                # possible
                pos_ends = test_ends
                top_right_end[1] = param_boundaries[1][1]
                test_ends = [pos_ends[0], top_right_end]
                test_top_boundary = endpoints_to_boundary(test_ends,
                                                          tolerance_a)
                if (not label(test_top_boundary)):
                    # if moving the right end as far up as possible makes it
                    # negative, find how far up it can go
                    top_right_limit = find_endpoint_bound(tolerance_a,
                                                          tolerance_b, 1,
                                                          pos_ends, test_ends)
                    have_found_top_right_end = True
                    # now, check if moving the left end as far as possible to
                    # the right makes it negative
                    top_left_end = [param_boundaries[0][1], pos_ends[0][1]]
                    test_ends = [top_left_end, top_right_limit]
                    test_top_boundary = endpoints_to_boundary(test_ends,
                                                              tolerance_a)
                    if (not label(test_top_boundary)):
                        # if that did make it negative, find the top left limit
                        # between those
                        top_left_limit = find_endpoint_bound(tolerance_a,
                                                             tolerance_b, 0,
                                                             pos_ends,
                                                             test_ends)
                    else:
                        # if moving the left end as far as possible to the right
                        # didn't make it negative, then that is now the top left
                        # limit
                        top_left_limit = top_left_end
                else:
                    # if moving the right end as far up as possible didn't make
                    # it negative, you've found the top right end.
                    top_right_limit = top_right_end
                    have_found_top_right_end = True
                    # you've also found the top left end
                    top_left_limit = pos_ends[0]

    # next, find the top right end if you haven't already
    if (not have_found_top_right_end):
        # to do this, first try moving the right end as far right as possible
        top_right_end = [param_boundaries[0][1], pos_ends[1][1]]
        test_ends = [pos_ends[0], top_right_end]
        test_top_boundary = endpoints_to_boundary(test_ends, tolerance_a)
        if (not label(test_top_boundary)):
            # if that makes the boundary a negative example, then find the top
            # right limit between pos_ends and test_ends
            top_right_limit = find_endpoint_bound(tolerance_a, tolerance_b, 1,
                                                  pos_ends, test_ends)
        else:
            # if that didn't make the boundary a negative example, move the
            # right end up to the height of the left end
            pos_ends = test_ends
            top_right_end = [pos_ends[1][0], pos_ends[0][1]]
            test_ends = [pos_ends[0], top_right_end]
            test_top_boundary = endpoints_to_boundary(test_ends, tolerance_a)
            if (not label(test_top_boundary)):
                # if that made the boundary a negative example, then find the
                # top right limit between pos_ends and test_ends
                top_right_limit = find_endpoint_bound(tolerance_a, tolerance_b,
                                                      1, pos_ends, test_ends)
            else:
                # if that didn't make the boundary a negative example, then that
                # must be the top right limit
                top_right_limit = top_right_end
    return (top_left_limit, top_right_limit)


def find_endpoint_bounds(positive_example, tolerance_a, tolerance_b,
                         param_boundaries):
    """
    Find upper and lower bounds for where the endpoints of the boundary can be

    Returns a tuple of 4 coordinates: (left_upper_bound, right_upper_bound,
    left_lower_bound, right_lower_bound)
    positive_example should be some boundary.
    tolerance_a and tolerance_b should be floats.
    param_boundaries should be a list of lists giving lower and upper bounds for
    the possible parameter values. Specifically, it should be a list of the form
    [[lower_bound_x, upper_bound_x], [lower_bound_y, upper_bound_y]]
    """
    # use helper function to find upper bounds
    upper_bounds = find_upper_endpoint_bounds(positive_example, tolerance_a,
                                              tolerance_b, param_boundaries)
    # negate everything, find upper bounds of the negated thing, then negate
    # that to find lower bounds
    neg_positive_example = [[(-1)*a for a in pair] for pair in positive_example]
    neg_param_boundaries = [[(-1)*a for a in pair] for pair in param_boundaries]
    neg_upper_bounds = find_upper_endpoint_bounds(neg_positive_example,
                                                  tolerance_a, tolerance_b,
                                                  neg_param_boundaries)
    return (upper_bounds[0], upper_bounds[1],
            neg_upper_bounds[0], neg_upper_bounds[1])
    
def find_endpoint_bound(tolerance_a, tolerance_b, vary_end, pos_ends, neg_ends):
    assert isinstance(tolerance_a, float), "first argument of find_endpoint_bound should be the tolerance to which you want to generate boundaries"
    assert tolerance_a > 0, "tolerance_a argument to find_endpoint_bound should be positive"
    assert isinstance(tolerance_b, float), "first argument of find_endpoint_bound should be the tolerance to which you want to find your endpoint bound"
    assert tolerance_b > 0, "tolerance_b argument to find_endpoint_bound should be positive"
    assert isinstance(vary_end, int), "vary_end should be an int in find_endpoint_bound"
    assert 0 <= vary_end and vary_end <= 1, "vary_end should be 0 or 1 in find_endpoint_bound"
    assert pos_ends[1 - vary_end] == neg_ends[1 - vary_end], "end you're not varying in find_endpoint_bound should be fixed between examples"
    # write more assertions later
    # write docstring later
    if pos_ends[vary_end][0] != neg_ends[vary_end][0]:
        vary_index = 0
    else:
        vary_index = 1
    distance = abs(pos_ends[vary_end][vary_index]
                   - neg_ends[vary_end][vary_index])
    num_iters = np.ceil(np.log2(distance / tolerance_b))
    num_iters = max(0, num_iters)
    for i in range(num_iters):
        test_ends = pos_ends
        test_val = 0.5*(pos_ends[vary_end][vary_index]
                        + neg_ends[vary_end][vary_index])
        test_ends[vary_end][vary_index] = test_val
        test_boundary = endpoints_to_boundary(test_ends, tolerance_a)
        if label(test_boundary):
            pos_ends[vary_end][vary_index] = test_val
        else:
            neg_ends[vary_end][vary_index] = test_val
    return pos_ends[vary_end]

def distance_point_to_line(point, a, b, c):
    """Distance from a point to the line ax + by + c = 0"""
    return abs(a*point[0] + b*point[1] + c) / np.sqrt(a**2 + b**2)

def maximally_extend_segment(endpoints, index, tolerance_a, tolerance_b,
                             is_positive):
    """
    Find the furthest out you can go between endpoints[index] and index+1.

    Arguments:
    endpoints -- list of endpoints of line segments connected together.
    index -- index of endpoints. We're going to return a point between 
             endpoints[index] and endpoints[index + 1].
    tolerance_a -- float, representing how much error we're OK with when 
                   generating boundaries from list of endpoints
    tolerance_b -- float, representing how much error we're OK with when 
                   approximating the convex set
    is_positive -- bool, representing whether we're extending a point out in the
                   positive direction or negative direction

    Returns:
    A point that is on the perpendicular bisector of the line segment between
    endpoints[index] and endpoints[index + 1], such that if that point were 
    further out and included in endpoints, the corresponding boundary would 
    barely be outside the convex set.
    """
    assert isinstance(endpoints, list), "first argument of maximally_extend_segment should be list of endpoints"
    assert len(endpoints) >= 2, "first argument of maximally_extend_segment should have at least 2 entries"
    assert index > 0, "second argument of maximally_extend_segment should be positive index"
    assert index < len(endpoints) - 1, "second argument of maximally_extend_segment should be index before end of list"
    assert isinstance(tolerance_a, float), "third argument of maximally_extend_segment should be a float representing how finely we're generating boundaries"
    assert isinstance(tolerance_b, float), "fourth argument of maximally_extend_segment should be a float representing how finely we're approximating the convex set"
    assert tolerance_a > 0, "tolerance_a should be positive in maximally_extend_segment"
    assert tolerance_b > 0, "tolerance_b should be positive in maximally_extend_segment"
    if is_positive:
        sign = 1
    else:
        sign = -1
    distance_min = 0
    mid = midpoint([endpoints[index], endpoints[index + 1]])
    distances = []
    if index > 0:
        # ensure we don't go beyond line connecting two previous endpoints
        a = endpoints[index][1] - endpoints[index - 1][1]
        b = endpoints[index - 1][0] - endpoints[index][0]
        c = (endpoints[index - 1][0] * endpoints[index][1]
             - endpoints[index][0] * endpoints[index - 1][1])
        distances.append(distance_point_to_line(mid, a, b, c))
    if index < len(endpoints) - 2:
        # ensure we don't go beyond line connecting two next endpoints
        a = endpoints[index + 2][1] - endpoints[index + 1][1]
        b = endpoints[index + 1][0] - endpoints[index + 2][0]
        c = (endpoints[index + 1][0] * endpoints[index + 2][1]
             - endpoints[index + 2][0] * endpoints[index + 1][1])
        distances.append(distance_point_to_line(mid, a, b, c))
    # ensure we don't go beyond endpoints in this interval
    length = np.sqrt((endpoints[index + 1][0] - endpoints[index][0])**2
                     + (endpoints[index + 1][1] - endpoints[index][1])**2)
    dist_to_top = ((-0.5)*length*(endpoints[index + 1][1] - endpoints[index][1])
                   / (endpoints[index + 1][0] - endpoints[index][0]))
    distances.append(dist_to_top)
    dist_to_right = ((-0.5)*length*(endpoints[index+1][0] - endpoints[index][0])
                     / (endpoints[index + 1][1] - endpoints[index][1]))
    distances.append(dist_to_right)
    distance_max = sign*min(distances)
    # find number of iterations we need
    num_iters = np.ceil(np.log2((distance_max - distance_min) / tolerance_b))
    num_iters = max(num_iters, 0)
    # do bisection on distance to move middle out
    for i in range(num_iters):
        test_distance = 0.5*(distance_max + distance_min)
        test_bump = move_middle_out([endpoints[index], endpoints[index + 1]],
                                    test_distance)
        test_endpoints = (endpoints[0:index+1] + [test_bump]
                          + endpoints[index + 1:])
        test_boundary = endpoints_to_boundary(test_endpoints, tolerance_a)
        if label(test_boundary):
            distance_min = test_distance
        else:
            distance_max = test_distance
    bump = move_middle_out([endpoints[index], endpoints[index + 1]],
                           distance_min)
    return bump

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

#b = ([[0, 100], [1, 75], [2, 50], [3, 25]], [[0, 75], [1, 50], [2, 25], [3, 0], [4, 0]])
#t = [[0.5, 75], [1.5, 50]]

def classify_trace(bounds, trace):
    # how to tell if a boundary is between two bounds?
    # for every point on the boundary, find the two points on the upper boundary
    # with x-coordinates that surround you, take the weighted average of their
    # y-coordinates, see if you're below that. same deal with lower boundary.

    # bounds is a tuple of lists of points. the lists have the same length.
    
    # ced works on this
    upper_bound, lower_bound = bounds[0], bounds[1]
    trace_length = len(trace)
    between = True

    def find_surround(x, bound):
        for i in range(len(bound)):
            if bound[i][0] == x:
                return i, i
            elif bound[i][0] < x and bound[i+1][0] > x:
                return i, i+1

    for i in range(trace_length):
        upper_index_1, upper_index_2 = find_surround(trace[i][0], upper_bound)
        upper_weighted_y = (upper_bound[upper_index_1][1] + upper_bound[upper_index_2][1]) / 2 * 1.0
        lower_index_1, lower_index_2 = find_surround(trace[i][0], lower_bound)
        lower_weighted_y = (lower_bound[lower_index_1][1] + lower_bound[lower_index_2][1]) / 2 * 1.0
        if trace[i][1] > upper_weighted_y or trace[i][1] < lower_weighted_y:
            between = False

    return between
