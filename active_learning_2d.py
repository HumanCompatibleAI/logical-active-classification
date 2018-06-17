import boundary_to_trace_v4 as bt
import trace_to_boundary as tb
import numpy as np
import copy
import z3
import utils
import math
# from oracle_positive_polygon import label, get_positive_example
from oracle_polygon import label, get_positive_example
# from oracle import label, get_positive_example

def find_upper_endpoint_bounds(positive_example, tolerance_a, tolerance_b,
                               param_boundaries, should_invert):
    """
    Find bounds for where the endpoints of the boundary can be

    Returns a tuple of 4 coordinates: (left_upper_bound, right_upper_bound,
    left_lower_bound, right_lower_bound)

    positive_example should be some boundary.
    tolerance_a and tolerance_b should be floats.
    param_boundaries should be a list of lists giving lower and upper bounds for
    the possible parameter values. Specifically, it should be a list of the form
    [[lower_bound_x, upper_bound_x], [lower_bound_y, upper_bound_y]]
    should_invert should be a bool telling you if you want to invert the 
    boundary before labelling it. Ideally, it is true when we are secretly 
    actually finding lower bounds.
    """
    assert isinstance(positive_example, list), "first argument of find_upper_endpoint_bounds should be a boundary"
    assert isinstance(tolerance_a, float), "second argument of find_upper_endpoint_bounds should be a float"
    assert isinstance(tolerance_b, float), "third argument of find_upper_endpoint_bounds should be a float"
    assert tolerance_a > 0, "second argument of find_upper_endpoint_bounds should be positive"
    assert tolerance_b > 0, "third argument of find_upper_endpoint_bounds should be positive"
    assert isinstance(param_boundaries, list) and len(param_boundaries) == 2 and isinstance(param_boundaries[0], list) and len(param_boundaries[0]) == 2 and isinstance(param_boundaries[1], list) and len(param_boundaries[1]) == 2, "fourth argument of find_upper_endpoint_bounds should be list of starting and ending points in parameter space"
    assert isinstance(should_invert, bool), "fifth argument of find_upper_endpoint_bounds should be bool denoting whether or not boundaries should be inverted in the labelling function"

    print("now inside find_upper_endpoint_bounds")
    print("should_invert is", should_invert)
    
    have_found_top_right_end = False
    
    example_ends = [positive_example[0], positive_example[-1]]

    print("positive ends to work with", example_ends)
    
    pos_ends = example_ends

    # first, move top left end of positive example to maximum y value, and
    # check if that makes it negative.
    top_left_end = [pos_ends[0][0], param_boundaries[1][1]]
    test_ends = [top_left_end, pos_ends[1]]
    print("first boundary tried in find_upper_endpoint_bounds", test_ends)
    test_top_boundary = utils.endpoints_to_boundary(test_ends, tolerance_a)
    # print("actual full boundary", test_top_boundary)
    # print("classification of that boundary (maybe inverted)",
    #       label(test_top_boundary, should_invert, param_boundaries))
    if (not label(test_top_boundary, should_invert, param_boundaries)):
        # if result is negative, then find the limit for the top left end
        # between our positive example and the result
        top_left_limit = find_endpoint_bound(tolerance_a, tolerance_b, 0,
                                             pos_ends, test_ends,
                                             should_invert, param_boundaries)
    else:
        # if not, move it as far right as the right end of the boundary and
        # check if it's negative now.
        pos_ends = test_ends
        top_left_end = [pos_ends[1][0], param_boundaries[1][1]]
        test_ends = [top_left_end, pos_ends[1]]
        print("next boundary ends tried in find_upper_endpoint_bounds",
              test_ends)
        test_top_boundary = utils.endpoints_to_boundary(test_ends, tolerance_a)
        # print("next full boundary tried in find_upper_endpoint_bounds",
        #       test_top_boundary)
        if (not label(test_top_boundary, should_invert, param_boundaries)):
            # if this new result is negative, find the limit for the top left
            # end between the old positive result and the new negative results
            top_left_limit = find_endpoint_bound(tolerance_a, tolerance_b, 0,
                                                 pos_ends, test_ends,
                                                 should_invert,
                                                 param_boundaries)
        else:
            # if the new result still isn't negative, find the limit of where
            # the right end can go
            # first, move the right end as far as possible to the right
            pos_ends = test_ends
            top_right_end = [param_boundaries[0][1], pos_ends[1][1]]
            test_ends = [pos_ends[0], top_right_end]
            print("next boundary ends tried in find_endpoint_bounds", test_ends)
            test_top_boundary = utils.endpoints_to_boundary(test_ends, tolerance_a)
            if (not label(test_top_boundary, should_invert, param_boundaries)):
                # if this makes the boundary negative, find the limit for the
                # top right end
                top_right_limit = find_endpoint_bound(tolerance_a, tolerance_b,
                                                      1, pos_ends, test_ends,
                                                      should_invert,
                                                      param_boundaries)
                have_found_top_right_end = True
                # then move the top left end that far right.
                pos_ends[1] = top_right_limit
                top_left_end = [top_right_limit[0], pos_ends[0][1]]
                test_ends = [top_left_end, pos_ends[1]]
                print("next boundary ends tried in find_endpoint_bounds",
                      test_ends)
                test_top_boundary = utils.endpoints_to_boundary(test_ends,
                                                                tolerance_a)
                if (not label(test_top_boundary, should_invert,
                              param_boundaries)):
                    # if that makes it negative, find the limit of where the top
                    # end can go
                    top_left_limit = find_endpoint_bound(tolerance_a,
                                                         tolerance_b, 0,
                                                         pos_ends, test_ends,
                                                         should_invert,
                                                         param_boundaries)
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
                print("next boundary ends tried in find_endpoint_bounds",
                      test_ends)
                test_top_boundary = utils.endpoints_to_boundary(test_ends,
                                                                tolerance_a)
                if (not label(test_top_boundary, should_invert,
                              param_boundaries)):
                    # if moving the right end as far up as possible makes it
                    # negative, find how far up it can go
                    top_right_limit = find_endpoint_bound(tolerance_a,
                                                          tolerance_b, 1,
                                                          pos_ends, test_ends,
                                                          should_invert,
                                                          param_boundaries)
                    have_found_top_right_end = True
                    # now, check if moving the left end as far as possible to
                    # the right makes it negative
                    top_left_end = [param_boundaries[0][1], pos_ends[0][1]]
                    test_ends = [top_left_end, top_right_limit]
                    print("next boundary ends tried in find_endpoint_bounds",
                          test_ends)                    
                    test_top_boundary = utils.endpoints_to_boundary(test_ends,
                                                                    tolerance_a)
                    if (not label(test_top_boundary, should_invert,
                                  param_boundaries)):
                        # if that did make it negative, find the top left limit
                        # between those
                        top_left_limit = find_endpoint_bound(tolerance_a,
                                                             tolerance_b, 0,
                                                             pos_ends,
                                                             test_ends,
                                                             should_invert,
                                                             param_boundaries)
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
        print("finding top right end")
        # to do this, first try moving the right end as far right as possible
        top_right_end = [param_boundaries[0][1], pos_ends[1][1]]
        test_ends = [pos_ends[0], top_right_end]
        print("next boundary ends tried in find_endpoint_bounds", test_ends)
        test_top_boundary = utils.endpoints_to_boundary(test_ends, tolerance_a)
        if (not label(test_top_boundary, should_invert, param_boundaries)):
            # if that makes the boundary a negative example, then find the top
            # right limit between pos_ends and test_ends
            top_right_limit = find_endpoint_bound(tolerance_a, tolerance_b, 1,
                                                  pos_ends, test_ends,
                                                  should_invert,
                                                  param_boundaries)
        else:
            # if that didn't make the boundary a negative example, move the
            # right end up to the height of the left end
            pos_ends = test_ends
            top_right_end = [pos_ends[1][0], pos_ends[0][1]]
            test_ends = [pos_ends[0], top_right_end]
            print("next boundary ends tried in find_endpoint_bounds", test_ends)
            test_top_boundary = utils.endpoints_to_boundary(test_ends,
                                                            tolerance_a)
            if (not label(test_top_boundary, should_invert, param_boundaries)):
                # if that made the boundary a negative example, then find the
                # top right limit between pos_ends and test_ends
                top_right_limit = find_endpoint_bound(tolerance_a, tolerance_b,
                                                      1, pos_ends, test_ends,
                                                      should_invert,
                                                      param_boundaries)
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
    print("In find_endpoint_bounds.")
    
    # use helper function to find upper bounds
    upper_bounds = find_upper_endpoint_bounds(positive_example, tolerance_a,
                                              tolerance_b, param_boundaries,
                                              False)
    # negate the example, find upper bounds of the negated thing, then negate
    # that to find lower bounds
    # but instead of negating, you want to take the limits minus the example,
    # so that you remain within the same parameter boundaries
    inv_pos_example = utils.invert(positive_example, param_boundaries)
    # print("in find_endpoint_bounds, the positive example is",
    #        positive_example)
    # print("in find_endpoint_bounds, inv_pos_example is", inv_pos_example)
    inv_upper_bounds = find_upper_endpoint_bounds(inv_pos_example,
                                                  tolerance_a, tolerance_b,
                                                  param_boundaries, True)
    lower_bounds = utils.invert(inv_upper_bounds, param_boundaries)
    return (upper_bounds[0], upper_bounds[1],
            lower_bounds[0], lower_bounds[1])
    
def find_endpoint_bound(tolerance_a, tolerance_b, vary_end, pos_ends, neg_ends,
                        should_invert, param_boundaries):
    """
    Find the boundary of the convex set between a positive and negative example.

    Arguments:
    tolerance_a -- the tolerance to which we want to generate boundaries from
                   endpoints
    tolerance_b -- the tolerance to which we want to find the edges of the 
                   convex set
    vary_end -- the end that we want to vary to find the edge
    pos_ends -- the ends of a positively classified boundary
    neg_ends -- the ends of a negatively classified boundary
    should_invert -- bool that encodes whether we should invert the boundary 
                     before finding its label
    param_boundaries -- list of lists that encodes the borders of parameter 
                        space
    """
    assert isinstance(tolerance_a, float), "first argument of find_endpoint_bound should be the tolerance to which you want to generate boundaries"
    assert tolerance_a > 0, "tolerance_a argument to find_endpoint_bound should be positive"
    assert isinstance(tolerance_b, float), "first argument of find_endpoint_bound should be the tolerance to which you want to find your endpoint bound"
    assert tolerance_b > 0, "tolerance_b argument to find_endpoint_bound should be positive"
    assert isinstance(vary_end, int), "vary_end should be an int in find_endpoint_bound"
    assert 0 <= vary_end and vary_end <= 1, "vary_end should be 0 or 1 in find_endpoint_bound"
    assert isinstance(pos_ends, list), "fourth argument of find_endpoint_bound should be a list of positive ends"
    assert len(pos_ends) == 2, "fourth argument of find_endpoint_bound should have length 2, being two ends of a boundary"
    assert isinstance(pos_ends[0], list) and isinstance(pos_ends[1], list), "fourth argument of find_endpoint_bound should be a list of two points"
    assert len(pos_ends[0]) == 2 and len(pos_ends[1]) == 2, "both elements of pos_ends should be of length 2 in find_endpoint_bound"
    assert isinstance(neg_ends, list), "fourth argument of find_endpoint_bound should be a list of negative ends"
    assert len(neg_ends) == 2, "fourth argument of find_endpoint_bound should have length 2, being two ends of a boundary"
    assert isinstance(neg_ends[0], list) and isinstance(neg_ends[1], list), "fourth argument of find_endpoint_bound should be a list of two points"
    assert len(neg_ends[0]) == 2 and len(neg_ends[1]) == 2, "both elements of neg_ends should be of length 2 in find_endpoint_bound"
    assert pos_ends[1 - vary_end] == neg_ends[1 - vary_end], "end you're not varying in find_endpoint_bound should be fixed between examples"
    assert (pos_ends[vary_end][0] == neg_ends[vary_end][0]) or (pos_ends[vary_end][1] == neg_ends[vary_end][1]), "in find_endpoint_bound, the ends to vary between the positive and negative example should only be different in one dimension."
    # TODO: add type asserts for should_invert and param_boundaries

    print("Inside find_endpoint_bound")
    
    # figure out which dimension to vary
    if pos_ends[vary_end][0] != neg_ends[vary_end][0]:
        vary_index = 0
    else:
        vary_index = 1
    # figure out how many iterations we need to do
    distance = abs(pos_ends[vary_end][vary_index]
                   - neg_ends[vary_end][vary_index])
    num_iters = np.ceil(np.log2(distance / tolerance_b))
    print("num_iters before max-ing with 0", num_iters)
    num_iters = int(max(0, num_iters))
    # repeatedly go half-way between the two ends, see if that's positive or
    # negative, update accordingly
    for i in range(num_iters):
        print("have just looped in find_endpoint_bound")
        print("pos_ends in loop:", pos_ends)
        print("neg_ends in loop:", neg_ends)
        test_ends = copy.deepcopy(pos_ends)
        test_val = 0.5*(pos_ends[vary_end][vary_index]
                        + neg_ends[vary_end][vary_index])
        test_ends[vary_end][vary_index] = test_val
        # print("pos_ends after we change test_ends:", pos_ends)
        print("boundary ends we're testing in loop",
              test_ends)
        test_boundary = utils.endpoints_to_boundary(test_ends, tolerance_a)
        if label(test_boundary, should_invert, param_boundaries):
            pos_ends[vary_end][vary_index] = test_val
        else:
            neg_ends[vary_end][vary_index] = test_val
    return pos_ends[vary_end]

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
    assert index >= 0, "second argument of maximally_extend_segment should be non-negative index, but it's actually %d" % index
    assert index <= len(endpoints) - 2, "second argument of maximally_extend_segment should be index before end of list"
    assert isinstance(tolerance_a, float), "third argument of maximally_extend_segment should be a float representing how finely we're generating boundaries"
    assert isinstance(tolerance_b, float), "fourth argument of maximally_extend_segment should be a float representing how finely we're approximating the convex set"
    assert tolerance_a > 0, "tolerance_a should be positive in maximally_extend_segment"
    assert tolerance_b > 0, "tolerance_b should be positive in maximally_extend_segment"

    print("inside maximally_extend_segment")
    
    if is_positive:
        sign = 1
    else:
        sign = -1
    distance_min = 0
    # print("index is", index)
    # print("first endpoint", endpoints[index])
    # print("second endpoint", endpoints[index+1])
    mid = utils.midpoint([endpoints[index], endpoints[index + 1]])
    # print("midpoint", mid)
    len_mid_end = np.sqrt((mid[1] - endpoints[index+1][1])**2
                          + (mid[0] - endpoints[index+1][0])**2)
    if len_mid_end == 0:
        return mid
    distances = []
    if index > 0:
        # ensure we don't go beyond line connecting two previous endpoints
        dist = utils.distance_to_hit_line(mid, endpoints[index],
                                          endpoints[index - 1])
        distances.append(dist)
        # print("distance to line connecting previous two endpoints:", dist)
    if index < len(endpoints) - 2:
        # ensure we don't go beyond line connecting two next endpoints
        dist = utils.distance_to_hit_line(mid, endpoints[index + 1],
                                          endpoints[index + 2])
        distances.append(dist)
        # print("distance to line connecting two next endpoints:", dist)
    # ensure we don't go beyond endpoints in this interval
    if endpoints[index+1][0] == endpoints[index][0]:
        dist_to_top = math.inf
    else:
        dist_to_top = ((-1.0)*len_mid_end*(endpoints[index + 1][1]
                                           - endpoints[index][1])
                       / (endpoints[index + 1][0] - endpoints[index][0]))
    # print("distance to top endpoint:", dist_to_top)
    distances.append(dist_to_top)
    if endpoints[index+1][1] == endpoints[index][1]:
        dist_to_right = math.inf
    else:
        dist_to_right = ((-1.0)*len_mid_end*(endpoints[index+1][0]
                                             - endpoints[index][0])
                         / (endpoints[index + 1][1] - endpoints[index][1]))
    distances.append(dist_to_right)
    # print("distance to right endpoint:", dist_to_right)
    distance_max = sign*min(distances)
    # find number of iterations we need
    # print("distance_max", distance_max)
    # print("distance_min", distance_min)
    # print("tolerance_b", tolerance_b)
    num_iters = np.ceil(np.log2(abs(distance_max - distance_min) / tolerance_b))
    # print("num_iters before maxing with 0", num_iters)
    num_iters = int(max(num_iters, 0))
    # do bisection on distance to move middle out
    for i in range(num_iters):
        test_distance = 0.5*(distance_max + distance_min)
        test_bump = utils.move_middle_out([endpoints[index],
                                           endpoints[index + 1]],
                                    test_distance)
        test_endpoints = (endpoints[0:index+1] + [test_bump]
                          + endpoints[index + 1:])
        test_boundary = utils.endpoints_to_boundary(test_endpoints,
                                                    tolerance_a)
        if label(test_boundary):
            distance_min = test_distance
        else:
            distance_max = test_distance
    distance_min = float(distance_min)
    bump = utils.move_middle_out([endpoints[index], endpoints[index + 1]],
                                 distance_min)
    print("bumped point is", bump)
    return bump

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
    
    # get a positive example
    positive_example = get_positive_example(param_boundaries)

    print("Positive example gained.")

    # find where the endpoints can be in the convex set
    endpoint_bounds = find_endpoint_bounds(positive_example, tolerance_a,
                                           tolerance_b, param_boundaries)
    upper_bound = [endpoint_bounds[0], endpoint_bounds[1]]
    lower_bound = [endpoint_bounds[2], endpoint_bounds[3]]

    print("after finding endpoint bounds, upper bound is", upper_bound)
    print("after finding endpoint bounds, lower bound is", lower_bound)
    
    # repeatedly find the midpoints of the line segments in the lower and upper
    # bounds, and move them out as far as possible
    for i in range(iterations):
        print("in loop number " + str(i) + " when moving out midpoints")
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
        upper_bound = utils.interleave(upper_extensions, upper_bound)
        lower_bound = utils.interleave(lower_extensions, lower_bound)
    utils.plot(upper_bound, 'tab:orange', param_boundaries)
    utils.plot(lower_bound, 'tab:orange', param_boundaries)
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

#boundary = [[19.0, 0.9], [19.0, 0.7000000000000001], [19.0, 0.5], [19.0, 0.29999999999999993]]
#label(boundary)
bounds = find_positive_set(4, 0.05, 0.00625, [[0.0, 20.0], [0.0, 1.0]])
print(bounds)
