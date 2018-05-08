import numpy as np
import matplotlib.pyplot as plt

def midpoint(endpoints):
    """Return the midpoint of a list of two endpoints"""
    assert isinstance(endpoints, list) or isinstance(endpoints, tuple), "endpoints of wrong type in function midpoint"
    assert len(endpoints) == 2, "you didn't give exactly two endpoints in the midpoint function"
    points = np.array(endpoints)
    return (0.5*(points[0] + points[1])).tolist()

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
        if total_dist == 0:
            boundary.append(start)
        n = int(np.ceil(total_dist / tolerance))
        for j in range(n):
            frac = j * tolerance / total_dist
            waypoint = [(1 - frac)*start[0] + frac*end[0],
                        (1 - frac)*start[1] + frac*end[1]]
            boundary.append(waypoint)
    return boundary

# print(endpoints_to_boundary([[0.0, 0.0], [0.0, 0.0]], 0.1))

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
    # this is just a formula that you can figure out if you really want to
    mpoint = midpoint(endpoints)
    length = np.sqrt((endpoints[0][0] - endpoints[1][0])**2
                     + (endpoints[0][1] - endpoints[1][1])**2)
    my_ends = np.array(endpoints)
    direction = np.array([-1,1]) * ((my_ends[1] - my_ends[0])[::-1])
    return (mpoint + (distance/length)*direction).tolist()

# print(move_middle_out([[0,1],[20,0]], 2.0))

def angle_line(point0, point1):
    """Takes two points, returns the angle of the line they define"""
    assert isinstance(point0, list), "first argument of angle_line should be list"
    assert len(point0) == 2, "first argument of angle_line should have length 2"
    assert isinstance(point1, list), "second argument of angle_line should be list"
    assert len(point1) == 2, "second argument of angle_line should have length 2"
    assert isinstance(point0[0], float) and isinstance(point0[1], float), "first argument of angle_line should be list of floats"
    assert isinstance(point1[0], float) and isinstance(point1[1], float), "second argument of angle_line should be list of floats"
    return np.arctan((point1[1] - point0[1]) / (point1[0] - point0[0]))

def distance_to_hit_line(point0, point1, point2):
    """
    Returns the distance that point0 can move perpendicularly to point 1 before
    hitting the line between point1 and point2.
    """
    len_0_1 = np.sqrt((point0[0] - point1[0])**2 + (point0[1] - point1[1])**2)
    t01 = angle_line(point0, point1)
    t12 = angle_line(point1, point2)
    t = abs(t12 - t01)
    return len_0_1*np.tan(t)

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

def invert(boundary, param_boundaries):
    """
    Takes a boundary, and rotates it by 180 degrees in parameter space

    Arguments:
    boundary -- a boundary
    param_boundaries -- list of lists giving upper and lower bounds for possible
                        parameter values. Specifically, should be of the form
                        [[lower_bound_x, upper_bound_x], 
                         [lower_bound_y, upper_bound_y]]

    Note that a 180 degree rotation is the same as flipping vertically and then
    horizontally, which is what we will actually do
    """
    new_boundary = []
    for point in boundary:
        new_x = param_boundaries[0][0] + param_boundaries[0][1] - point[0]
        new_y = param_boundaries[1][0] + param_boundaries[1][1] - point[1]
        new_boundary.append([new_x, new_y])
    return new_boundary[::-1]

def plot(data):
    """
    Convert a trace or a boundary to two lists, x and y axis. And plot it.
    """

    x_axis = [data[i][0] for i in range(len(data))]
    y_axis = [data[i][1] for i in range(len(data))]
    plt.plot(x_axis, y_axis)
    plt.show()
