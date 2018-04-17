import z3
import boundary_to_trace_v2

def label(boundary):
    # takes a boundary, and returns its proper label
    # implementation will depend on context, and in the extreme case will
    # require computation done by the human
    pass

def grid_point_to_params(grid_point, grid_dims):
    # takes a point in the discretisation of parameter space, returns a point
    # in parameter space
    # grid_point is a 2-tuple of natural numbers.
    # grid_dims is a dict containing variables below.
    # x_spacing is the width of bins of the first parameter.
    # y_spacing is analogous.
    x_spacing = grid_dims['x_spacing']
    y_spacing = grid_dims['y_spacing']
    return ((grid_point[0] + 0.5)*x_spacing, (grid_point[1] + 0.5)*y_spacing)
    
def grid_to_boundary(grid_boundary, grid_dims):
    # this function takes a boundary in the discretisation of parameter space,
    # and returns a boundary in parameter space
    # grid_boundary is a list of 2-tuples of natural numbers, representing a
    # boundary in a discretised version of parameter space.
    # grid_dims is a dict characterising the scale of discretisation of the grid
    # representing parameter space
    boundary = []
    for point in grid_boundary:
        boundary.append(grid_point_to_params(point, grid_dims))
    return boundary

def classify_n_kinks(n, grid_dims):
    # this function iterates through discretised boundaries with n 'kinks',
    # classifies them, and returns an artefact recording these classifications
    # TODO: decide if that artefact is a list or a dict or something
    # n is the integer number of 'kinks' that boundaries should have
    # grid_dims is a dict containing variables below.
    # x_spacing is the width of bins of the first parameter.
    # num_x is the number of bins of the first parameter.
    # y_spacing and num_y are analogous
    x_spacing = grid_dims['x_spacing']
    num_x = grid_dims['num_x']
    y_spacing = grid_dims['y_spacing']
    num_y = grid_dims['num_y']
    pass

def classify_grid_trajectories(grid_dims):
    # this function iterates through numbers of kinks, classifies all
    # discretised boundaries with that number of kinks, and returns an artefact
    # recording these classifications
    # grid_dims is a dict containing variables below.
    # x_spacing is the width of bins of the first parameter.
    # num_x is the number of bins of the first parameter.
    # y_spacing and num_y are analogous
    x_spacing = grid_dims['x_spacing']
    num_x = grid_dims['num_x']
    y_spacing = grid_dims['y_spacing']
    num_y = grid_dims['num_y']
    pass
