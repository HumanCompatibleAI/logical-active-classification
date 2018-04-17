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
    # spacings is an array containing the width of bins of each parameter.
    spacings = grid_dims['spacings']
    assert (len(grid_point) == len(spacings),
            "dimensions don't match in grid_point_to_params")
    params = []
    for i in range(len(grid_point)):
        params.append(spacings[i]*(grid_point[i] + 0.5))
    return params
    
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

def classify_grid_trajectories(grid_dims):
    # this function iterates through discretised boundaries, classifies them,
    # and returns an artefact recording these classifications.
    # TODO: what should this artefact be?
    # grid_dims is a dict containing variables below.
    # spacings is an array containing the width of bins of each parameter.
    # nums is an array containing the number of bins for each parameter.
    spacings = grid_dims['spacings']
    nums = grid_dims['nums']
    pass
