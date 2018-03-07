# D TODO: write pseudocode for this algorithm

import data # the unlabelled data. maybe we augment it to add labels in places
import hierarchical_clustering # a hierarchical clustering algorithm

# relevant data structures:
# Tree. A tree is like a normal tree, but each node corresponds to a subset of
# the unlabelled data. The root note corresponds to all the data, and the
# children of a node correspond to subsets of their parent node that are
# exhaustive and non-intersecting. We're going to have one tree, and it's going
# to be output by our hierarchical clustering algorithm. Tree has the usual
# structure with nodes etc.
# Pruning. A pruning is a selection of nodes in a tree that correspond to subsets
# that don't intersect and cover all of the data. This is like a cut in the tree
# between the root and the bottom.
# Partial pruning. A subset of a pruning.
# Labelling: dict associating nodes with labels.

labels = ['label1', 'label2'] # our list of labels

def empirical_prob_node_label(node, label, data):
    # return fraction of labelled points in node that have the right label
    pass

def empirical_error_node(node, data, labels):
    # return 1 - (max over labels of
    #             empirical_prob_node_label(node, label, data))
    pass

def conf_int_node_label(node, label, data, labels):
    # return tuple of two numbers that we have high confidence that the interval
    # contains the true probability of a data point coming from that node having
    # that label
    pass

beta = 2

def admissible_label(node, label, data, beta, labels, was_admissible):
    # returns a bool. true if we're sure that the label makes at most beta times
    # as much classification error on this node as any other label, or if we've
    # ever been sure of that. In original paper, beta = 2.
    label_prob = conf_int_node_label(node, label, data)[0]
    if was_admissible:
        return True
    else:
        my_bool = True
        for label_ in labels:
            other_label_prob = conf_int_node_label(node, label_, data)[1]
            if (label_ != label):
                my_bool = my_bool & label_prob > 2 * other_label_prob - 1
        return my_bool

def admissible_pruning_labelling(pruning, labelling, data):
    # returns a bool.
    return # labelling[node] defined for node in pruning and ancestors of pruning
           # AND for any node a strict ancestor of the pruning, labelling[node]
           # is admissible AND for any node in the pruning ( labelling[node] is
           # admissible OR no label is admissible and labelling[node] is
           # admissible for node's parent)

def weight_node(node, data):
    # return proportion of data contained in node
    pass

def weight_part_pruning(part_pruning, data):
    # unit test: for a pruning, weight_part_pruning(pruning, data) should be 1
    weight = 0
    for node in part_pruning:
        weight += weight_node(node)
    return weight

def empirical_error_part_pruning_labelling(part_pruning, labelling, data):
    weighted_error = 0
    for node in part_pruning:
        weighted_error += weight_node(node, data)
                          * empirical_error_node(node, data, labels)
    normaliser = weight_part_pruning(part_pruning, data)
    return weighted_error / normaliser
