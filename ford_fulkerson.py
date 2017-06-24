#!/usr/bin/env python3

"Implementation of the greedy Ford-Fukerson algorithm for network flows"

__author__ = "João Francisco B. S. Martins"
__email__ = "joaofbsm@dcc.ufmg.br"
__license__ = "GPL"
__version__ = "3.0"

import sys
import numpy as np

def create_adjacency_matrix(c, N):
    """Creates matrix with weights based in arch capacity and arch connections.
    
    Arguments:
        c -- Arch capacity vector.
        N -- Incidence matrix.
    """

    vertex_num = N.shape[0]
    # Weighted matrix
    A = np.zeros((vertex_num, vertex_num))

    for i, vertex in enumerate(N):
        for j, edge in enumerate(vertex):
            if edge == -1:
                # If edge points outward the vertex, add its weight
                for k in range(vertex_num):
                    if N[k][j] == 1:
                        # k is vertex to where the edge points
                        A[i][k] = c[j]
                        break

    return A


def find_path(s, t, A, N):
    """Finds path using recursive depth first search.

    Returns path if t is reachable and set of reachable vertices otherwise. The
    path is the sequence of the indexes from visited vertices.

    Arguments:
        s -- Source vertex index.
        t -- Target(Sink) vertex index.
        A -- Weighted adjacency matrix.
        N -- Incidence matrix.
    """

    # Visited vertex
    visited = np.zeros(N.shape[0], dtype="int")
    # 
    stack = []
    path = np.zeros(N.shape[0], dtype="int")

    # Begin by visiting source
    stack.append((s, [s]))

    while stack:
        # Gets last element from the stack
        u, path = stack.pop()
        visited[u] = 1
        for v, capacity in enumerate(A[u]):
            if visited[v] == 0 and capacity > 0:
                if v == t:
                    # If found sink, return path to it
                    return path + [v], None
                else:
                    # Add connected vertex to the stack
                    stack.append((v, path + [v]))

    # If t is unreachable, return set of reachable vertices
    return None, visited


def find_edge(flow, A, N, u, v):
    """Find edge number for the edge between u and v.

    Returns the edge index and its directions: 1 if from u to v and -1 
    otherwise.
    
    Arguments:
        flow -- Array of flows per edge.
        A -- Weighted adjacency matrix.
        N -- Incidence matrix.
        u -- Vertex index for the first end of the edge.
        v -- Vertex index for the second end of the edge.
    """

    for edge, incidence in enumerate(N[u]):
        if incidence == -1 and N[v][edge] == 1 and A[u][v] > 0:
            return edge, N[v][edge]
        elif incidence == 1 and N[v][edge] == -1 and flow[edge] > 0:
            return edge, N[v][edge]


def convert_path(path, flow, A, N):
    """Converts path from vertex list to list of edge presence on it.
    
    New path has 1 in edge index if edge is present and 0 if not.
    
    Arguments:
        flow -- Array of flows per edge.
        path -- Path as returned from find_path().
        A -- Weighted adjacency matrix.
        N -- Incidence matrix.
    """

    new_path = np.zeros(N.shape[1], dtype="int")
    for i in range(len(path) - 1):
        u = path[i]
        v = path[i + 1]
        edge, direction = find_edge(flow, A, N, u, v)
        new_path[edge] = direction

    return new_path


def ind_to_adj(N, A, edge, direction):
    """Finds edge in incidence matrix and return adjacency coordinates for it.
    
    Arguments:
        N -- Incidence matrix.
        A -- Adjacency matrix.
        edge -- Edge index.
        direction -- Edge direction: 1 if from u to v and -1 otherwise.
    """

    u = -1
    v = -1
    for i in range(N.shape[0]):
        if N[i][edge] == direction * -1:
            u = i
        elif N[i][edge] == direction:
            v = i
    return u, v
                

def find_min_flow(path, A, N):
    """Returns the minimum flow in current path to t.

    Arguments:
        path -- Path as returned from convert_path().
        A -- Weighted adjacency matrix.
        N -- Incidence matrix.
    """

    min_flow = float("inf")
    for edge, direction in enumerate(path):
        u, v = ind_to_adj(N, A, edge, direction)
        if A[u][v] != 0:
            flow = A[u][v]
            min_flow = min(min_flow, flow)
    return int(min_flow)


def update_graph(N, A, flow, path, min_flow):
    """Updates flows and capacities on the graph.
    
    Arguments:
        N -- Incidence matrix.
        A -- Weighted adjacency matrix.
        flow -- Array of flows per edge.
        path -- Path as returned from convert_path().
        min_flow -- Minimum flow in the path.
    """

    for edge, direction in enumerate(path):
        if direction != 0:
            flow[edge] += min_flow * direction
            u, v = ind_to_adj(N, A, edge, direction)
            A[u][v] -= min_flow
            A[v][u] += min_flow


def get_t_set(s_set):
    """Returns complementary set from s_set in the graph.
    
    Arguments:
        s_set -- Array of indexes of vertices reachable from s.
    """

    t_set = np.ones(s_set.shape[0], dtype="int")
    for v, present in enumerate(s_set):
        if present:
            t_set[v] = 0

    return t_set


def find_minimum_cut(N, s_set):
    """Find minimum st-cut in the graph.
    
    Arguments:
        N -- Incidence matrix.
        s_set -- Array of indexes of vertices reachable from s.
    """

    t_set = get_t_set(s_set)
    cut = np.zeros(N.shape[1], dtype="int")
    for v, present in enumerate(s_set):
        if present:
            cut = get_fringe_edges(N, v, t_set, cut)

    return cut


def get_fringe_edges(N, v, t_set, cut):
    """Get edges that departs from v and end in t_set and add then to the cut.
    
    Returns updated cut.

    Arguments:
        N -- Incidence matrix.
        v -- Vertex inside s_set.
        t_set -- Complement of s_set.
        cut -- Array of edge indexes that has 1 if edge is part of the cut 
               and 0 otherwise.
    """

    for edge, incidence in enumerate(N[v]):
        if incidence:
            for vertex in range(N.shape[0]):
                if N[vertex][edge] == 1 and t_set[vertex]:
                    cut[edge] = 1

    return cut


def solve(c, N):
    """Solve Ford-Fulkerson for a given vector c and matrix N.
    
    Arguments:
        c -- Arch capacity vector.
        N -- Incidence matrix.
    """

    # Weighted adjacency matrix
    A = create_adjacency_matrix(c, N)
    # Number of archs(edges)
    n_archs = N.shape[1]
    # Source vertex
    s = 0
    # Sink vertex
    t = N.shape[0] - 1
    # Max flow in network
    max_flow = 0
    # Flow vector. At the beginning flow is not present.
    flow = np.zeros(n_archs, dtype="int")

    path, s_set = find_path(s, t, A, N)

    while path is not None:
        #print("path", path)
        path = convert_path(path, flow, A, N)
        # Find minimum flow for path
        min_flow = find_min_flow(path, A, N)
        # Update total flow
        max_flow += min_flow
        # Update flows and capacities
        update_graph(N, A, flow, path, min_flow)

        print(path)
        print(flow)
        print(c, "\n", sep="")

        # Find next path
        path, s_set = find_path(s, t, A, N)

    print("Maximum flow:", max_flow)
    min_cut = find_minimum_cut(N, s_set)
    print("Minimum cut:", min_cut)
