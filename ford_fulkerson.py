#!/usr/bin/env python3

"Implementation of the greedy Ford-Fukerson algorithm for network flows"

__author__ = "João Francisco B. S. Martins"
__email__ = "joaofbsm@dcc.ufmg.br"
__license__ = "GPL"
__version__ = "3.0"

import sys
import numpy as np

# TODO
# - Implement return edge flow = -1 * edge flow

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

    Arguments:
        s -- Source vertex index.
        t -- Target(Sink) vertex index.
    """

    visited = np.zeros(N.shape[0], dtype="int")
    stack = []
    path = np.zeros(N.shape[0], dtype="int")

    # Begin by visiting source
    stack.append((s, [s]))

    while stack:
        u, path = stack.pop()
        visited[u] = 1
        for v, capacity in enumerate(A[u]):
            if visited[v] == 0 and capacity > 0:
                if v == t:
                    return path + [v], None
                else:
                    stack.append((v, path + [v]))

    return None, visited


def find_edge(N, u, v):
    """Convert adjacency edge to incidence edge number.
    
    [description]
    
    Arguments:
        A -- [description]
        N -- [description]
        u -- [description]
        v -- [description]
    """

    for edge, incidence in enumerate(N[u]):
        if incidence != 0 and N[v][edge] != 0:
            #print("Find edge -> u", u, "v", v, "edge", edge, "direction", N[v][edge])
            return edge, N[v][edge]


def convert_path(path, N):
    new_path = np.zeros(N.shape[1], dtype="int")
    for i in range(len(path) - 1):
        u = path[i]
        v = path[i + 1]
        edge, direction = find_edge(N, u, v)
        new_path[edge] = direction

    return new_path

def ind_to_adj(N, A, edge, direction):
    u = -1
    v = -1
    for i in range(N.shape[0]):
        if N[i][edge] == direction * -1:
            u = i
        elif N[i][edge] == direction:
            v = i
    return u, v
                

def find_min_flow(path, A, N):
    min_flow = float("inf")
    for edge, direction in enumerate(path):
        u, v = ind_to_adj(N, A, edge, direction)
        if A[u][v] != 0:
            flow = A[u][v]
            min_flow = min(min_flow, flow)
    return int(min_flow)


def update_graph(N, A, flow, path, min_flow):
    for edge, direction in enumerate(path):
        if direction != 0:
            flow[edge] += min_flow * direction
            u, v = ind_to_adj(N, A, edge, direction)
            A[u][v] -= min_flow
            A[v][u] += min_flow






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
        path = convert_path(path, N)
        # Find minimun flow for path
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

    print(max_flow)

