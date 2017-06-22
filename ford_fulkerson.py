#!/usr/bin/env python3

"Implementation of the greedy Ford-Fukerson algorithm for network flows"

__author__ = "João Francisco B. S. Martins"
__email__ = "joaofbsm@dcc.ufmg.br"
__license__ = "GPL"
__version__ = "3.0"

import sys
import numpy as np

def create_weighted_matrix(c, N):
    """Creates matrix with weights based in arch capacity and arch connections.
    
    Arguments:
        c -- Arch capacity vector.
        N -- Incidence matrix.
    """

    vertex_num = N.shape[0]
    # Weighted matrix
    G = np.zeros((vertex_num, vertex_num))

    for i, vertex in enumerate(N):
        for j, edge in enumerate(vertex):
            if edge == -1:
                # If edge points outward the vertex, add its weight
                for k in range(vertex_num):
                    if N[k][j] == 1:
                        # k is vertex to where the edge points
                        G[i][k] = c[j]
                        break

    return G


def find_path(s, t, G):
    """Finds path using depth first search.

    Arguments:
        s -- Source vertex index.
        t -- Target(Sink) vertex index.
        G -- Weighted matrix.
    """

    visited = np.zeros(G.shape[0], dtype="int")
    stack = []
    path = np.zeros(G.shape[0], dtype="int")

    # At the beginning, visit source
    stack.append(s)
    visited[s] = 1

    while stack:
        vertex = stack.pop(0)
        for i, capacity in enumerate(G[vertex]):
            if not visited[i] and capacity > 0:
                stack.append(i)
                visited[i] = 1
                path[i] = vertex
    return path if visited[t] else 0


def solve(c, N):
    """Solve Ford-Fulkerson for a given vector c and matrix N.
    
    Arguments:
        c -- Arch capacity vector.
        N -- Incidence matrix.
    """

    # G combines the forward and backward residual capacities
    G = create_weighted_matrix(c, N)
    orig_G = np.copy(G)
    flow = 0  # Total flow
    s = 0  # Source
    t = N.shape[0] - 1  # Target

    path = find_path(s, t, G)
    while np.count_nonzero(path):
        path_flow = float("inf")
        v = t
        while v != s:
            if G[path[v]][v] < path_flow:
                path_flow = G[path[v]][v]
            v = path[v]

        flow += int(path_flow)

        v = t
        while v != s:
            u = path[v]
            G[u][v] -= path_flow
            G[v][u] += path_flow
            v = u

        path = find_path(s, t, G)

    print(flow)

    for i in range(G.shape[0]):
        for j in range(G.shape[1]):
            if G[i][j] == 0 and orig_G[i][j] != 0:
                print(i, "-", j, sep="")
