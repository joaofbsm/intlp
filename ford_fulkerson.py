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

def find_dest(e, s, N):
    """Find destination for edge e.
    
    Arguments:
        e -- [description]
        N -- [description]
    """

    for v in range(N.shape[0]):
        if N[v][e] != 0 and v != s:
            return v

    return None


def find_path(s, t, cf, N, path=None):
    """Finds path using recursive depth first search.

    Arguments:
        s -- Source vertex index.
        t -- Target(Sink) vertex index.
    """
    #print("Source vertex:", s)

    if path is None:
        path = np.zeros(N.shape[1], dtype="int")
    if s == t:
        #print("FOUND END")
        return path

    for edge, incidence in enumerate(N[s]):
        if incidence != 0:
            # Edge connects to s
            next_vertex = find_dest(edge, s, N)

            if incidence == -1 and path[edge] == 0 and cf[edge] > 0:
                path[edge] = 1
                results = find_path(next_vertex, t, cf, N, path)
                if results is not None:
                    return results
            elif incidence == 1 and path[edge] == 0 and cf[edge] == 0:
                path[edge] = -1
                results = find_path(next_vertex, t, cf, N, path)
                if results is not None:
                    return results
    # Couldn't get to t
    return None


def find_min_flow(path, c):
    min_flow = float("inf")

    for i, edge in enumerate(path):
        if edge != 0 and 0 < c[i] < min_flow:
            min_flow = c[i]

    return min_flow


def update_flow(flow, path, min_flow):
    for i, edge in enumerate(path):
        if edge != 0:
            flow[i] += min_flow * edge

    return flow


"""def find_cut(c, cf):
    n_archs = c.shape[0]
    cut = np.zeros(n_archs)
    for edge in range(n_archs):
        if cf[edge] == 0 and c[edge] != 0:
            cut[edge] = 1

    return cut
"""


def solve(c, N):
    """Solve Ford-Fulkerson for a given vector c and matrix N.
    
    Arguments:
        c -- Arch capacity vector.
        N -- Incidence matrix.
    """

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
    # Residual capacities
    cf = np.copy(c)

    path = find_path(s, t, cf, N)

    while path is not None:
        # Find minimun flow for path
        min_flow = find_min_flow(path, cf)
        # Update total flow
        max_flow += min_flow
        # Update flows
        flow = update_flow(flow, path, min_flow)
        # Update residual capacities
        cf = np.subtract(c, flow)

        print(path)
        print(flow)
        print(c, "\n", sep="")

        # Find next path
        path = find_path(s, t, cf, N)

    print(max_flow)
