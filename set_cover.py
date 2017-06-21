#!/usr/bin/env python3

"Implementation of set cover approximation using the primal-dual method"

__author__ = "João Francisco B. S. Martins"
__email__ = "joaofbsm@dcc.ufmg.br"
__license__ = "GPL"
__version__ = "3.0"

import sys
import numpy as np

def find_fit(vertex, c, N, x):
    """Solve a restriction with equality and return set number and var value.
    
    Arguments:
        vertex -- Variable(vertex) index.
        c -- Set cost vector.
        N -- Correlation matrix.
        x -- Dual solution.
    """

    # Correct matrix for dual
    N = np.transpose(N)
    # Tranpose c so hstack works
    c = np.transpose(c[np.newaxis])
    # Add c to the end of N
    N = np.hstack((N, c))
    # Min value slack is the maximum value possible for the variable
    min_value = float("inf")
    set_num = 0
    i = 0
    for restriction in N:
        if restriction[vertex] != 0:
            # If vertex takes part in the restriction
            slack = restriction[-1] - np.dot(restriction[:-1], x)
            if slack < min_value:
                min_value = slack
                set_num = i
        i += 1
        
    return set_num, min_value


def get_covered(S, N, n_verts):
    """Update the covered vertices by looking at the selected sets.

    Arguments:
        S -- Collection of sets.
        N -- Correlation matrix.
        n_verts -- Number of vertices in the problem.
    """
    N = np.transpose(N)
    covered = np.zeros(n_verts)
    for set_num in S:
        for i in range(n_verts):
            if N[set_num][i] == 1:
                covered[i] = 1

    return covered

def calculate_cost(S, c):
    """Return the cost for the set cover.
    
    Arguments:
        S -- Collection of sets.
        c -- Set cost vector.
    """

    cost = 0
    for set_num in S:
        cost += c[set_num]

    return cost


def solve(c, N):
    """Solve set cover problem for a given cost vector c and matrix N
    
    Arguments:
        c -- Set cost vector.
        N -- Correlation matrix.
    """

    # Initialize vectors
    n_verts = N.shape[0]  # Number of vertices
    n_sets = N.shape[1]  # Number of sets

    S = []  # Sets that cover the elements
    covered = np.zeros(n_verts)  # Covered vertices
    x = np.zeros(n_verts)  # Dual solution
    y = np.zeros(n_sets)  # Primal solution

    # Set precision on float arrays to 3
    np.set_printoptions(precision=3)
    
    while np.count_nonzero(covered) != n_verts:
        # Get uncovered element
        uncovered_vertex = np.where(covered==0)[0][0]
        # Find the max value it can have according to the restrictions
        set_num, var_value = find_fit(uncovered_vertex, c, N, x)
        # Update collection of sets
        S.append(set_num)
        # Update primal solution
        y[set_num] = 1
        # Update dual solution
        x[uncovered_vertex] = var_value
        # Finds which vertex have been covered by the current collection
        covered = get_covered(S, N, n_verts)
        # Print solutions after every iteration
        print(y, "\n", x, "\n\n", sep="")
    cover_cost = calculate_cost(S, c)
    print("Set cover cost:", cover_cost)

