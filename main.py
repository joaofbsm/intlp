#!/usr/bin/env python3

"Execute set cover and Ford-Fulkerson integral LP implementations"

__author__ = "João Francisco B. S. Martins"
__email__ = "joaofbsm@dcc.ufmg.br"
__license__ = "GPL"
__version__ = "3.0"

import sys
import numpy as np
import set_cover as sc
import ford_fulkerson as ff

def read_file(file_path, algorithm):
    f = open(file_path, "r")
    results = None
    if algorithm == "set":
        # Number of vertices
        n_verts = int(f.readline())  
        # Number of sets
        n_sets = int(f.readline())  

        # Get c vector with n_sets entries as a numpy array
        string_list = f.readline().split(" ")[:n_sets]
        int_list = [int(x) for x in string_list]
        c = np.asarray(int_list)

        # Correlation matrix N
        N = np.empty(shape=(n_verts, n_sets))
        for i in range(n_verts):
            string_list = f.readline().split(" ")[:n_sets]
            int_list = [int(x) for x in string_list]
            N[i] = np.asarray(int_list)
        results = (c, N)

    elif algorithm == "flux":
        pass
    else:
        print("[ERROR] Algorithm type does not match any.")
    
    f.close()
    return results


def main(argv):
    if not argv:
        print("usage: main.py set|flux input_file")
    else:
        algorithm = argv[0]
        input_file = argv[1]
        results = read_file(input_file, algorithm)
        if algorithm == "set":
            sc.solve(*results)
        elif algorithm == "flux":
            ff.solve(*results)
        else:
            print("Error")


if __name__ == "__main__":
    main(sys.argv[1:])