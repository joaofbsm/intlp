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

def main(argv):
    if not argv:
        print("usage: main.py set|flow input_file")
    else:
        algorithm = argv[0]
        input_file = argv[1]
        results = read_file(input_file, algorithm)
        if algorithm == "set":
            sc.solve(*results)
        elif algorithm == "flow":
            ff.solve(*results)
            pass
        else:
            print("Error")


def read_file(file_path, algorithm):
    f = open(file_path, "r")
    results = None

    # Number of rows
    m = int(f.readline())  
    # Number of columns
    n = int(f.readline())  

    # Get c vector with n entries as a numpy array
    string_list = f.readline().split()[:n]
    if algorithm == "set":
        number_list = [float(x) for x in string_list]
    elif algorithm == "flow":
        number_list = [int(x) for x in string_list]
    c = np.asarray(number_list)

    # Correlation/Incidence matrix N
    N = np.empty(shape=(m, n))
    for i in range(m):
        string_list = f.readline().split()[:n]
        int_list = [int(x) for x in string_list]
        N[i] = np.asarray(int_list)
    results = (c, N)   
    f.close()
    return results


if __name__ == "__main__":
    main(sys.argv[1:])
    