import math
import matplotlib.pyplot as plt
import numpy as np


def plot_table(function):
    x = np.arange(0.0, 4.0, 0.01)
    y = [function(n) for n in np.arange(0.0, 4.0, 0.01)]

    fig, ax = plt.subplots()
    ax.plot(x, y)

    ax.set(xlabel='x', ylabel='y')
    ax.grid()

    plt.show()

def get_function_table(partion: int=10):
    """returns table with values (x, f(x))
    """
    step = 4 / (partion - 1)
    f_table_dict = dict()
    f_table_list = list()

    x = 0
    for i in range(partion - 1):
        f_table_dict[x] = x * math.sqrt(x)
        f_table_list.append((x, x * math.sqrt(x)))
        x += step

    f_table_dict[4] = 4 * math.sqrt(4)
    f_table_list.append((4, 4 * math.sqrt(4)))

    return f_table_dict, f_table_list

def lagrangian_interpolation(x: float, f_table: dict):
    """Calculates function value
    """
    f_coef = list()

    rez_x = 0
    for key_i in f_table:
        temp = f_table[key_i]
        
        # upper part
        for key_j in f_table:
            if key_j != key_i:
                temp = temp * (x - key_j)

        # lower part
        for key_j in f_table:
            if key_j != key_i:
                temp = temp / (key_i - key_j)
        
        rez_x += temp
    
    return rez_x

def _gaus_forward(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    A_cur = np.array(A, copy=True)
    b_cur = np.array(b, copy=True)
    
    n = len(A_cur)

    for j in range(n):
        # if diagonal element is zero
        if A_cur[j][j] == 0:
            big = 0
            k_row = j
        
            for k in range(j + 1, n):
                if abs(A_cur[k][j]) > big:
                    big = abs(A_cur[k][j])
                    k_row = k
            
            for l in range(j, n):
                A_cur[j][l], A_cur[k_row][l] = A_cur[k_row][l], A_cur[j][l]

            b_cur[j], b_cur[k_row] = b_cur[k_row], b_cur[j]

        pivot = A_cur[j][j]

        # error case
        if pivot == 0:
            raise ValueError("Given matrix is singular")

        # main part
        for i in range(j + 1, n):
            mult = A_cur[i][j] / pivot

            for l in range(j, n):
                A_cur[i][l] = A_cur[i][l] - mult * A_cur[j][l]
            
            b_cur[i] = b_cur[i] - mult * b_cur[j]
        
    return A_cur, b_cur

def _gaus_backward(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    n = len(A)
    X = np.zeros((n, 1))

    for i in range(n-1, -1, -1):
        sum = 0

        for j in range(i+1, n):
            sum = sum + X[j] * A[i][j]
        
        X[i] = 1 / A[i][i] * (b[i] - sum)
    
    return X

def gaus(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    A, b = _gaus_forward(A=A, b=b)
    X = _gaus_backward(A=A, b=b)

    return X

def get_cubic_spline_coefs(f_table: list) -> list:
    n = len(f_table) - 1

    # coefs value: [..., [a_i, b_i, c_i, d_i], ...]
    coefs = [[None for j in range(4)] for i in range(n)]
    
    # Calculate a_i
    for i in range(n):
        coefs[i][0] = f_table[i][1]


    # calculate c_i from Ac = b
    A = np.zeros((n, n))
    A[0][0] = 1
    A[n-1][n-1] = 1
    b = np.zeros((n, 1))

    for i in range(1, n - 1):
        A[i][i - 1] = (f_table[i][0] - f_table[i - 1][0])
        A[i][i] = 2 * (f_table[i][0] - f_table[i - 1][0] + f_table[i+1][0] - f_table[i][0])
        A[i][i + 1] = (f_table[i + 1][0] - f_table[i][0])

        b[i][0] = 3 * ((f_table[i + 1][1] - f_table[i][1]) / (f_table[i + 1][0] - f_table[i][0])\
                     - (f_table[i][1] - f_table[i - 1][1]) / (f_table[i][0] - f_table[i - 1][0]))

    c = gaus(A=A, b=b)

    for i in range(n):
        coefs[i][2] = c[i][0]

    # Calculate d_i
    for i in range(n - 1):
        coefs[i][3] = (coefs[i + 1][2] - coefs[i][2]) / (3 * (f_table[i+1][0] - f_table[i][0]))

    coefs[n-1][3] = - coefs[n-1][2] / (3 * (f_table[n-1][0] - f_table[n-2][0]))

    # Calculate b_i
    for i in range(n - 1):
        coefs[i][1] = (f_table[i+1][1] - f_table[i][1]) / (f_table[i+1][0] - f_table[i][0])\
                    - (f_table[i+1][0] - f_table[i][0]) / 3 * (coefs[i + 1][2] + 2 * coefs[i][2])

    coefs[n-1][1] = (f_table[n - 1][1] - f_table[n - 2][1]) / (f_table[n - 1][0] - f_table[n - 2][0])\
                  - 2 * (f_table[n - 1][0] - f_table[n - 2][0]) * coefs[n - 1][2] / 3

    return coefs

def cubic_spline(x: float, f_table: list, coefs: list) -> float:
    if x <  f_table[0][0] or x > f_table[-1][0]:
        raise ValueError()
    
    n = len(f_table)

    # find the interval index: x < x_i
    i = 1
    while x > f_table[i][0]:
        i += 1

    # Calculate value
    i -= 1
    rez = coefs[i][0] + coefs[i][1] * (x - f_table[i][0])\
            + coefs[i][2] * pow(x - f_table[i][0], 2) + coefs[i][3] * pow(x - f_table[i][0], 3)
      
    return rez

def main():
    f_table_dict, f_table_list = get_function_table(partion=20)

    # plot function with lagrangian interpolation
    plot_table(lambda x: lagrangian_interpolation(x=x, f_table=f_table_dict))

    # plot function with spline interpolation
    coefs = get_cubic_spline_coefs(f_table=f_table_list)
    plot_table(lambda x: cubic_spline(x=x, f_table=f_table_list, coefs=coefs))

    # plot fault graph
    plot_table(lambda x: lagrangian_interpolation(x=x, f_table=f_table_dict) - x * math.sqrt(x))
    plot_table(lambda x: cubic_spline(x=x, f_table=f_table_list, coefs=coefs) - x * math.sqrt(x))

if __name__ == "__main__":
    main()

