import math
import numpy as np

def prime_iteration(A, b, e):
    n = len(A)
    
    d = np.copy(b)
    for i in range(n):
        max_in_line_abs = 0
        max_in_line = 0
        for j in range(n):
            if abs(A[i][j]) > max_in_line_abs:
                max_in_line_abs = abs(A[i][j])
                max_in_line = A[i][j]
        
        d[i] = d[i] / max_in_line
    
    C = (np.zeros((n, n)))
    for i in range(n):
        for j in range(n):
            if i != j:
                C[i][j] = -1 * (A[i][j] / A[i][i])

    q = 0
    for i in range(n):
        temp_sum = 0
        for j in range(n):
            temp_sum += abs(C[i][j])
        
        if temp_sum > q:
            q = temp_sum
    
    for j in range(n):
        temp_sum = 0
        for i in range(n):
            temp_sum += abs(C[i][j])
        
        if temp_sum > q:
            q = temp_sum
    
    X = [d]

    end_value = e
    curent_iteration = 1
    while end_value >= e:
        X.append(np.matmul(C, X[-1]) + d)
        
        end_value = 0
        for j in range(n):
            temp_end_value = 1/(1 - q) * abs(X[-1][j] - X[-2][j])

            if temp_end_value > end_value:
                end_value = temp_end_value
        
        print(f"Iteration: {curent_iteration}")
        curent_iteration += 1
        print(f"x = {np.resize(X[-1], (1, n))}")
        print(f"Ax - b = {np.resize(np.matmul(A, X[-1]) - b, (1, n))}")
        print(f"end_value = {end_value}")
        print("")


    return X[-1]

def main():
    A = np.array([
        [5.02, -0.61, 1.04, 2.18],
        [0.42, 3.95, 1.87, 0.43],
        [3.6, 1.66, 7.07, 0.95],
        [0.88, 0.43, 0.46, 4.44],
    ])
    b = np.array([
        [11.172],
        [0.115],
        [0.009],
        [9.349],
    ])

    # A = np.array([
    #     [10.0, 2.0, -1.0, 2.0],
    #     [1.0, 5.0, 1.0, 0.0],
    #     [1.0, -2.0, -5.0, 1.0],
    #     [3.0, 0.0, 0.0, -9.0],
    # ])
    # b = np.array([
    #     [-4.0],
    #     [1.0],
    #     [2.0],
    #     [10.0],
    # ])

    x = prime_iteration(A=A, b=b, e=0.0001)

    print(f"Ax - b = {np.resize(np.matmul(A, x) - b, (1, len(A)))}")

if __name__ == "__main__":
    np.set_printoptions(precision=4, suppress=True)
    main()
