import math


def print_matrix(M, pre_text=""):
    print(f"{pre_text}[")
    for l in M:
        print(f"    [", end="")
        for i in range(len(l) - 1):
            print(f"{l[i]:0.6f}", end=", ")
        
        print(f"{l[-1]:0.6f}],")
    
    print("]")

def print_vector(v, pre_text=""):
    print(f"{pre_text}[", end="")
    for i in range(len(v) - 1):
        print(f"{v[i]:0.6f}", end=", ")
    
    print(f"{v[-1]:0.6f}]")
    

def transpose(M):
    n = len(M)
    rez = [[0 for i in range(n)] for i in range(n)]
    for l in range(len(M)):
        for c in range(len(M)):
            rez[l][c] = M[c][l]

    return rez

def forward_Gaus(A):
    n = len(A)
    T = [[0 for i in range(n)] for i in range(n)]

    # 1)
    T[0][0] = math.sqrt(A[0][0])

    print_matrix(T)

    # 2)
    for l in range(1, n):
        T[0][l] = A[0][l] / T[0][0]
    
    print_matrix(T)

    # 3), 4)
    for j in range(1, n):
        for i in range(1, j):
            temp_sum = 0
            for k in range(0, i):
                temp_sum += T[k][i] * T[k][j]

            T[i][j] = (A[i][j] - temp_sum) / T[i][i]

            print_matrix(T)
        
        temp_sum = 0
        for k in range(0, j):
            temp_sum += T[k][j] * T[k][j]

        T[j][j] = math.sqrt(A[j][j] - temp_sum)

        print_matrix(T)

    return T

def sqrt_Gaus(T_upper, b):
    n = len(T_upper)
    y = [None for i in range(n)]

    y[0] = b[0] / T_upper[0][0]
    for i in range(1, n):
        temp_sum = 0
        for k in range(i):
            temp_sum += T_upper[k][i] * y[k]

        y[i] = (b[i] - temp_sum) / T_upper[i][i]

    x = [None for i in range(n)]
    x[n-1] = y[n-1] / T_upper[n-1][n-1]

    for i in range(n-2, -1, -1):
        temp_sum = 0
        for k in range(i+1, n):
            temp_sum += T_upper[i][k] * x[k]

        x[i] = (y[i] - temp_sum) / T_upper[i][i]

    return x

def check_answer(A, b, x):
    n = len(A)
    rez =[0 for i in range(n)]

    for i in range(n):
        element_sum = 0
        for k in range(n):
            element_sum += A[i][k] * x[k]
        
        rez[i] = element_sum - b[i]

    return rez

def square_root_method(A, b):
    T = forward_Gaus(A)
    print_matrix(M=T, pre_text="T' = ")

    x = sqrt_Gaus(T_upper=T, b=b)
    print_vector(x, pre_text="x = ")

    rez = check_answer(A=A, b=b, x=x)
    print_vector(rez, pre_text="AX-b = ")
    # print(rez)


def main():
    A = [
        [2.12, 0.42, 1.34, 0.88],
        [0.42, 3.95, 1.87, 0.43],
        [1.34, 1.87, 2.98, 0.46],
        [0.88, 0.43, 0.46, 4.44],
    ]
    b = [
        11.172,
        0.115,
        0.009,
        9.349,
    ]
    
    square_root_method(A=A, b=b)

if __name__ == "__main__":
    main()


