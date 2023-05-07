import math
import numpy as np

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

def _get_system_of_equations(A: np.ndarray, y_0: np.ndarray) -> tuple:
    """Retrun system (B, b): B x = b"""
    n = len(y_0)

    B = np.zeros((n, n))

    y_k = np.array(y_0, copy=True)
    for k in range(n):
        for i in range(n):
            B[i, n-1-k] = y_k[i]

        y_k = np.matmul(A, y_k)  # The last multiplication is nesesary
        print(f"y_{k+1}:")
        print(y_k)
    
    # y_k is b vector
    if n % 2 == 0:
        y_k = -1 * y_k

    return B, y_k

def _func(x: float, p: np.ndarray) -> float:
    n = len(p)
    rez = 0

    temp = 1
    for i in range(n-1, -1, -1):
        rez += temp * p[i][0]
        temp *= x
        # print(f">> p_{i}: {p[i][0]}, temp: {temp}, rez: {rez}")
    
    rez += temp

    return rez

def Krilov_method(A: np.ndarray, y_0: np.ndarray) -> np.ndarray:
    n = len(y_0)

    B, b = _get_system_of_equations(A=A, y_0=y_0)
    print(f"B:\n{B}")
    print(f"b:\n{b}")

    p = gaus(A=B, b=b)
    print(f"Polinom coeficients:\n{p}")

    return p

def main():
    A = np.array([
        [6.30, 1.07, 0.99, 1.20],
        [1.07, 4.12, 1.30, 0.16],
        [0.99, 1.30, 5.48, 2.10],
        [1.20, 0.16, 2.10, 6.06],
    ])

    y_0 = np.array([
        [1],
        [1],
        [1],
        [1],
    ])

    p = Krilov_method(A=A, y_0=y_0)


def main_2():
    import lab_1

    p=np.array([
        [ -21.96  ],
        [ 169.721 ],
        [-550.4405],
        [ 631.388 ],
    ])

    point_1 = lab_1.bisection(a=0, b=10, e=0.0001, f=lambda x: _func(x=x, p=p))[0]
    print(f"rez = {point_1:0.4f}")

    point_2 = lab_1.bisection(a=point_1 + 0.0001, b=10, e=0.0001, f=lambda x: _func(x=x, p=p))[0]
    print(f"rez = {point_2:0.4f}")

    point_3 = lab_1.bisection(a=0, b=point_1 - 0.0001, e=0.0001, f=lambda x: _func(x=x, p=p))[0]
    print(f"rez = {point_3:0.4f}")

    point_4 = lab_1.bisection(a=point_3 + 0.0001, b=point_1 - 0.0001, e=0.0001, f=lambda x: _func(x=x, p=p))[0]
    print(f"rez = {point_4:0.4f}")
    

if __name__ == "__main__":
    np.set_printoptions(precision=4, suppress=True)
    main()
    print("///////////////////////////////")
    main_2()

