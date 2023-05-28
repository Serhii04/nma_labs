import math
import matplotlib.pyplot as plt
import numpy as np


# ************************************
#           gieven functions
# ************************************

def exact_solution(x):
    return x * math.cos(x)

def equation(x, y):
    return (1 - pow(x, 2)) * y + math.cos(x) - x * math.sin(x) - (1 - pow(x, 2)) * (x * math.cos(x))

# ************************************
#               Methods
# ************************************

def Renge_Kut(f, x_0, y_0, h, r):
    x = list()
    y = list()

    x_i = x_0
    y_i = y_0
    x.append(x_i)
    y.append(y_i)

    for i in range(r):
        k_1 = h * f(x_i, y_i)
        k_2 = h * f(x_i + 0.5 * h, y_i + 0.5 * k_1)
        k_3 = h * f(x_i + 0.5 * h, y_i + 0.5 * k_2)
        k_4 = h * f(x_i + h, y_i + k_3)

        x_i += h
        y_i = y_i + 1/6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4)
        x.append(x_i)
        y.append(y_i)

    return x, y

def Adamas_Boshford(f, x_0, y_0, h, r, m=4):
    x, y = Renge_Kut(f, x_0, y_0, h, r=(m - 1))

    x_i = x[-1]
    y_i = y[-1]
    for i in range(r - m):
        y_i = y_i + h * (55 * f(x[-1], y[-1]) - 59 * f(x[-2], y[-2]) + 37 * f(x[-3], y[-3]) - 9 * f(x[-4], y[-4])) / 24
        x_i = x_i + h

        x.append(x_i)
        y.append(y_i)

    return x, y

# ************************************
#              Example
# ************************************

def plot_all_graphs():
    from_ = -1 * np.pi
    to_ = 1 * np.pi
    h = 0.1

    # exact graph points
    x = np.arange(from_, to_, 0.01)
    y = [exact_solution(n) for n in np.arange(from_, to_, 0.01)]

    # aproximate Renge-Kut points
    x_rk, y_rk = Renge_Kut(
        x_0=from_,
        y_0=exact_solution(from_),
        f=equation,
        h=h,
        r=int((to_ - from_) / h)
    )

    # aproximate Adamas-Boshford points
    x_ab, y_ab = Adamas_Boshford(
        x_0 = from_,
        y_0=exact_solution(from_),
        f=equation,
        h=h,
        r=int((to_ - from_) / h)
    )

    # plot graphs
    fig, ax = plt.subplots()
    ax.plot(x, y, 'b-', label='real')
    ax.plot(x_rk, y_rk, 'g-', label='line 1')
    ax.plot(x_ab, y_ab, 'c-', label='line 2')

    ax.set(xlabel='x', ylabel='y')
    ax.grid()

    plt.show()

    # Calculation error
    fig, ax = plt.subplots()

    dy_rk = list()
    for i in range(len(x_rk)):
        dy_rk.append(y_rk[i] - exact_solution(x_rk[i]))

    dy_ab = list()
    for i in range(len(x_ab)):
        dy_ab.append(y_ab[i] - exact_solution(x_ab[i]))


    ax.plot(x_rk, dy_rk, 'g-', label='line 1')
    ax.plot(x_ab, dy_ab, 'c-', label='line 2')

    ax.set(xlabel='x', ylabel='y')
    ax.grid()

    plt.show()

if __name__ == "__main__":
    plot_all_graphs()