def f(x):
    return (- pow(x, 4) + 3 * pow(x, 3) - 2 * x + 4)

def f_(x):
    return (- 4*pow(x, 3) + 9*pow(x, 2) - 2)

def bisection(f, a, b, e):
    iteration = 0
    while abs(a-b) > e:
        iteration += 1
        c = (a+b)/2
        
        if f(a)*f(c) <= 0:
            b = c
        else:
            a = c

    return (a+b)/2, iteration

def hord_method(f, a, b, e):
    iteration = 0

    c_prew = None
    c = c = (a*f(b) - b*f(a))/((f(b) - f(a)))

    if f(a)*f(c) <= 0:
        b = c
    else:
        a = c

    while abs(f(c)) > e:
        iteration += 1
        c_prew = c
        c = (a*f(b) - b*f(a))/((f(b) - f(a)))

        if abs(c - c_prew) < e:
            return (a+b)/2, iteration
        
        if f(a)*f(c) <= 0:
            b = c
        else:
            a = c

    return (a+b)/2, iteration

def newton_method(f, f_, a, b, e):
    iteration = 0
    x = (a+b)/2

    while abs(f(x)) > e:
        iteration += 1
        
        x = x - f(x) / f_(x)

    return (a+b)/2, iteration

def main():
    print("\nBisection:")
    rez = bisection(f, -5, 0, 0.001)
    print(f"res: {rez[0]}, iteration: {rez[1]}")
    
    rez = bisection(f, 0, 5, 0.001)
    print(f"res: {rez[0]}, iteration: {rez[1]}")
    
    
    print("\nHord method:")
    rez = hord_method(f, -5, 0, 0.001)
    print(f"res: {rez[0]}, iteration: {rez[1]}")
    
    rez = hord_method(f, 0, 5, 0.001)
    
    print(f"res: {rez[0]}, iteration: {rez[1]}")
    
    
    print("\nNewton method:")
    rez = newton_method(f, f_, -5, 0, 0.000001)
    print(f"res: {rez[0]}, iteration: {rez[1]}")
    
    rez = newton_method(f, f_, 0, 5, 0.000001)
    print(f"res: {rez[0]}, iteration: {rez[1]}")


if __name__ == "__main__":
    main()