from sympy import Symbol, pi, S, sin, cos, Function, solve, latex

def get_pi_shift(arg):
    """Parse arg to determine pi shifts.  Returns x, n, N where:
    arg = x + n*pi/N
    """
    if arg.is_Add:
        x, m = arg.as_independent(S.Pi)
        if m:
            pi_coef = m/pi
            if pi_coef.is_Rational:
                return x, pi_coef.p, pi_coef.q
            else:
                return x, pi_coef, 1
        else:
            return x, 0

x = Symbol('x')
n = Symbol('n')
N = Symbol('N')
e = x + 2*pi/12
print get_pi_shift(e)

#stop







x = Symbol('x')
print 'sin(pi/2 - x) =',sin(pi/2 - x)
