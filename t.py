from sympy import (Symbol, pi, S, Basic, Function, solve, latex, sqrt, sympify,
        var, Wild)

C0 = (sqrt(3)-1)/(2*sqrt(2))
C1 = S(1)/2
C2 = sqrt(2)/2
C3 = sqrt(3)/2
C4 = (sqrt(3)+1)/(2*sqrt(2))


# sin_table[n] represents the value of sin(2*pi*n/24) for n = 0..23
sin_table = [
        0,  C0,  C1,  C2,  C3,  C4,  1,  C4,  C3,  C2,  C1,  C0,
        0, -C0, -C1, -C2, -C3, -C4, -1, -C4, -C3, -C2, -C1, -C0
        ]

class TrigFunction(Basic):

    def __new__(cls, arg):
        arg = sympify(arg)
        r = cls.eval(arg)
        if r is not None:
            return r
        else:
            return Basic.__new__(cls, arg)

    @classmethod
    def eval(cls, arg):
        x, n = get_pi_shift(arg)
        if n.is_integer:
            m = n % (12*cls.period)
            if x == 0:
                # if x == 0, it means we can immediatelly simplify
                return cls.eval_direct(m)

class Sin(TrigFunction):
    odd = True
    period = 2

    @classmethod
    def eval_direct(cls, m):
        """
        Returns the value of sin(2*pi*m/24) where m is an integer.
        """
        return sin_table[m % 24]


class Cos(TrigFunction):
    odd = False
    period = 2

    @classmethod
    def eval_direct(cls, m):
        """
        Returns the value of cos(2*pi*m/24) where m is an integer.
        """
        # we use the relation cos(2*pi*m/24) = sin(2*pi*(m+6)/24)
        return Sin.eval_direct(m+6)

class Tan(TrigFunction):
    odd = True
    period = 1

    @classmethod
    def eval_direct(cls, m):
        """
        Returns the value of tan(2*pi*m/24) where m is an integer.
        """
        # we use the relation tan(x) = sin(x)/cos(x)
        return Sin.eval_direct(m)/Cos.eval_direct(m)

class Cot(TrigFunction):
    odd = True
    period = 1

    @classmethod
    def eval_direct(cls, m):
        """
        Returns the value of cot(2*pi*m/24) where m is an integer.
        """
        # we use the relation cot(x) = cos(x)/sin(x)
        return Cos.eval_direct(m)/Sin.eval_direct(m)

def get_pi_shift(arg):
    """
    If arg = x + n*pi/12, returns (x, n), otherwise None.
    """
    x = Wild("x", exclude=[pi])
    n = Wild("n", exclude=[pi])
    r = arg.match(x+n*pi/12)
    # I think it should always match:
    assert r is not None
    return r[x], r[n]

var("x n N y")
print get_pi_shift(x+2*pi/12)
print get_pi_shift(y+n*pi)
print get_pi_shift(pi/2)
print get_pi_shift(y)

print Sin(0)
print Sin(pi/2)
print Sin(pi)
print Cos(0)
print Cos(pi/2)
print Cos(pi)

print Sin(pi/5)
