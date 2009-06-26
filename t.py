from sympy import Symbol, pi, S, Basic, Function, solve, latex, sqrt

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
        r = cls.eval(arg)
        if r:
            return r
        else:
            return Basic.__new__(cls, arg)

    @classmethod
    def eval(cls, arg):
        pass

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

print Sin(0)
print Sin(pi/2)
print Sin(pi)
print Cos(0)
print Cos(pi/2)
print Cos(pi)

