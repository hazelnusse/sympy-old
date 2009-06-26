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

    def __new__(cls, arg, eval=True):
        arg = sympify(arg)
        if not eval:
            return Basic.__new__(cls, arg)
        r = cls.eval(arg)
        if r is not None:
            return r
        else:
            return Basic.__new__(cls, arg)

    @classmethod
    def handle_minus(cls, x):
        """
        Returns cls(x), but takes into account "-".
        """
        if x.could_extract_minus_sign():
            if cls.odd:
                return -cls(-x, eval=False)
            else:
                return cls(-x, eval=False)
        else:
            return cls(x, eval=False)

    @classmethod
    def eval(cls, arg):
        x, n = get_pi_shift(arg)
        if n.is_integer:
            m = n % (12*cls.period)
            if x == 0:
                # if x == 0, it means we can immediatelly simplify
                return cls.eval_direct(m)
            # Full-period symmetry
            if not m % (12*cls.period):
                return cls.handle_minus(x)
            else:
                # Half-period symmetry
                if not m % 12:
                    return -cls.handle_minus(x)
                # Quarter-period symmetry
                elif not m % 6:
                    f = conjugates[cls]
                    sign = (-1)**((((m-6)//12) % cls.period) + f.odd)
                    return sign * f(x)

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

# pi/2-x symmetry:
conjugates = {
    Sin: Cos,
    Cos: Sin,
    Tan: Cot,
    Cot: Tan,
    }

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

sin = Sin
cos = Cos
tan = Tan
cot = Cot


var("x n N y")

def test_get_pi_shift():
    assert get_pi_shift(x+2*pi/12) == (x, 2)
    assert get_pi_shift(y+n*pi) == (y, 12*n)
    assert get_pi_shift(pi/2) == (0, 6)
    assert get_pi_shift(y) == (y, 0)

def test_sin():
    assert sin(0) == 0
    assert sin(y+0) == sin(y)
    assert sin(y+2*pi) == sin(y)
    assert sin(y-2*pi) == sin(y)
    assert sin(pi/2 - x) == cos(x)

    assert sin(pi/2) == 1
    assert sin(pi/3) == sqrt(3)/2

def test_cos():
    assert cos(0) == 1
    assert cos(y+0) == cos(y)
    assert cos(y+2*pi) == cos(y)
    assert cos(y-2*pi) == cos(y)
    assert cos(pi/2 - x) == sin(x)

    assert cos(pi/2) == 0
    assert cos(pi/3) == S(1)/2
