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
                # if x == 0, it means we can immediately simplify
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

rule_dict = { #'TR1': ((sec(a), S(1)/cos(a)), (csc(a), S(1)/sin(a))),
                  'TR2': ((tan(a), sin(a)/cos(a)), (cot(a), cos(a)/sin(a))),
                  'TR5': ((sin(a)**2, 1 - cos(a)**2)),
                  'TR6': ((cos(a)**2, 1 - sin(a)**2)),
                  'TR7': ((cos(a)**2, (1 + cos(2*a))/2)),
                  'TR8': ((sin(a)*cos(b), (sin(a + b) + sin(a - b))/2),
                          (cos(a)*sin(b), (sin(a + b) - sin(a - b))/2),
                          (cos(a)*cos(b), (cos(a + b) + cos(a - b))/2),
                          (sin(a)*sin(b), -(cos(a + b) - cos(a - b))/2)),
                  'TR9': ((sin(a) + sin(b), 2*sin((a + b)/2)*cos((a - b)/2)),
                          (sin(a) - sin(b), 2*cos((a + b)/2)*sin((a - b)/2)),
                          (cos(a) + cos(b), 2*cos((a + b)/2)*cos((a - b)/2)),
                          (cos(a) - cos(b), -2*sin((a + b)/2)*sin((a - b)/2))),
                  'TR10': ((sin(a + b), sin(a)*cos(b) + cos(a)*sin(b)),
                           (sin(a - b), sin(a)*cos(b) - cos(a)*sin(b)),
                           (cos(a + b), cos(a)*cos(b) - sin(a)*sin(b)),
                           (cos(a - b), cos(a)*cos(b) + sin(a)*sin(b))),
                  'TR11': ((sin(2*a), 2*sin(a)*cos(a)),
                           (cos(2*a), cos(a)**2 - sin(a)**2)),
                  'TR12': ((tan(a + b), (tan(a) + tan(b))/(1 - tan(a)*tan(b))),
                           (tan(a - b), (tan(a) - tan(b))/(1 +
                               tan(a)*tan(b)))),
                  'TR13': ((tan(a)*tan(b), 1 - (tan(a) + tan(b))*cot(a + b)),
                           (cot(a)*cot(b), 1 + (cot(a) + cot(b))*cot(a + b)))}


var("x n N y")
ex1 = 1 - 1/4*sin(2*x)**2 - sin(y)**2 - cos(x)**4
ex2 = cos(pi/9)*cos(2*pi/9)*cos(3*pi/9)*cos(4*pi/9)
ex3 = tan(7*pi/18)+tan(5*pi/18)-sqrt(3)*tan(5*pi/18)*tan(7*pi/18)

def test_get_pi_shift():
    assert get_pi_shift(x+2*pi/12) == (x, 2)
    assert get_pi_shift(y+n*pi) == (y, 12*n)
    assert get_pi_shift(pi/2) == (0, 6)
    assert get_pi_shift(y) == (y, 0)

def test_sin():
    assert sin(-y) == -Sin(y)
    assert sin(pi - y) == Sin(y)
    assert sin(pi + y) == -Sin(y)
    assert sin(2*pi - y) == -Sin(y)
    assert sin(pi/2 + y) == cos(y)
    assert sin(pi/2 - y) == cos(y)
    assert sin(0) == 0
    assert sin(pi/6) == S(1)/2
    assert sin(pi/4) == 1/sqrt(2)
    assert sin(pi/3) == sqrt(3)/2
    assert sin(pi/2) == 1

def test_cos():
    assert cos(-y) == Cos(y)
    assert cos(pi - y) == -Cos(y)
    assert cos(pi + y) == -Cos(y)
    assert cos(2*pi - y) == Cos(y)
    assert cos(pi/2 + y) == -sin(y)
    assert cos(pi/2 - y) == sin(y)
    assert cos(0) == 1
    assert cos(pi/6) == sqrt(3)/2
    assert cos(pi/4) == 1/sqrt(2)
    assert cos(pi/3) == 1/S(2)
    assert cos(pi/2) == 0

def test_tan():
    assert tan(-y) == -Tan(y)
    assert tan(pi - y) == -Tan(y)
    assert tan(pi + y) == Tan(y)
    assert tan(2*pi - y) == -Tan(y)
    assert tan(pi/2 + y) == -cot(y)
    assert tan(pi/2 - y) == cot(y)
    assert tan(0) == 0
    assert tan(pi/6) == 1/sqrt(3)
    assert tan(pi/4) == 1
    assert tan(pi/3) == sqrt(3)
    #assert tan(pi/2) == oo

def test_cot():
    assert cot(-y) == -Cot(y)
    assert cot(pi - y) == -Cot(y)
    assert cot(pi + y) == Cot(y)
    assert cot(2*pi - y) == -Cot(y)
    assert cot(pi/2 + y) == -tan(y)
    assert cot(pi/2 - y) == tan(y)
    #assert cot(0) == 0
    assert cot(pi/6) == sqrt(3)
    assert cot(pi/4) == 1
    assert cot(pi/3) == 1/sqrt(3)
    assert cot(pi/2) == 0
