from sympy import (Symbol, pi, S, Basic, Function, solve, latex, sqrt, sympify,
        var, Wild, symbols)

def sec(a):
    pass

def csc(a):
    pass

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
        Returns cls(x), but takes into account the sign of "x".

        e.g. sin(-x) -> -sin(x), but cos(-x) -> cos(x)
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
        #x, n = get_pi_shift(arg)
        x, r = get_pi_shift2(arg)
        if r == S(0):
            return None
        else:
            return cls.eval_direct2(x, r)
        """
        if n.is_integer:
            m = n % (12 * cls.period)
            if x == 0:
                # if x == 0, it means we can immediately simplify
                return cls.eval_direct(m)
            # Full-period symmetry (2*pi for sin/cos and pi for tan/cot)
            if m % (12 * cls.period) == 0:
                return cls.handle_minus(x)
            else:
                # pi-period symmetry (e.g. only sin/cos, since tan/cot were
                # already handled above for this case)
                if m % 12 == 0:
                    return -cls.handle_minus(x)
                # pi/2-period symmetry (m=6, 18 for sin/cos and m=6 for tan/cot)
                elif m % 6 == 0:
                    f = conjugates[cls]
                    if m == 6:
                        sign = 1
                    else:
                        sign = -1
                    if f.odd:
                        sign = -sign
                    return sign * f.handle_minus(x)
        """

class Sin(TrigFunction):
    odd = True
    period = 2

    @classmethod
    def eval_direct(cls, m):
        """
        Returns the value of sin(2*pi*m/24) where m is an integer.
        """
        return sin_table[m % 24]

    @classmethod
    def eval_direct2(cls, x, m):
        """
        Returns the value of sin(2*pi*m/24) where m is an integer.
        """
        if x == S(0):
            return sin_table[m/S(12) % 24]
        else:


    def as_Cos(self):
        return Cos(pi/2 - self.args[0])


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

    def as_Sin(self):
            return Sin(pi/2 - self.args[0])

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

def get_pi_shift(arg): """
    If arg = x + n*pi/12, returns (x, n), otherwise None.
    """
    x = Wild("x", exclude=[pi])
    n = Wild("n", exclude=[pi])
    r = arg.match(x+n*pi/12)
    # I think it should always match:
    if r is None:
        return arg, S(0)
    else:
        return r[x], r[n]

def get_pi_shift2(arg): """
    If arg = x + c*pi, returns (x, c), otherwise None.
    """
    x = Wild("x", exclude=[pi])
    c = Wild("c", exclude=[pi])
    r = arg.match(x+c*pi)
    # I think it should always match:
    if r is None:
        return arg, S(0)
    else:
        return r[x], r[c] % S(2)*pi

sin = Sin
cos = Cos
tan = Tan
cot = Cot

a,b = map(Wild, 'ab')

rule_dict = { 'TR1': ((sec(a), S(1)/cos(a)), (csc(a), S(1)/sin(a))),
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
                           (tan(a - b), (tan(a) - tan(b))/(1 + tan(a)*tan(b)))),
                  'TR13': ((tan(a)*tan(b), 1 - (tan(a) + tan(b))*cot(a + b)),
                           (cot(a)*cot(b), 1 + (cot(a) + cot(b))*cot(a + b)))}


var("x n N y")
"""
### Three examples from Fu et. al
print '*'*20, '  Example 1  ', '*'*20
print 'Original expression: 1 - 1/S(4)*sin(2*x)**2 - sin(y)**2 - cos(x)**4'
ex1 = 1 - 1/S(4)*sin(2*x)**2 - sin(y)**2 - cos(x)**4
print ex1
ex1_1 = ex1.subs(sin(2*x), 2*sin(x)*cos(x))
print ex1_1
ex1_2 = ex1_1.subs({sin(x)**2: 1-cos(x)**2, sin(y)**2: 1-cos(y)**2}).expand()
print ex1_2
ex1_3 = ex1_2.subs({cos(x)**2: 1 + cos(2*x)/2, cos(y)**2: 1 +
    cos(2*y)/2}).expand()
print ex1_3
ex1_4 = (ex1_3*2).expand().subs(cos(2*y) - cos(2*x), -2*sin(x+y)*sin(y - x))/2
print ex1_4
ex1_5 = ex1_4.subs(-sin(y - x), sin(y - x))
print ex1_5

### Example 2
print '*'*20, '  Example 2  ', '*'*20
print 'Original expression: cos(pi/9)*cos(2*pi/9)*cos(3*pi/9)*cos(4*pi/9)'
ex2 = cos(pi/9)*cos(2*pi/9)*cos(3*pi/9)*cos(4*pi/9)
print ex2
ex2_1 = ex2.subs(cos(4*pi/9), cos(4*pi/9).as_Sin())
print ex2_1

"""
print 'sin(17*pi/18)=',sin(17*pi/18)
print cos(-4*pi/9)
print sin(-4*pi/9)

#ex3 = tan(7*pi/18)+tan(5*pi/18)-sqrt(3)*tan(5*pi/18)*tan(7*pi/18)

def test_get_pi_shift():
    assert get_pi_shift(x+2*pi/12) == (x, 2)
    assert get_pi_shift(y+n*pi) == (y, 12*n)
    assert get_pi_shift(pi/2) == (0, 6)
    assert get_pi_shift(y) == (y, 0)
    assert get_pi_shift(x + y) == (x+y, 0)

def test_sin():
    assert sin(-y) == -sin(y)
    assert sin(pi - y) == sin(y)
    assert sin(pi + y) == -sin(y)
    assert sin(2*pi - y) == -sin(y)
    assert sin(pi/2 + y) == cos(y)
    assert sin(pi/2 - y) == cos(y)
    assert sin(0) == 0
    assert sin(pi/6) == S(1)/2
    assert sin(pi/4) == 1/sqrt(2)
    assert sin(pi/3) == sqrt(3)/2
    assert sin(pi/2) == 1

    assert sin(-y + 2*pi) == -sin(y)
    assert sin(pi - y + 2*pi) == sin(y)
    assert sin(pi + y + 2*pi) == -sin(y)
    assert sin(2*pi - y + 2*pi) == -sin(y)
    assert sin(pi/2 + y + 2*pi) == cos(y)
    assert sin(pi/2 - y + 2*pi) == cos(y)
    assert sin(0 + 2*pi) == 0
    assert sin(pi/6 + 2*pi) == S(1)/2
    assert sin(pi/4 + 2*pi) == 1/sqrt(2)
    assert sin(pi/3 + 2*pi) == sqrt(3)/2
    assert sin(pi/2 + 2*pi) == 1

    assert sin(-y + 4*pi) == -sin(y)
    assert sin(pi - y + 4*pi) == sin(y)
    assert sin(pi + y + 4*pi) == -sin(y)
    assert sin(2*pi - y + 4*pi) == -sin(y)
    assert sin(pi/2 + y + 4*pi) == cos(y)
    assert sin(pi/2 - y + 4*pi) == cos(y)
    assert sin(0 + 4*pi) == 0
    assert sin(pi/6 + 4*pi) == S(1)/2
    assert sin(pi/4 + 4*pi) == 1/sqrt(2)
    assert sin(pi/3 + 4*pi) == sqrt(3)/2
    assert sin(pi/2 + 4*pi) == 1

    assert sin(3*pi/2) == -1
    assert sin(5*pi/2) == 1

def test_cos():
    assert cos(-y) == cos(y)
    assert cos(pi - y) == -cos(y)
    assert cos(pi + y) == -cos(y)
    assert cos(2*pi - y) == cos(y)
    assert cos(pi/2 + y) == -sin(y)
    assert cos(pi/2 - y) == sin(y)
    assert cos(0) == 1
    assert cos(pi/6) == sqrt(3)/2
    assert cos(pi/4) == 1/sqrt(2)
    assert cos(pi/3) == 1/S(2)
    assert cos(pi/2) == 0

    assert cos(-y + 2*pi) == cos(y)
    assert cos(pi - y + 2*pi) == -cos(y)
    assert cos(pi + y + 2*pi) == -cos(y)
    assert cos(2*pi - y + 2*pi) == cos(y)
    assert cos(pi/2 + y + 2*pi) == -sin(y)
    assert cos(pi/2 - y + 2*pi) == sin(y)
    assert cos(0 + 2*pi) == 1
    assert cos(pi/6 + 2*pi) == sqrt(3)/2
    assert cos(pi/4 + 2*pi) == 1/sqrt(2)
    assert cos(pi/3 + 2*pi) == 1/S(2)
    assert cos(pi/2 + 2*pi) == 0

    assert cos(pi) == -1
    assert cos(8*pi) == 1
    assert cos(-9*pi) == -1
    assert cos(3*pi/2) == 0
    assert cos(11*pi/2) == 0
    assert cos(pi/12) == (1 + sqrt(3)) / (2 * sqrt(2))


def test_tan():
    assert tan(-y) == -tan(y)
    assert tan(pi - y) == -tan(y)
    assert tan(pi + y) == tan(y)
    assert tan(2*pi - y) == -tan(y)
    assert tan(pi/2 + y) == -cot(y)
    assert tan(pi/2 - y) == cot(y)
    assert tan(0) == 0
    assert tan(pi/6) == 1/sqrt(3)
    assert tan(pi/4) == 1
    assert tan(pi/3) == sqrt(3)
    #assert tan(pi/2) == oo

    assert tan(-y + pi) == -tan(y)
    assert tan(pi - y + pi) == -tan(y)
    assert tan(pi + y + pi) == tan(y)
    assert tan(2*pi - y + pi) == -tan(y)
    assert tan(pi/2 + y + pi) == -cot(y)
    assert tan(pi/2 - y + pi) == cot(y)
    assert tan(0 + pi) == 0
    assert tan(pi/6 + pi) == 1/sqrt(3)
    assert tan(pi/4 + pi) == 1
    assert tan(pi/3 + pi) == sqrt(3)

    assert tan(7*pi/12) == sin(7*pi/12)/cos(7*pi/12)

def test_cot():
    assert cot(-y) == -cot(y)
    assert cot(pi - y) == -cot(y)
    assert cot(pi + y) == cot(y)
    assert cot(2*pi - y) == -cot(y)
    assert cot(pi/2 + y) == -tan(y)
    assert cot(pi/2 - y) == tan(y)
    #assert cot(0) == 0
    assert cot(pi/6) == sqrt(3)
    assert cot(pi/4) == 1
    assert cot(pi/3) == 1/sqrt(3)
    assert cot(pi/2) == 0

    assert cot(-y + pi) == -cot(y)
    assert cot(pi - y + pi) == -cot(y)
    assert cot(pi + y + pi) == cot(y)
    assert cot(2*pi - y + pi) == -cot(y)
    assert cot(pi/2 + y + pi) == -tan(y)
    assert cot(pi/2 - y + pi) == tan(y)
    assert cot(pi/6 + pi) == sqrt(3)
    assert cot(pi/4 + pi) == 1
    assert cot(pi/3 + pi) == 1/sqrt(3)
    assert cot(pi/2 + pi) == 0

def test_symmetry():
    assert sin(-x) == -sin(x)
    assert cos(-x) == cos(x)
    assert tan(-x) == -tan(x)
    assert cot(-x) == -cot(x)
    assert sin(x+pi) == -sin(x)
    assert sin(x+2*pi) == sin(x)
    assert sin(x+3*pi) == -sin(x)
    assert sin(x+4*pi) == sin(x)
    assert sin(x-5*pi) == -sin(x)
    assert cos(x+pi) == -cos(x)
    assert cos(x+2*pi) == cos(x)
    assert cos(x+3*pi) == -cos(x)
    assert cos(x+4*pi) == cos(x)
    assert cos(x-5*pi) == -cos(x)
    assert tan(x+pi) == tan(x)
    assert tan(x-3*pi) == tan(x)
    assert cot(x+pi) == cot(x)
    assert cot(x-3*pi) == cot(x)
    assert sin(pi/2-x) == cos(x)
    assert sin(3*pi/2-x) == -cos(x)
    assert sin(5*pi/2-x) == cos(x)
    assert cos(pi/2-x) == sin(x)
    assert cos(3*pi/2-x) == -sin(x)
    assert cos(5*pi/2-x) == sin(x)
    assert tan(pi/2-x) == cot(x)
    assert tan(3*pi/2-x) == cot(x)
    assert tan(5*pi/2-x) == cot(x)
    assert cot(pi/2-x) == tan(x)
    assert cot(3*pi/2-x) == tan(x)
    assert cot(5*pi/2-x) == tan(x)
    assert sin(pi/2+x) == cos(x)
    assert cos(pi/2+x) == -sin(x)
    assert tan(pi/2+x) == -cot(x)
    assert cot(pi/2+x) == -tan(x)
