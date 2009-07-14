from sympy import (Symbol, pi, S, Basic, Function, solve, latex, sqrt, sympify,
        var, Wild, symbols, floor, Rational, I, E, zoo, exp)

C0 = (sqrt(3)-1)/(2*sqrt(2))
C1 = S(1)/2
C2 = 1/sqrt(2)
C3 = sqrt(3)/2
C4 = (sqrt(3)+1)/(2*sqrt(2))

# sin_table[n] represents the value of sin(2*pi*n/24) for n = 0..23
sin_table = [
        0,  C0,  C1,  C2,  C3,  C4,  1,  C4,  C3,  C2,  C1,  C0,
        0, -C0, -C1, -C2, -C3, -C4, -1, -C4, -C3, -C2, -C1, -C0
        ]
# sin_table2[n] represents the value of sin(n*pi/10) for n = 0...19
C02 = (sqrt(5)-1)/4
C12 = sqrt(5/S(8) - sqrt(5)/8)
C22 = (sqrt(5)+1)/4
C32 = sqrt(5/S(8) + sqrt(5)/8)

sin_table2 = [
        0,  C02,  C12,  C22,  C32,   1,  C32,  C22,  C12,  C02,
        0, -C02, -C12, -C22, -C32,  -1, -C32, -C22, -C12, -C02]

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
        # Match arg = a + b*pi
        a, b = get_pi_shift(arg)

        # Case 1:  a == 0 and b == 0
        if a == 0 and b == 0:
            return cls.eval_direct(0)
        # Case 2:  a != 0 and b == 0
        elif a != 0 and b == 0:
            if a.is_imaginary:
                return cls.hyper(a.coeff(I))
            else:
                return cls.handle_minus(a)
        # Case 3: a == 0 and b != 0
        elif a == 0 and b != 0:
            if a.is_imaginary:
                return cls.hyper(b.coeff(I))
            else:
                return cls.eval_direct(b)
        # Case 4:  a != 0 and b != 0
        elif a != 0 and b != 0:
            if a.is_imaginary and b.is_imaginary:
                return cls.hyper(a.coeff(I) + b.coeff(I)*pi)
            elif b.is_rational and b.is_real:
                # Bring it to inside of the period
                b = b % 2
                oct = cls.determine_octant(b)
                b_mod = (b % (1/S(4))) if b != 1/S(4) else 1/S(4)
                return cls.eval_indirect(a, b, b_mod, oct)
            else:
                return cls(a + b*pi, eval=False)

    @classmethod
    def determine_octant(cls, b):
        # Determine octant
        if 0 <= b <= 1/S(4):
            oct = 1
        elif 1/S(4) < b <= 1/S(2):
            oct = 2
        elif 1/S(2) < b <= 3/S(4):
            oct = 3
        elif 3/S(4) < b <= S(1):
            oct = 4
        elif S(1) < b <= 5/S(4):
            oct = 5
        elif 5/S(4) < b <= S(3)/2:
            oct = 6
        elif S(3)/2 < b <= 7/S(4):
            oct = 7
        else:
            oct = 8
        return oct

class Sin(TrigFunction):
    odd = True
    period = 2

    @classmethod
    def eval_direct(cls, b):
        """
        Returns the value of sin(b*pi) when direct evaluation is possible.
        """
        b = sympify(b)
        if b.is_rational:
            if b.is_integer:
                return 0
            elif (b*10).is_integer:
                return sin_table2[int(b*10) % 20]
            elif (b*12).is_integer:
                return sin_table[int(b*12) % 24]
            else:
                b = b % 2
                oct = cls.determine_octant(b)
                if oct == 1:
                    return cls.handle_minus(b*pi)
                elif oct == 2:
                    return Cos.handle_minus((1/S(2) - b)*pi)
                elif oct == 3:
                    return Cos.handle_minus((b % (1/S(2)))*pi)
                elif oct == 4:
                    return cls.handle_minus((1 - b)*pi)
                elif oct == 5:
                    return -cls.handle_minus((b % 1)*pi)
                elif oct == 6:
                    return -Cos.handle_minus((3/S(2) - b)*pi)
                elif oct == 7:
                    return -Cos.handle_minus((b % (3/S(2)))*pi)
                elif oct == 8:
                    return -cls.handle_minus((2 - b)*pi)
        else:
            return cls.handle_minus(b*pi)

    @classmethod
    def eval_indirect(cls, a, b, b_mod, oct):
        """
        Puts any pi-shifts into the interval (0, pi/4)
        """
        if oct == 1:
            return cls.handle_minus(a + b_mod*pi)
        elif oct == 2:
            return Cos.handle_minus(-a + b_mod*pi)
        elif oct == 3:
            return Cos.handle_minus(a + b_mod*pi)
        elif oct == 4:
            return cls.handle_minus(-a + b_mod*pi)
        elif oct == 5:
            return -cls.handle_minus(a + b_mod*pi)
        elif oct == 6:
            return -Cos.handle_minus(-a + b_mod*pi)
        elif oct == 7:
            return -Cos.handle_minus(a + b_mod*pi)
        elif oct == 8:
            return -cls.handle_minus(-a + b_mod*pi)

class Cos(TrigFunction):
    odd = False
    period = 2

    @classmethod
    def eval_direct(cls, b):
        """
        Returns the value of cos(b*pi) when direct evaluation is possible.
        """
        b = sympify(b)
        if b.is_rational:
            if b.is_integer:
                if b.is_even:
                    return 1
                elif b.is_odd:
                    return -1
            elif (b*10).is_integer:
                return sin_table2[int(b*10 + 5) % 20]
            elif (b*12).is_integer:
                return sin_table[int(b*12 + 6) % 24]
            else:
                b = b % 2
                oct = cls.determine_octant(b)
                if oct == 1:
                    return cls.handle_minus(b*pi)
                elif oct == 2:
                    return Cos.handle_minus((1/S(2) - b)*pi)
                elif oct == 3:
                    return Cos.handle_minus((b % (1/S(2)))*pi)
                elif oct == 4:
                    return -cls.handle_minus((1 - b)*pi)
                elif oct == 5:
                    return -cls.handle_minus((b % 1)*pi)
                elif oct == 6:
                    return -Cos.handle_minus((3/S(2) - b)*pi)
                elif oct == 7:
                    return -Cos.handle_minus((b % (3/S(2)))*pi)
                elif oct == 8:
                    return cls.handle_minus((2 - b)*pi)
        else:
            return cls.handle_minus(b*pi)

    @classmethod
    def eval_indirect(cls, a, b, b_mod, oct):
        return Sin.eval(a + (b + 1/S(2))*pi)

    def as_Sin(self):
            return Sin(pi/2 - self.args[0])

class Tan(TrigFunction):
    odd = True
    period = 1

    @classmethod
    def eval_direct(cls, b):
        """
        Returns the value of cos(b*pi) when direct evaluation is possible.
        """
        return Sin.eval_direct(b)/Cos.eval_direct(b)
        """
        b = sympify(b)
        if b.is_rational:
            if (b*12).is_integer:
                return sin_table[int(b*12) % 24]/sin_table[int(b*12 + 6) % 24]
            elif b.q == 5:
                sign = 1 if (b.p % 2 == 0) else -1
                return Sin.eval_direct(b)/Cos.eval_direct(b)
            else:
                b = b % 2
                oct = cls.determine_octant(b)
                if oct == 1:
                    return cls.handle_minus(b*pi)
                elif oct == 2:
                    return Cos.handle_minus((1/S(2) - b)*pi)
                elif oct == 3:
                    return Cos.handle_minus((b % (1/S(2)))*pi)
                elif oct == 4:
                    return cls.handle_minus((1 - b)*pi)
                elif oct == 5:
                    return -cls.handle_minus((b % 1)*pi)
                elif oct == 6:
                    return -Cos.handle_minus((3/S(2) - b)*pi)
                elif oct == 7:
                    return -Cos.handle_minus((b % (3/S(2)))*pi)
                elif oct == 8:
                    return -cls.handle_minus((2 - b)*pi)
        else:
            return cls.handle_minus(b*pi)
        """

    @classmethod
    def eval_indirect(cls, a, b, b_mod, oct):
        if oct == 1 or oct == 5:
            return cls.handle_minus(a + b_mod*pi)
        elif oct == 2 or oct == 6:
            return Cot.handle_minus(-a + b_mod*pi)
        elif oct == 3 or oct == 7:
            return -Cot.handle_minus(a + b_mod*pi)
        elif oct == 4 or oct == 8:
            return -cls.handle_minus(-a + b_mod*pi)

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

    @classmethod
    def eval_indirect(cls, a, b, b_mod, oct):
        if oct == 1 or oct == 5:
            return cls.handle_minus(a + b_mod*pi)
        elif oct == 2 or oct == 6:
            return Tan.handle_minus(-a + b_mod*pi)
        elif oct == 3 or oct == 7:
            return -Tan.handle_minus(a + b_mod*pi)
        elif oct == 4 or oct == 8:
            return -cls.handle_minus(-a + b_mod*pi)

def Sec(a):
    pass

def Csc(a):
    pass

def Sinh(a):
    pass

def Cosh(a):
    pass

def Tanh(a):
    pass

def Coth(a):
    pass

def Sech(a):
    pass

def Csch(a):
    pass
# pi/2-x symmetry:
conjugates = {
    Sin: Cos,
    Cos: Sin,
    Tan: Cot,
    Cot: Tan,
    }

def get_pi_shift(arg):
    """
    If arg = a + b*pi, returns (a, b), otherwise None.
    """
    a = Wild("a", exclude=[pi])
    b = Wild("b", exclude=[pi])
    r = arg.match(a+b*pi)
    # I think it should always match:
    if r is None:
        return arg, S(0)
    else:
        return sympify(r[a]), sympify(r[b])



sin = Sin
cos = Cos
tan = Tan
cot = Cot
sec = Sec
csc = Csc
sinh = Sinh
cosh = Cosh
tanh = Tanh
coth = Coth
sech = Sech
csch = Csch



def test_get_pi_shift():
    x, y, n = symbols('x y n')
    assert get_pi_shift(x+2*pi/12) == (x, 1/S(6))
    assert get_pi_shift(y+n*pi) == (y, n)
    assert get_pi_shift(pi/2) == (0, 1/S(2))
    assert get_pi_shift(y) == (y, 0)
    assert get_pi_shift(x + y) == (x+y, 0)

def test_sin():
    x, y, n = symbols('x y n')
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

    assert sin(x - 15*pi/8) == sin(x + pi/8)
    assert sin(x + 3*pi/8) == cos(pi/8 - x)
    assert sin(x - 13*pi/8) == sin(x - 13*pi/8)
    assert sin(x + 5*pi/8) == cos(x + pi/8)
    assert sin(x - 11*pi/8) == cos(x + pi/8)
    assert sin(x + 7*pi/8) == sin(pi/8 - x)
    assert sin(x - 9*pi/8) == sin(pi/8 - x)
    assert sin(x + 9*pi/8) == -sin(x + pi/8)
    assert sin(x - 7*pi/8) == -sin(x + pi/8)
    assert sin(x + 11*pi/8) == -cos(pi/8 - x)
    assert sin(x - 5*pi/8) == -cos(pi/8 - x)
    assert sin(x + 13*pi/8) == -cos(x + pi/8)
    assert sin(x - 3*pi/8) == -cos(x + pi/8)
    assert sin(x + 15*pi/8) == -sin(pi/8 - x)
    assert sin(x - pi/8) == -sin(pi/8 - x)

    x, y = symbols('xy')

    r = Symbol('r', real=True)

    k = Symbol('k', integer=True)

    #assert sin(nan) == nan

    #assert sin(oo*I) == oo*I
    #assert sin(-oo*I) == -oo*I

    assert sin(0) == 0

    assert sin(1) == sin(1)
    assert sin(-1) == -sin(1)

    assert sin(x) == sin(x)
    assert sin(-x) == -sin(x)

    #assert sin(asin(x)) == x
    #assert sin(atan(x)) == x / sqrt(1 + x**2)
    #assert sin(acos(x)) == sqrt(1 - x**2)
    #assert sin(acot(x)) == 1 / (sqrt(1 + 1 / x**2) * x)

    #assert sin(pi*I) == sinh(pi)*I
    #assert sin(-pi*I) == -sinh(pi)*I

    assert sin(2**1024 * E) == sin(2**1024 * E)
    assert sin(-2**1024 * E) == -sin(2**1024 * E)

    assert sin(pi) == 0
    assert sin(-pi) == 0
    assert sin(2*pi) == 0
    assert sin(-2*pi) == 0
    assert sin(-3*10**73*pi) == 0
    assert sin(7*10**103*pi) == 0

    assert sin(pi/2) == 1
    assert sin(-pi/2) == -1
    assert sin(5*pi/2) == 1
    assert sin(7*pi/2) == -1

    assert sin(pi/3) == S.Half*sqrt(3)
    assert sin(-2*pi/3) == -S.Half*sqrt(3)

    assert sin(pi/4) == S.Half*sqrt(2)
    assert sin(-pi/4) == -S.Half*sqrt(2)
    assert sin(17*pi/4) == S.Half*sqrt(2)
    assert sin(-3*pi/4) == -S.Half*sqrt(2)

    assert sin(pi/6) == S.Half
    assert sin(-pi/6) == -S.Half
    assert sin(7*pi/6) == -S.Half
    assert sin(-5*pi/6) == -S.Half

    assert sin(1*pi/5) == sqrt((5 - sqrt(5)) / 8)
    assert sin(2*pi/5) == sqrt((5 + sqrt(5)) / 8)
    assert sin(3*pi/5) == sin(2*pi/5)
    assert sin(4*pi/5) == sin(1*pi/5)
    assert sin(6*pi/5) == -sin(1*pi/5)
    assert sin(8*pi/5) == -sin(2*pi/5)

    assert sin(-1273*pi/5) == -sin(2*pi/5)

    assert sin(pi/105) == sin(pi/105)
    assert sin(-pi/105) == -sin(pi/105)

    assert sin(104*pi/105) == sin(pi/105)
    assert sin(106*pi/105) == -sin(pi/105)

    assert sin(-104*pi/105) == -sin(pi/105)
    assert sin(-106*pi/105) == sin(pi/105)

    assert sin(2 + 3*I) == sin(2 + 3*I)

    #assert sin(x*I) == sinh(x)*I

    assert sin(k*pi) == 0
    assert sin(17*k*pi) == 0

    #assert sin(k*pi*I) == sinh(k*pi)*I

    #assert sin(r).is_real == True

    assert sin(exp(10)-1) == sin(-1+exp(10))

    assert sin(1*pi/9) == sin(pi/9)
    assert sin(2*pi/9) == sin(2*pi/9)
    assert sin(3*pi/9) == sqrt(3)/2
    assert sin(4*pi/9) == cos(pi/18)
    assert sin(5*pi/9) == cos(pi/18)
    assert sin(6*pi/9) == sqrt(3)/2
    assert sin(7*pi/9) == sin(2*pi/9)
    assert sin(8*pi/9) == sin(pi/9)
    assert sin(9*pi/9) == 0
    assert sin(10*pi/9) == -sin(pi/9)
    assert sin(11*pi/9) == -sin(2*pi/9)
    assert sin(12*pi/9) == -sqrt(3)/2
    assert sin(13*pi/9) == -cos(pi/18)
    assert sin(14*pi/9) == -cos(pi/18)
    assert sin(15*pi/9) == -sqrt(3)/2
    assert sin(16*pi/9) == -sin(2*pi/9)
    assert sin(17*pi/9) == -sin(pi/9)
    assert sin(18*pi/9) == 0

def test_cos():
    n, x, y = symbols('n x y')

    r = Symbol('r', real=True)

    k = Symbol('k', integer=True)
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

    assert cos(x - 15*pi/8) == cos(x + pi/8)
    assert cos(x + 3*pi/8) == sin(pi/8 - x)
    assert cos(x - 13*pi/8) == sin(pi/8 - x)
    assert cos(x + 5*pi/8) == -sin(x + pi/8)
    assert cos(x - 11*pi/8) == -sin(x + pi/8)
    assert cos(x + 7*pi/8) == -cos(pi/8 - x)
    assert cos(x - 9*pi/8) == -cos(pi/8 - x)
    assert cos(x + 9*pi/8) == -cos(x + pi/8)
    assert cos(x - 7*pi/8) == -cos(x + pi/8)
    assert cos(x + 11*pi/8) == sin(x - pi/8)
    assert cos(x - 5*pi/8) == sin(x - pi/8)
    assert cos(x + 13*pi/8) == sin(x + pi/8)
    assert cos(x - 3*pi/8) == sin(x + pi/8)
    assert cos(x + 15*pi/8) == cos(pi/8 - x)
    assert cos(x - pi/8) == cos(pi/8 - x)

    #assert cos(nan) == nan

    #assert cos(oo*I) == oo
    #assert cos(-oo*I) == oo

    assert cos(0) == 1

    assert cos(1) == cos(1)
    assert cos(-1) == cos(1)

    assert cos(x) == cos(x)
    assert cos(-x) == cos(x)

    #assert cos(acos(x)) == x
    #assert cos(atan(x)) == 1 / sqrt(1 + x**2)
    #assert cos(asin(x)) == sqrt(1 - x**2)
    #assert cos(acot(x)) == 1 / sqrt(1 + 1 / x**2)

    #assert cos(pi*I) == cosh(pi)
    #assert cos(-pi*I) == cosh(pi)

    assert cos(2**1024 * E) == cos(2**1024 * E)
    assert cos(-2**1024 * E) == cos(2**1024 * E)

    assert cos(pi/2) == 0
    assert cos(-pi/2) == 0
    assert cos(pi/2) == 0
    assert cos(-pi/2) == 0
    assert cos((-3*10**73+1)*pi/2) == 0
    assert cos((7*10**103+1)*pi/2) == 0

    assert cos(pi) == -1
    assert cos(-pi) == -1
    assert cos(2*pi)==1
    assert cos(5*pi) == -1
    assert cos(8*pi) == 1

    assert cos(pi/3) == S.Half
    assert cos(-2*pi/3) == -S.Half

    assert cos(pi/4) == S.Half*sqrt(2)
    assert cos(-pi/4) == S.Half*sqrt(2)
    assert cos(11*pi/4) == -S.Half*sqrt(2)
    assert cos(-3*pi/4) == -S.Half*sqrt(2)

    assert cos(pi/6) == S.Half*sqrt(3)
    assert cos(-pi/6) == S.Half*sqrt(3)
    assert cos(7*pi/6) == -S.Half*sqrt(3)
    assert cos(-5*pi/6) == -S.Half*sqrt(3)

    assert cos(1*pi/5) == (sqrt(5) + 1)/4
    assert cos(2*pi/5) == (sqrt(5) - 1)/4
    assert cos(3*pi/5) == -cos(2*pi/5)
    assert cos(4*pi/5) == -cos(1*pi/5)
    assert cos(6*pi/5) == -cos(1*pi/5)
    assert cos(8*pi/5) == cos(2*pi/5)

    assert cos(-1273*pi/5) == -cos(2*pi/5)

    assert cos(pi/105) == cos(pi/105)
    assert cos(-pi/105) == cos(pi/105)

    assert cos(104*pi/105) == -cos(pi/105)
    assert cos(106*pi/105) == -cos(pi/105)

    assert cos(-104*pi/105) == -cos(pi/105)
    assert cos(-106*pi/105) == -cos(pi/105)

    assert cos(2 + 3*I) == cos(2 + 3*I)

    assert cos(x*I) == cosh(x)

    assert cos(k*pi) == cos(k*pi)
    assert cos(17*k*pi) == cos(17*k*pi)

    assert cos(k*pi*I) == cosh(k*pi)

    assert cos(r).is_real == True

    assert cos(exp(10)-1) == cos(-1+exp(10))


def test_tan():
    var("x y n")
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

    assert tan(x - 15*pi/8) == tan(x + pi/8)
    assert tan(x + 3*pi/8) == cot(pi/8 - x)
    assert tan(x - 13*pi/8) == cot(pi/8 - x)
    assert tan(x + 5*pi/8) == -cot(x + pi/8)
    assert tan(x - 11*pi/8) == -cot(x + pi/8)
    assert tan(x + 7*pi/8) == -tan(pi/8 - x)
    assert tan(x - 9*pi/8) == -tan(pi/8 - x)
    assert tan(x + 9*pi/8) == tan(x + pi/8)
    assert tan(x - 7*pi/8) == tan(x + pi/8)
    assert tan(x + 11*pi/8) == cot(pi/8 - x)
    assert tan(x - 5*pi/8) == cot(pi/8 - x)
    assert tan(x + 13*pi/8) == -cot(x + pi/8)
    assert tan(x - 3*pi/8) == -cot(x + pi/8)
    assert tan(x + 15*pi/8) == -tan(pi/8 - x)
    assert tan(x - pi/8) == -tan(pi/8 - x)


def test_cot():
    var("x y n")
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

    assert cot(x - 15*pi/8) == cot(x + pi/8)
    assert cot(x + 3*pi/8) == tan(pi/8 - x)
    assert cot(x - 13*pi/8) == tan(pi/8 - x)
    assert cot(x + 5*pi/8) == -tan(x + pi/8)
    assert cot(x - 11*pi/8) == -tan(x + pi/8)
    assert cot(x + 7*pi/8) == -cot(pi/8 - x)
    assert cot(x - 9*pi/8) == -cot(pi/8 - x)
    assert cot(x + 9*pi/8) == cot(x + pi/8)
    assert cot(x - 7*pi/8) == cot(x + pi/8)
    assert cot(x + 11*pi/8) == tan(pi/8 - x)
    assert cot(x - 5*pi/8) == tan(pi/8 - x)
    assert cot(x + 13*pi/8) == -tan(x + pi/8)
    assert cot(x - 3*pi/8) == -tan(x + pi/8)
    assert cot(x + 15*pi/8) == -cot(pi/8 - x)
    assert cot(x - pi/8) == -cot(pi/8 - x)

def test_symmetry():
    var("x y n")
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
