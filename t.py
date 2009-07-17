from sympy import (Symbol, pi, S, Basic, Function, solve, latex, sqrt, sympify,
        var, Wild, symbols, floor, Rational, I, E, zoo, oo, exp, nan)
from sympy.core.cache import cacheit
from sympy import sinh, cosh, tanh, coth, asinh, acosh, atanh, acoth
Sinh = sinh
Cosh = cosh
Tanh = tanh
Coth = coth

C0 = (sqrt(3)-1)/(2*sqrt(2))
C1 = S(1)/2
C2 = 1/sqrt(2)
C3 = sqrt(3)/2
C4 = (sqrt(3)+1)/(2*sqrt(2))

# sin_table[n] represents the value of sin(2*pi*n/24) for n = 0..23
sin_table = [
        S.Zero,  C0,  C1,  C2,  C3,  C4,  S.One,  C4,  C3,  C2,  C1,  C0,
        S.Zero, -C0, -C1, -C2, -C3, -C4, S.NegativeOne, -C4, -C3, -C2, -C1, -C0
        ]
# sin_table2[n] represents the value of sin(n*pi/10) for n = 0...19
C02 = (sqrt(5)-1)/4
C12 = sqrt(5/S(8) - sqrt(5)/8)
C22 = (sqrt(5)+1)/4
C32 = sqrt(5/S(8) + sqrt(5)/8)

sin_table2 = [
        S.Zero,  C02,  C12,  C22,  C32,   S.One,  C32,  C22,  C12,  C02,
        S.Zero, -C02, -C12, -C22, -C32,  S.NegativeOne, -C32, -C22, -C12, -C02]

class TrigFunction(Basic):
    """
    Base class for all trigonometric functions.
    """

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
            # Handle inverse trig functions as the last step before evaluation
            if isinstance(-x, (ASin, ACos, ATan, ACot, ASec, ACsc)):
                if cls.odd:
                    return -cls.handle_inverse_trig(-x, type(-x))
                else:
                    return cls.handle_inverse_trig(-x, type(-x))
            elif cls.odd:
                return -cls(-x, eval=False)
            else:
                return cls(-x, eval=False)
        else:
            return cls(x, eval=False)

    @classmethod
    def eval(cls, arg):

        if isinstance(arg, (ASin, ACos, ATan, ACot, ASec, ACsc)):
            return cls.handle_inverse_trig(arg, type(arg))
        # Match arg = a + b*pi
        a, b = get_pi_shift(arg)

        # Case 1:  a == 0 and b == 0
        if a == 0 and b == 0:
            return cls.eval_direct(0)

        # Case 2:  a != 0 and b == 0
        elif a != 0 and b == 0:
            if a.is_imaginary:
                return cls.hyper(a.coeff(I))
            elif a is S.NaN:
                return S.NaN
            else:
                return cls.handle_minus(a)

        # Case 3: a == 0 and b != 0
        elif a == 0 and b != 0:
            if b.is_imaginary:
                return cls.hyper(b.coeff(I)*pi)
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

    @classmethod
    def still_trig_function(cls, expr):
        """
        Used to determine if the result of eval_direct returned an expression
        with a trig function or not.
        """
        if expr.atoms(TrigFunction) == set():
            return False
        else:
            return True

class InverseTrigFunction(Basic):
    """
    Base class for all inverse trig functions.
    """

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
    def eval(cls, arg):
        # Handle inverse trig functions
        if isinstance(arg, (Sin, Cos, Tan, Cot, Sec, Csc)):
            return cls.handle_inverse_trig(arg, type(arg))


class Sin(TrigFunction):
    odd = True
    period = 2

    @classmethod
    def eval_direct(cls, b):
        """
        Returns the value of Sin(b*pi) when direct evaluation is possible.
        """
        b = sympify(b)
        if b.is_rational:
            if b.is_integer:
                return S.Zero
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
            # Case when no direct evaluation is possible
            return cls.handle_minus(b*pi)

    @classmethod
    def eval_indirect(cls, a, b, b_mod, oct):
        """
        Puts any pi-shifts into the interval (0, pi/4)
        """
        # Re-evaluate if the pi shift can be zero using another trig function
        if b_mod == 0:
            if oct == 1:
                return cls.eval(a)
            elif oct == 2:
                return Cos.eval(a)
            elif oct == 3:
                return Cos.eval(a)
            elif oct == 4:
                return -cls.eval(a)
            elif oct == 5:
                return -cls.eval(a)
            elif oct == 6:
                return -Cos.eval(a)
            elif oct == 7:
                return -Cos.eval(a)
            elif oct == 8:
                return cls.eval(a)
        else:
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
    def fdiff(self, argindex=1):
        if argindex == 1:
            return Cos(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def handle_inverse_trig(cls, arg, argtype):
        x = arg.args[0]
        if argtype == ASin:
            return x
        elif argtype == ACos:
            return sqrt(1 - x**2)
        elif argtype == ATan:
            return x / sqrt(1 + x**2)
        elif argtype == ACsc:
            return 1/x
        elif argtype == ASec:
            return sqrt(1 - 1 / x**2)
        elif argtype == ACot:
            return 1 / (sqrt(1 + 1 / x**2) * x)

    @classmethod
    def reciprocal(cls):
        return Csc

    @classmethod
    def hyper(cls, coef):
        """
        Called when argument is imaginary.
        """
        return I*Sinh(coef)

class Cos(TrigFunction):
    odd = False
    period = 2

    @classmethod
    def eval_direct(cls, b):
        """
        Returns the value of Cos(b*pi) when direct evaluation is possible.
        """
        return Sin.eval_direct(b + 1/S(2))

    @classmethod
    def eval_indirect(cls, a, b, b_mod, oct):
        """
        Puts any pi-shifts into the interval (0, pi/4)
        """
        return Sin.eval(a + (b + 1/S(2))*pi)

    def as_Sin(self):
            return Sin(pi/2 - self.args[0])

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -Sin(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def handle_inverse_trig(cls, arg, argtype):
        x = arg.args[0]
        if argtype == ASin:
            return sqrt(1 - x**2)
        elif argtype == ACos:
            return x
        elif argtype == ATan:
            return 1 / sqrt(1 + x**2)
        elif argtype == ACsc:
            return sqrt(1 - 1 / x**2)
        elif argtype == ASec:
            return 1 / x
        elif argtype == ACot:
            return 1 / (sqrt(1 + 1 / x**2))

    @classmethod
    def reciprocal(cls):
        return Sec

    @classmethod
    def hyper(cls, coef):
        """
        Called when argument is imaginary.
        """
        return Cosh(coef)

class Tan(TrigFunction):
    odd = True
    period = 1

    @classmethod
    def eval_direct(cls, b):
        """
        Returns the value of Tan(b*pi) when direct evaluation is possible.
        """
        num = Sin.eval_direct(b)
        den = Cos.eval_direct(b)
        ln = len(num.args)
        ld = len(den.args)
        if num.atoms(Sin):
            tf = Tan
        else:
            tf = Cot

        if den == S.Zero:
            return zoo
        else:
            if all(map(cls.still_trig_function, (num, den))):
                if ln == 1 and ld == 1:
                    assert num.args == den.args
                    return tf.handle_minus(num.args[0])
                elif ln == 1 and ld == 2:
                    assert num.args[0] == den.args[1].args[0]
                    return tf.handle_minus(num.args[0]) / den.args[0]
                elif ln == 2 and ld == 1:
                    assert num.args[1].args[0] == den.args[0]
                    return num.args[0] * tf.handle_minus(num.args[1].args[0])
                elif ln == 2 and ld == 2:
                    assert num.args[1].args[0] == den.args[1].args[0]
                    return num.args[0] * tf.handle_minus(num.args[1].args[0]) /\
                        den.args[0]
            else:
                return num / den

    @classmethod
    def eval_indirect(cls, a, b, b_mod, oct):
        """
        Puts any pi-shifts into the interval (0, pi/4)
        """
        if b_mod == 0:
            if oct == 1 or oct == 5:
                return cls.eval(a)
            elif oct == 2 or oct == 6:
                return -Cot.eval(a)
            elif oct == 3 or oct == 7:
                return -Cot.eval(a)
            elif oct == 4 or oct == 8:
                return cls.eval(a)
        else:
            if oct == 1 or oct == 5:
                return cls.handle_minus(a + b_mod*pi)
            elif oct == 2 or oct == 6:
                return Cot.handle_minus(-a + b_mod*pi)
            elif oct == 3 or oct == 7:
                return -Cot.handle_minus(a + b_mod*pi)
            elif oct == 4 or oct == 8:
                return -cls.handle_minus(-a + b_mod*pi)

    def fdiff(self, argindex=1):
        if argindex == 1:
            return Sec(self.args[0])**2
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def handle_inverse_trig(cls, arg, argtype):
        x = arg.args[0]
        if argtype == ASin:
            return x / sqrt(1 - x**2)
        elif argtype == ACos:
            return sqrt(1 - x**2) / x
        elif argtype == ATan:
            return x
        elif argtype == ACsc:
            return 1 / (sqrt(1 - 1 / x**2) * x)
        elif argtype == ASec:
            return sqrt(1 - 1 / x**2) * x
        elif argtype == ACot:
            return 1 / x

    @classmethod
    def reciprocal(cls):
        return Cot

    @classmethod
    def hyper(cls, coef):
        """
        Called when argument is imaginary.
        """
        return I*Tanh(coef)

class Csc(TrigFunction):
    odd = True
    period = 2

    @classmethod
    def eval_direct(cls, m):
        """
        Returns the value of Csc(b*pi).
        """
        den = Sin.eval_direct(m)
        if den == 0:
            return zoo
        ld = len(den.args)
        if cls.still_trig_function(den):
            if ld == 1:
                return den.reciprocal.handle_minus(den.args[0])
            else:
                return den.args[1].reciprocal.handle_minus(den.args[1].args) / \
                        den.args[0]
        else:
            return 1 / den

    @classmethod
    def eval_indirect(cls, a, b, b_mod, oct):
        """
        Puts any pi-shifts into the interval (0, pi/4)
        """
        den = Sin.eval_indirect(m)
        ld = len(den.args)
        if ld == 1:
            return den.reciprocal.handle_minus(den.args[0])
        else:
            return den.args[1].reciprocal.handle_minus(den.args[1])\
                    / den.args[0]

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -Cot(self.args[0])*Csc(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def handle_inverse_trig(cls, arg, argtype):
        x = arg.args[0]
        if argtype == ASin:
            return 1 / x
        elif argtype == ACos:
            return 1 / sqrt(1 - x**2)
        elif argtype == ATan:
            return sqrt(1 + x**2) / x
        elif argtype == ACsc:
            return x
        elif argtype == ASec:
            return sqrt(1 - 1 / x**2)
        elif argtype == ACot:
            return sqrt(1 + 1 / x**2)

    @classmethod
    def reciprocal(cls):
        return Sin

    @classmethod
    def hyper(cls, coef):
        """
        Called when argument is imaginary.
        """
        return -I*Csch(coef)

class Sec(TrigFunction):
    odd = False
    period = 2

    @classmethod
    def eval_direct(cls, b):
        """
        Returns the value of Sec(b*pi).
        """
        den = Cos.eval_direct(m)
        if den == 0:
            return zoo
        ld = len(den.args)
        if cls.still_trig_function(den):
            if ld == 1:
                return den.reciprocal.handle_minus(den.args[0])
            else:
                return den.args[1].reciprocal.handle_minus(den.args[1].args) \
                        / den.args[0]
        else:
            return 1 / den

    @classmethod
    def eval_indirect(cls, a, b, b_mod, oct):
        """
        Puts any pi-shifts into the interval (0, pi/4)
        """
        den = Cos.eval_indirect(m)
        ld = len(den.args)
        if ld == 1:
            return den.reciprocal.handle_minus(den.args[0])
        else:
            return den.args[1].reciprocal.handle_minus(den.args[1]) / \
                    den.args[0]

    def fdiff(self, argindex=1):
        if argindex == 1:
            return Tan(self.args[0])*Sec(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def handle_inverse_trig(cls, arg, argtype):
        x = arg.args[0]
        if argtype == ASin:
            return 1 / sqrt(1 - x**2)
        elif argtype == ACos:
            return 1 / x
        elif argtype == ATan:
            return sqrt(1 + x**2)
        elif argtype == ACsc:
            return 1 / sqrt(1 - 1 / x**2)
        elif argtype == ASec:
            return x
        elif argtype == ACot:
            return sqrt(1 + 1 / x**2)

    @classmethod
    def reciprocal(cls):
        return Cos

    @classmethod
    def hyper(cls, coef):
        """
        Called when argument is imaginary.
        """
        return Sech(coef)

class Cot(TrigFunction):
    odd = True
    period = 1

    @classmethod
    def eval_direct(cls, b):
        """
        Returns the value of Cot(b*pi) when direct evaluation is possible.
        """
        num = Cos.eval_direct(b)
        den = Sin.eval_direct(b)
        # Sometimes eval_direct will return a -Sin/-Cos ...
        ln = len(num.args)
        ld = len(den.args)
        # Check if the numerator has a Sin
        if num.atoms(Sin):
            tf = Tan
        else:
            tf = Cot

        if den == S.Zero:
            return zoo
        else:
            if all(map(cls.still_trig_function, (num, den))):
                if ln == 1 and ld == 1:
                    assert num.args == den.args
                    return tf.handle_minus(num.args[0])
                elif ln == 1 and ld == 2:
                    assert num.args[0] == den.args[1].args[0]
                    return tf.handle_minus(num.args[0]) / den.args[0]
                elif ln == 2 and ld == 1:
                    assert num.args[1].args[0] == den.args[0]
                    return num.args[0] * tf.handle_minus(num.args[1].args[0])
                elif ln == 2 and ld == 2:
                    assert num.args[1].args[0] == den.args[1].args[0]
                    return num.args[0] * tf.handle_minus(num.args[1].args[0]) /\
                        den.args[0]
            else:
                return num / den

    @classmethod
    def eval_indirect(cls, a, b, b_mod, oct):
        """
        Puts any pi-shifts into the interval (0, pi/4)
        """
        if oct == 1 or oct == 5:
            return cls.handle_minus(a + b_mod*pi)
        elif oct == 2 or oct == 6:
            return Tan.handle_minus(-a + b_mod*pi)
        elif oct == 3 or oct == 7:
            return -Tan.handle_minus(a + b_mod*pi)
        elif oct == 4 or oct == 8:
            return -cls.handle_minus(-a + b_mod*pi)

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -Csc(self.args[0])**2
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def handle_inverse_trig(cls, arg, argtype):
        x = arg.args[0]
        if argtype == ASin:
            return sqrt(1 - x**2) / x
        elif argtype == ACos:
            return x / sqrt(1 - x**2)
        elif argtype == ATan:
            return 1 / x
        elif argtype == ACsc:
            return sqrt(1 - 1 / x**2) * x
        elif argtype == ASec:
            return 1 / (sqrt(1 - 1 / x**2) * x)
        elif argtype == ACot:
            return x

    @classmethod
    def reciprocal(cls):
        return Tan

    @classmethod
    def hyper(cls, coef):
        """
        Called when argument is imaginary.
        """
        return -I*Coth(coef)

"""
def Sinh(a):
    pass

def Cosh(a):
    pass

def Tanh(a):
    pass

def Csch(a):
    pass

def Sech(a):
    pass

def Coth(a):
    pass

def ASinh(a):
    pass

def ACosh(a):
    pass

def ATanh(a):
    pass

def ACsch(a):
    pass

def ASech(a):
    pass

def ACoth(a):
    pass
"""

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
        return sympify(r.get(a, 0)), sympify(r.get(b, 0))

def get_imag_coef(arg):
    """
    If arg = a + b*I, returns (a, b), otherwise None.
    """
    a = Wild("a", exclude=[I])
    b = Wild("b", exclude=[I])
    r = arg.match(a+b*I)
    # I think it should always match:
    if r is None:
        return arg, S(0)
    else:
        return sympify(r.get(a, 0)), sympify(r.get(b, 0))





###############################################################################
########################### TRIGONOMETRIC INVERSES ############################
###############################################################################

class ASin(InverseTrigFunction):
    """
    Usage
    =====
      asin(x) -> Returns the arc sine of x (measured in radians)
    """

    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return sqrt(1 - self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def _eval_apply_subs(self, *args):
        return

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.NegativeInfinity * S.ImaginaryUnit
            elif arg is S.NegativeInfinity:
                return S.Infinity * S.ImaginaryUnit
            elif arg is S.Zero:
                return S.Zero
            elif arg is S.One:
                return S.Pi / 2
            elif arg is S.NegativeOne:
                return -S.Pi / 2

        if arg.is_number:
            cst_table = {
                S.Half     : 6,
                -S.Half    : -6,
                sqrt(2)/2  : 4,
                -sqrt(2)/2 : -4,
                1/sqrt(2)  : 4,
                -1/sqrt(2) : -4,
                sqrt(3)/2  : 3,
                -sqrt(3)/2 : -3,
                }

            if arg in cst_table:
                return S.Pi / cst_table[arg]
            elif arg.is_negative:
                return -cls(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * C.asinh(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -cls(-arg)


    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return p * (n-2)**2/(k*(k-1)) * x**2
            else:
                k = (n - 1) // 2

                R = C.RisingFactorial(S.Half, k)
                F = C.Factorial(k)

                return R / F * x**n / n

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if C.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real and (self.args[0]>=-1 and self.args[0]<=1)

    def _sage_(self):
        import sage.all as sage
        return sage.asin(self.args[0]._sage_())

class ACos(InverseTrigFunction):
    """
    Usage
    =====
      acos(x) -> Returns the arc cosine of x (measured in radians)
    """

    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -sqrt(1 - self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def _eval_apply_subs(self, *args):
        return

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Infinity * S.ImaginaryUnit
            elif arg is S.NegativeInfinity:
                return S.NegativeInfinity * S.ImaginaryUnit
            elif arg is S.Zero:
                return S.Pi / 2
            elif arg is S.One:
                return S.Zero
            elif arg is S.NegativeOne:
                return S.Pi

        if arg.is_number:
            cst_table = {
                S.Half     : S.Pi/3,
                -S.Half    : 2*S.Pi/3,
                sqrt(2)/2  : S.Pi/4,
                -sqrt(2)/2 : 3*S.Pi/4,
                1/sqrt(2)  : S.Pi/4,
                -1/sqrt(2) : 3*S.Pi/4,
                sqrt(3)/2  : S.Pi/6,
                -sqrt(3)/2 : 5*S.Pi/6,
                }

            if arg in cst_table:
                return cst_table[arg]


    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n == 0:
            return S.Pi / 2
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return p * (n-2)**2/(k*(k-1)) * x**2
            else:
                k = (n - 1) // 2

                R = C.RisingFactorial(S.Half, k)
                F = C.Factorial(k)

                return -R / F * x**n / n

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if C.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real and (self.args[0]>=-1 and self.args[0]<=1)

    def _sage_(self):
        import sage.all as sage
        return sage.acos(self.args[0]._sage_())

class ATan(InverseTrigFunction):
    """
    Usage
    =====
      atan(x) -> Returns the arc tangent of x (measured in radians)
    """

    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/(1+self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def _eval_apply_subs(self, *args):
        return

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Pi / 2
            elif arg is S.NegativeInfinity:
                return -S.Pi / 2
            elif arg is S.Zero:
                return S.Zero
            elif arg is S.One:
                return S.Pi / 4
            elif arg is S.NegativeOne:
                return -S.Pi / 4

        if arg.is_number:
            cst_table = {
                sqrt(3)/3  : 6,
                -sqrt(3)/3 : -6,
                1/sqrt(3)  : 6,
                -1/sqrt(3) : -6,
                sqrt(3)    : 3,
                -sqrt(3)   : -3,
                }

            if arg in cst_table:
                return S.Pi / cst_table[arg]
            elif arg.is_negative:
                return -cls(-arg)

        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * C.atanh(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -cls(-arg)


    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)
            return (-1)**((n-1)//2) * x**n / n

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if C.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real

    def _sage_(self):
        import sage.all as sage
        return sage.atan(self.args[0]._sage_())

class ACsc(InverseTrigFunction):
    """
    Usage
    =====
      ACsc(x) -> Returns the arc cosecant of x (measured in radians)
    """

    nargs = 1
    # done
    def fdiff(self, argindex=1):
        if argindex == 1:
            return -1 / (sqrt(1 - 1 / x**2) * x**2)
        else:
            raise ArgumentIndexError(self, argindex)
    # done
    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Pi / 2
            elif arg is S.NegativeInfinity:
                return S.Pi / 2
            elif arg is S.Zero:
                return zoo
            elif arg is S.One:
                return S.Zero
            elif arg is S.NegativeOne:
                return S.Pi

        if arg.is_number:
            cst_table = {
                sqrt(3)/3  : 3,
                -sqrt(3)/3 : -3,
                1/sqrt(3)  : 3,
                -1/sqrt(3) : -3,
                sqrt(3)    : 6,
                -sqrt(3)   : -6,
                }

            if arg in cst_table:
                return S.Pi / cst_table[arg]
            elif arg.is_negative:
                return -cls(-arg)

        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return -S.ImaginaryUnit * C.acoth(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -cls(-arg)


    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n == 0:
            return S.Pi / 2 # FIX THIS
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)
            return (-1)**((n+1)//2) * x**n / n

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if C.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real

    def _sage_(self):
        import sage.all as sage
        return sage.acot(self.args[0]._sage_())

class ASec(InverseTrigFunction):
    """
    Usage
    =====
      ASec(x) -> Returns the arc secant of x (measured in radians)
    """

    nargs = 1
    # done
    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1 / (sqrt(1 - 1 / x**2) * x**2)
        else:
            raise ArgumentIndexError(self, argindex)
    # done
    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Pi / 2
            elif arg is S.NegativeInfinity:
                return S.Pi / 2
            elif arg is S.Zero:
                return zoo
            elif arg is S.One:
                return S.Zero
            elif arg is S.NegativeOne:
                return S.Pi

        if arg.is_number:
            cst_table = {
                sqrt(3)/3  : 3,
                -sqrt(3)/3 : -3,
                1/sqrt(3)  : 3,
                -1/sqrt(3) : -3,
                sqrt(3)    : 6,
                -sqrt(3)   : -6,
                }

            if arg in cst_table:
                return S.Pi / cst_table[arg]
            elif arg.is_negative:
                return -cls(-arg)

        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return -S.ImaginaryUnit * C.acoth(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -cls(-arg)


    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n == 0:
            return S.Pi / 2 # FIX THIS
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)
            return (-1)**((n+1)//2) * x**n / n

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if C.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real

    def _sage_(self):
        import sage.all as sage
        return sage.acot(self.args[0]._sage_())

class ACot(InverseTrigFunction):
    """
    Usage
    =====
      acot(x) -> Returns the arc cotangent of x (measured in radians)
    """

    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -1 / (1+self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Zero
            elif arg is S.NegativeInfinity:
                return S.Zero
            elif arg is S.Zero:
                return S.Pi/ 2
            elif arg is S.One:
                return S.Pi / 4
            elif arg is S.NegativeOne:
                return -S.Pi / 4

        if arg.is_number:
            cst_table = {
                sqrt(3)/3  : 3,
                -sqrt(3)/3 : -3,
                1/sqrt(3)  : 3,
                -1/sqrt(3) : -3,
                sqrt(3)    : 6,
                -sqrt(3)   : -6,
                }

            if arg in cst_table:
                return S.Pi / cst_table[arg]
            elif arg.is_negative:
                return -cls(-arg)

        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return -S.ImaginaryUnit * C.acoth(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -cls(-arg)


    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n == 0:
            return S.Pi / 2 # FIX THIS
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)
            return (-1)**((n+1)//2) * x**n / n

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if C.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real

    def _sage_(self):
        import sage.all as sage
        return sage.acot(self.args[0]._sage_())



sin = Sin
cos = Cos
tan = Tan
cot = Cot
sec = Sec
csc = Csc
asin = ASin
acos = ACos
atan = ATan
acsc = ACsc
asec = ASec
acot = ACot
#sinh = Sinh
#cosh = Cosh
#tanh = Tanh
#coth = Coth
#sech = Sech
#acsch = ACsch
#asinh = ASinh
#acosh = ACosh
#atanh = ATanh
#acoth = ACoth
#asech = ASech
#acsch = ACsch

###### Tests #####
def test_get_pi_shift():
    x, y, n = symbols('x y n')
    assert get_pi_shift(x+2*pi/12) == (x, 1/S(6))
    assert get_pi_shift(y+n*pi) == (y, n)
    assert get_pi_shift(pi/2) == (0, 1/S(2))
    assert get_pi_shift(y) == (y, 0)
    assert get_pi_shift(x + y) == (x+y, 0)


# Need to test 4 cases of for all trig functions, e.g. sin(a + b*pi):
# Case 1) a == 0 and b == 0
def test_sin_case1():
    assert sin(0) == 0

# Case 2) a != 0 and b == 0
def test_sin_case2():
    x, y, n = symbols('x y n')
    r = Symbol('r', real=True)
    k = Symbol('k', integer=True)
    #assert sin(r).is_real == True

    assert sin(exp(10)-1) == sin(-1+exp(10))
    assert sin(S.NaN) == S.NaN
    # See Issue 1540 for why these two fail:
    #assert sin(oo*I) == oo*I
    #assert sin(-oo*I) == -oo*I

    assert sin(x) == sin(x)
    assert sin(-x) == -sin(x)

    assert sin(-y) == -sin(y)

    assert sin(asin(x)) == x
    assert sin(acos(x)) == sqrt(1 - x**2)
    assert sin(atan(x)) == x/sqrt(1 + x**2)
    assert sin(acsc(x)) == 1/x
    assert sin(asec(x)) == sqrt(1 - 1/x**2)
    assert sin(acot(x)) == 1 / (x*sqrt(1 + 1 / x**2))
    assert sin(-asin(x)) == -x
    assert sin(asin(-x)) == -x
    assert sin(-acos(x)) == -sqrt(1 - x**2)
    assert sin(-atan(x)) == -x/sqrt(1 + x**2)
    assert sin(-acsc(x)) == -1/x
    assert sin(-asec(x)) == -sqrt(1 - 1/x**2)
    assert sin(-acot(x)) == -1 / (x*sqrt(1 + 1 / x**2))

    assert sin(pi*I) == sinh(pi)*I
    assert sin(-pi*I) == -sinh(pi)*I

    assert sin(2**1024 * E) == sin(2**1024 * E)
    assert sin(-2**1024 * E) == -sin(2**1024 * E)


    assert sin(2 + 3*I) == sin(2 + 3*I)
    x = Symbol('x', real=True)
    assert sin(x*I) == sinh(x)*I

# Case 3) a == 0 and b != 0
def test_sin_case3():
    x, y, n = symbols('x y n')
    r = Symbol('r', real=True)
    k = Symbol('k', integer=True)
    assert sin(k*pi) == 0
    assert sin(17*k*pi) == 0

    assert sin(k*pi*I) == sinh(k*pi)*I
    assert sin(pi/6) == S(1)/2
    assert sin(pi/4) == 1/sqrt(2)
    assert sin(pi/3) == sqrt(3)/2
    assert sin(pi/2) == 1
    assert sin(0 + 4*pi) == 0
    assert sin(pi/6 + 4*pi) == S(1)/2
    assert sin(pi/4 + 4*pi) == 1/sqrt(2)
    assert sin(pi/3 + 4*pi) == sqrt(3)/2
    assert sin(pi/2 + 4*pi) == 1

    assert sin(0 + 2*pi) == 0
    assert sin(pi/6 + 2*pi) == S(1)/2
    assert sin(pi/4 + 2*pi) == 1/sqrt(2)
    assert sin(pi/3 + 2*pi) == sqrt(3)/2
    assert sin(pi/2 + 2*pi) == 1
    assert sin(3*pi/2) == -1
    assert sin(5*pi/2) == 1
    assert sin(1) == sin(1)
    assert sin(-1) == -sin(1)
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

    assert sin(pi/10) == (sqrt(5)-1)/4
    assert sin(3*pi/10) == (sqrt(5)+1)/4
    assert sin(7*pi/10) == (sqrt(5)+1)/4
    assert sin(9*pi/10) == (sqrt(5)-1)/4
    assert sin(13*pi/10) == (-1-sqrt(5))/4
    assert sin(17*pi/10) == (-1-sqrt(5))/4
    assert sin(19*pi/10) == (1-sqrt(5))/4
    assert sin(-pi/10) == (-sqrt(5)+1)/4
    assert sin(-3*pi/10) == -(sqrt(5)+1)/4
    assert sin(-7*pi/10) == -(sqrt(5)+1)/4
    assert sin(-9*pi/10) == (-sqrt(5)+1)/4
    assert sin(-13*pi/10) == (1+sqrt(5))/4
    assert sin(-17*pi/10) == (1+sqrt(5))/4
    assert sin(-19*pi/10) == (-1+sqrt(5))/4

    assert sin(-1273*pi/5) == -sin(2*pi/5)

    assert sin(pi/105) == sin(pi/105)
    assert sin(-pi/105) == -sin(pi/105)

    assert sin(104*pi/105) == sin(pi/105)
    assert sin(106*pi/105) == -sin(pi/105)

    assert sin(-104*pi/105) == -sin(pi/105)
    assert sin(-106*pi/105) == sin(pi/105)

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

    assert sin(-1*pi/9) == -sin(pi/9)
    assert sin(-2*pi/9) == -sin(2*pi/9)
    assert sin(-3*pi/9) == -sqrt(3)/2
    assert sin(-4*pi/9) == -cos(pi/18)
    assert sin(-5*pi/9) == -cos(pi/18)
    assert sin(-6*pi/9) == -sqrt(3)/2
    assert sin(-7*pi/9) == -sin(2*pi/9)
    assert sin(-8*pi/9) == -sin(pi/9)
    assert sin(-9*pi/9) == 0
    assert sin(-10*pi/9) == sin(pi/9)
    assert sin(-11*pi/9) == sin(2*pi/9)
    assert sin(-12*pi/9) == sqrt(3)/2
    assert sin(-13*pi/9) == cos(pi/18)
    assert sin(-14*pi/9) == cos(pi/18)
    assert sin(-15*pi/9) == sqrt(3)/2
    assert sin(-16*pi/9) == sin(2*pi/9)
    assert sin(-17*pi/9) == sin(pi/9)
    assert sin(-18*pi/9) == 0

# Case 4) a != 0 and b != 0
def test_sin_case4():
    x, y, n = symbols('x y n')
    r = Symbol('r', real=True)
    k = Symbol('k', integer=True)
    assert sin(-y + 2*pi) == -sin(y)
    assert sin(pi - y + 2*pi) == sin(y)
    assert sin(pi + y + 2*pi) == -sin(y)
    assert sin(2*pi - y + 2*pi) == -sin(y)
    assert sin(pi/2 + y + 2*pi) == cos(y)
    assert sin(pi/2 - y + 2*pi) == cos(y)

    assert sin(-y + 4*pi) == -sin(y)
    assert sin(pi - y + 4*pi) == sin(y)
    assert sin(pi + y + 4*pi) == -sin(y)
    assert sin(2*pi - y + 4*pi) == -sin(y)
    assert sin(pi/2 + y + 4*pi) == cos(y)
    assert sin(pi/2 - y + 4*pi) == cos(y)

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
    assert sin(pi - y) == sin(y)
    assert sin(pi + y) == -sin(y)
    assert sin(2*pi - y) == -sin(y)
    assert sin(pi/2 + y) == cos(y)
    assert sin(pi/2 - y) == cos(y)

# Case 1) a == 0 and b == 0
def test_cos_case1():
    assert cos(0) == 1

# Case 2) a != 0 and b == 0
def test_cos_case2():
    n, x, y = symbols('n x y')
    r = Symbol('r', real=True)
    k = Symbol('k', integer=True)
    assert cos(2 + 3*I) == cos(2 + 3*I)
    #assert cos(r).is_real == True
    assert cos(exp(10)-1) == cos(-1+exp(10))
    # Issue 1540 causes these two to fail:
    #assert cos(oo*I) == oo
    #assert cos(-oo*I) == oo

    assert cos(asin(x)) == sqrt(1 - x**2)
    assert cos(acos(x)) == x
    assert cos(atan(x)) == 1 / sqrt(1 + x**2)
    assert cos(acsc(x)) == sqrt(1 - 1 / x**2)
    assert cos(asec(x)) == 1 / x
    assert cos(acot(x)) == 1 / sqrt(1 + 1 / x**2)
    assert cos(-asin(x)) == sqrt(1 - x**2)
    assert cos(asin(-x)) == sqrt(1 - x**2)
    assert cos(-asin(x)) == sqrt(1 - x**2)
    assert cos(-acos(x)) == x
    assert cos(-atan(x)) == 1 / sqrt(1 + x**2)
    assert cos(-acsc(x)) == sqrt(1 - 1 / x**2)
    assert cos(-asec(x)) == 1 / x
    assert cos(-acot(x)) == 1 / sqrt(1 + 1 / x**2)

    assert cos(pi*I) == cosh(pi)
    assert cos(-pi*I) == cosh(pi)

    assert cos(2**1024 * E) == cos(2**1024 * E)
    assert cos(-2**1024 * E) == cos(2**1024 * E)
    assert cos(1) == cos(1)
    assert cos(-y) == cos(y)
    assert cos(nan) == nan
    assert cos(-1) == cos(1)
    assert cos(x) == cos(x)
    assert cos(-x) == cos(x)
    x = Symbol('x', real=True)
    assert cos(x*I) == cosh(x)

# Case 3) a == 0 and b == 0
def test_cos_case3():
    n, x, y = symbols('n x y')
    r = Symbol('r', real=True)
    k = Symbol('k', integer=True)
    assert cos(k*pi) == cos(k*pi)
    assert cos(17*k*pi) == cos(17*k*pi)

    assert cos(k*pi*I) == cosh(k*pi)
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

    assert cos(pi/10) == sqrt(5/S(8) + sqrt(5)/8)
    assert cos(3*pi/10) == sqrt(5/S(8) - sqrt(5)/8)
    assert cos(7*pi/10) == -sqrt(5/S(8) - sqrt(5)/8)
    assert cos(9*pi/10) == -sqrt(5/S(8) + sqrt(5)/8)
    assert cos(11*pi/10) == -sqrt(5/S(8) + sqrt(5)/8)
    assert cos(13*pi/10) == -sqrt(5/S(8) - sqrt(5)/8)
    assert cos(17*pi/10) == sqrt(5/S(8) - sqrt(5)/8)
    assert cos(-pi/10) == sqrt(5/S(8) + sqrt(5)/8)
    assert cos(-3*pi/10) == sqrt(5/S(8) - sqrt(5)/8)
    assert cos(-7*pi/10) == -sqrt(5/S(8) - sqrt(5)/8)
    assert cos(-9*pi/10) == -sqrt(5/S(8) + sqrt(5)/8)
    assert cos(-11*pi/10) == -sqrt(5/S(8) + sqrt(5)/8)
    assert cos(-13*pi/10) == -sqrt(5/S(8) - sqrt(5)/8)
    assert cos(-17*pi/10) == sqrt(5/S(8) - sqrt(5)/8)
    assert cos(-1273*pi/5) == -cos(2*pi/5)

    assert cos(pi/105) == cos(pi/105)
    assert cos(-pi/105) == cos(pi/105)

    assert cos(104*pi/105) == -cos(pi/105)
    assert cos(106*pi/105) == -cos(pi/105)

    assert cos(-104*pi/105) == -cos(pi/105)
    assert cos(-106*pi/105) == -cos(pi/105)
    assert cos(pi) == -1
    assert cos(8*pi) == 1
    assert cos(-9*pi) == -1
    assert cos(3*pi/2) == 0
    assert cos(11*pi/2) == 0
    assert cos(pi/12) == (1 + sqrt(3)) / (2 * sqrt(2))
    assert cos(pi/6) == sqrt(3)/2
    assert cos(pi/4) == 1/sqrt(2)
    assert cos(pi/3) == 1/S(2)
    assert cos(pi/2) == 0
    assert cos(0 + 2*pi) == 1
    assert cos(pi/6 + 2*pi) == sqrt(3)/2
    assert cos(pi/4 + 2*pi) == 1/sqrt(2)
    assert cos(pi/3 + 2*pi) == 1/S(2)
    assert cos(pi/2 + 2*pi) == 0

# Case 4) a != 0 and b != 0
def test_cos_case4():
    n, x, y = symbols('n x y')
    r = Symbol('r', real=True)
    k = Symbol('k', integer=True)
    assert cos(pi - y) == -cos(y)
    assert cos(pi + y) == -cos(y)
    assert cos(2*pi - y) == cos(y)
    assert cos(pi/2 + y) == -sin(y)
    assert cos(pi/2 - y) == sin(y)
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
    assert cos(-y + 2*pi) == cos(y)
    assert cos(pi - y + 2*pi) == -cos(y)
    assert cos(pi + y + 2*pi) == -cos(y)
    assert cos(2*pi - y + 2*pi) == cos(y)
    assert cos(pi/2 + y + 2*pi) == -sin(y)
    assert cos(pi/2 - y + 2*pi) == sin(y)

# Case 1) a == 0 and b == 0
def test_tan_case1():
    assert tan(0) == 0

# Case 2) a != 0 and b == 0
def test_tan_case2():
    var("x y n")
    x, y = symbols('xy')
    r = Symbol('r', real=True)
    k = Symbol('k', integer=True)
    assert tan(2**1024 * E) == tan(2**1024 * E)
    assert tan(-2**1024 * E) == -tan(2**1024 * E)
    assert tan(2 + 3*I) == tan(2 + 3*I)
    #assert tan(r).is_real == True
    assert tan(nan) == nan

    # Bug in match is causing this to fail
    # see Issue 1540
    #assert tan(oo*I) == I
    #assert tan(-oo*I) == -I

    assert tan(1) == tan(1)
    assert tan(-1) == -tan(1)

    assert tan(x) == tan(x)
    assert tan(-x) == -tan(x)

    assert tan(asin(x)) == x / sqrt(1 - x**2)
    assert tan(acos(x)) == sqrt(1 - x**2) / x
    assert tan(atan(x)) == x
    assert tan(acsc(x)) == 1 / (sqrt(1- 1 / x**2) * x)
    assert tan(asec(x)) == sqrt(1- 1 / x**2) * x
    assert tan(acot(x)) == 1 / x
    assert tan(-asin(x)) == -x / sqrt(1 - x**2)
    assert tan(asin(-x)) == -x / sqrt(1 - x**2)
    assert tan(-acos(x)) == -sqrt(1 - x**2) / x
    assert tan(-atan(x)) == -x
    assert tan(-acsc(x)) == -1 / (sqrt(1- 1 / x**2) * x)
    assert tan(-asec(x)) == -sqrt(1- 1 / x**2) * x
    assert tan(-acot(x)) == -1 / x
    assert tan(-y) == -tan(y)

    x = Symbol('x', real=True)
    assert tan(x*I) == tanh(x)*I

# Case 3) a == 0 and b != 0
def test_tan_case3():
    x, y = symbols('xy')
    r = Symbol('r', real=True)
    k = Symbol('k', integer=True)
    assert tan(k*pi) == 0
    assert tan(pi*I) == tanh(pi)*I
    assert tan(-pi*I) == -tanh(pi)*I
    assert tan(17*k*pi) == 0
    assert tan(k*pi*I) == tanh(k*pi)*I
    assert tan(pi/6) == 1/sqrt(3)
    assert tan(pi/4) == 1
    assert tan(pi/3) == sqrt(3)
    assert tan(pi/2) == zoo
    assert tan(0 + pi) == 0
    assert tan(pi/6 + pi) == 1/sqrt(3)
    assert tan(pi/4 + pi) == 1
    assert tan(pi/3 + pi) == sqrt(3)
    assert tan(7*pi/12) == sin(7*pi/12)/cos(7*pi/12)
    assert tan(-16*pi/16) == 0
    assert tan(-15*pi/16) == tan(pi/16)
    assert tan(-14*pi/16) == tan(pi/8)
    assert tan(-13*pi/16) == tan(3*pi/16)
    assert tan(-12*pi/16) == 1
    assert tan(-11*pi/16) == cot(3*pi/16)
    assert tan(-10*pi/16) == cot(pi/8)
    assert tan(-9*pi/16) == cot(pi/16)
    assert tan(-8*pi/16) == zoo
    assert tan(-7*pi/16) == -cot(pi/16)
    assert tan(-6*pi/16) == -cot(pi/8)
    assert tan(-5*pi/16) == -cot(3*pi/16)
    assert tan(-4*pi/16) == -1
    assert tan(-3*pi/16) == -tan(3*pi/16)
    assert tan(-2*pi/16) == -tan(pi/8)
    assert tan(-1*pi/16) == -tan(pi/16)
    assert tan(0*pi/16) == 0
    assert tan(1*pi/16) == tan(pi/16)
    assert tan(2*pi/16) == tan(pi/8)
    assert tan(3*pi/16) == tan(3*pi/16)
    assert tan(4*pi/16) == 1
    assert tan(5*pi/16) == cot(3*pi/16)
    assert tan(6*pi/16) == cot(pi/8)
    assert tan(7*pi/16) == cot(pi/16)
    assert tan(8*pi/16) == zoo
    assert tan(9*pi/16) == -cot(pi/16)
    assert tan(10*pi/16) == -cot(pi/8)
    assert tan(11*pi/16) == -cot(3*pi/16)
    assert tan(12*pi/16) == -1
    assert tan(13*pi/16) == -tan(3*pi/16)
    assert tan(14*pi/16) == -tan(pi/8)
    assert tan(15*pi/16) == -tan(pi/16)
    assert tan(16*pi/16) == 0
    assert tan(pi) == 0
    assert tan(-pi) == 0
    assert tan(2*pi) == 0
    assert tan(-2*pi) == 0
    assert tan(-3*10**73*pi) == 0
    assert tan(7*10**103*pi) == 0
    assert tan(pi/2) == tan(pi/2)
    assert tan(-pi/2) == -tan(pi/2)
    assert tan(5*pi/2) == tan(5*pi/2)
    assert tan(7*pi/2) == tan(7*pi/2)

    assert tan(pi/3) == sqrt(3)
    assert tan(-2*pi/3) == sqrt(3)

    assert tan(pi/4) == S.One
    assert tan(-pi/4) == -S.One
    assert tan(17*pi/4) == S.One
    assert tan(-3*pi/4) == S.One

    assert tan(pi/6) == 1/sqrt(3)
    assert tan(-pi/6) == -1/sqrt(3)
    assert tan(7*pi/6) == 1/sqrt(3)
    assert tan(-5*pi/6) == 1/sqrt(3)

    assert tan(pi/105) == tan(pi/105)
    assert tan(-pi/105) == -tan(pi/105)

# Case 4) a != 0 and b != 0
def test_tan_case4():
    x, y = symbols('xy')
    r = Symbol('r', real=True)
    k = Symbol('k', integer=True)

    assert tan(pi - y) == -tan(y)
    assert tan(pi + y) == tan(y)
    assert tan(2*pi - y) == -tan(y)
    assert tan(pi/2 + y) == -cot(y)
    assert tan(pi/2 - y) == cot(y)
    assert tan(-y + pi) == -tan(y)
    assert tan(pi - y + pi) == -tan(y)
    assert tan(pi + y + pi) == tan(y)
    assert tan(2*pi - y + pi) == -tan(y)
    assert tan(pi/2 + y + pi) == -cot(y)
    assert tan(pi/2 - y + pi) == cot(y)
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


# Case 1) a == 0 and b == 0
def test_csc_case1():
    assert csc(0) == zoo

# Case 2) a != 0 and b == 0
def test_csc_case2():
    x, y, n = symbols('x y n')
    r = Symbol('r', real=True)
    k = Symbol('k', integer=True)
    #assert csc(r).is_real == True

    assert csc(exp(10)-1) == csc(-1+exp(10))
    assert sin(S.NaN) == S.NaN
    # See Issue 1540 for why these two fail:
    #assert sin(oo*I) == oo*I
    #assert sin(-oo*I) == -oo*I

    assert csc(x) == csc(x)
    assert csc(-x) == -csc(x)

    assert csc(-y) == -csc(y)

    assert csc(asin(x)) == 1 / x
    assert csc(acos(x)) == 1 / sqrt(1 - x**2)
    assert csc(atan(x)) == sqrt(1 + x**2) / x
    assert csc(acsc(x)) == x
    assert csc(asec(x)) == sqrt(1 - 1/x**2)
    assert csc(acot(x)) == sqrt(1 + 1 / x**2)
    assert csc(-asin(x)) == -1 / x
    assert csc(-acos(x)) == -1 / sqrt(1 - x**2)
    assert csc(-atan(x)) == -sqrt(1 + x**2) / x
    assert csc(-acsc(x)) == -x
    assert csc(-asec(x)) == -sqrt(1 - 1/x**2)
    assert csc(-acot(x)) == -sqrt(1 + 1 / x**2)
    # Need to implement csch and sech
    #assert csc(pi*I) == csch(pi)*I
    #assert csc(-pi*I) == -csch(pi)*I

    assert csc(2**1024 * E) == csc(2**1024 * E)
    assert csc(-2**1024 * E) == -csc(2**1024 * E)


    assert csc(2 + 3*I) == csc(2 + 3*I)
    x = Symbol('x', real=True)
    #assert csc(x*I) == csch(x)*I

# Case 3) a == 0 and b != 0
def test_csc_case3():
    x, y, n = symbols('x y n')
    r = Symbol('r', real=True)
    k = Symbol('k', integer=True)
    assert csc(k*pi) == zoo
    assert csc(17*k*pi) == zoo

    #assert csc(k*pi*I) == csch(k*pi)*I
    assert csc(pi/6) == 2
    assert csc(pi/4) == sqrt(2)
    assert csc(pi/3) == 2/sqrt(3)
    assert csc(pi/2) == 1
    assert csc(0 + 4*pi) == zoo
    assert csc(pi/6 + 4*pi) == 2
    assert csc(pi/4 + 4*pi) == sqrt(2)
    assert csc(pi/3 + 4*pi) == 2/sqrt(3)
    assert csc(pi/2 + 4*pi) == 1

    assert csc(0 + 2*pi) == zoo
    assert csc(pi/6 + 2*pi) == 2
    assert csc(pi/4 + 2*pi) == sqrt(2)
    assert csc(pi/3 + 2*pi) == 2/sqrt(3)
    assert csc(pi/2 + 2*pi) == 1
    assert csc(3*pi/2) == -1
    assert csc(5*pi/2) == 1
    assert csc(1) == csc(1)
    assert csc(-1) == -csc(1)
    assert csc(pi) == zoo
    assert csc(-pi) == zoo
    assert csc(2*pi) == zoo
    assert csc(-2*pi) == zoo
    assert csc(-3*10**73*pi) == zoo
    assert csc(7*10**103*pi) == zoo

    assert csc(pi/2) == 1
    assert csc(-pi/2) == -1
    assert csc(5*pi/2) == 1
    assert csc(7*pi/2) == -1

    assert csc(pi/3) == 2/sqrt(3)
    assert csc(-2*pi/3) == -2/sqrt(3)
    assert csc(pi/4) == 2/sqrt(2)
    assert csc(-pi/4) == -2/sqrt(2)
    assert csc(17*pi/4) == 2/sqrt(2)
    assert csc(-3*pi/4) == -2/sqrt(2)

    assert csc(pi/6) == 2
    assert csc(-pi/6) == -2
    assert csc(7*pi/6) == -2
    assert csc(-5*pi/6) == -2

    assert csc(1*pi/5) == 1/sqrt((5 - sqrt(5)) / 8)
    assert csc(2*pi/5) == 1/sqrt((5 + sqrt(5)) / 8)
    assert csc(3*pi/5) == csc(2*pi/5)
    assert csc(4*pi/5) == csc(1*pi/5)
    assert csc(6*pi/5) == -csc(1*pi/5)
    assert csc(8*pi/5) == -csc(2*pi/5)

    assert csc(pi/10) == 1/sin(pi/10)
    assert csc(3*pi/10) == (sqrt(5)+1)/4
    assert csc(7*pi/10) == (sqrt(5)+1)/4
    assert csc(9*pi/10) == (sqrt(5)-1)/4
    assert csc(13*pi/10) == (-1-sqrt(5))/4
    assert csc(17*pi/10) == (-1-sqrt(5))/4
    assert csc(19*pi/10) == (1-sqrt(5))/4
    assert csc(-pi/10) == (-sqrt(5)+1)/4
    assert csc(-3*pi/10) == -(sqrt(5)+1)/4
    assert csc(-7*pi/10) == -(sqrt(5)+1)/4
    assert csc(-9*pi/10) == (-sqrt(5)+1)/4
    assert csc(-13*pi/10) == (1+sqrt(5))/4
    assert csc(-17*pi/10) == (1+sqrt(5))/4
    assert csc(-19*pi/10) == (-1+sqrt(5))/4

    assert csc(-1273*pi/5) == -csc(2*pi/5)

    assert csc(pi/105) == csc(pi/105)
    assert csc(-pi/105) == -csc(pi/105)

    assert csc(104*pi/105) == csc(pi/105)
    assert csc(106*pi/105) == -csc(pi/105)

    assert csc(-104*pi/105) == -csc(pi/105)
    assert csc(-106*pi/105) == csc(pi/105)

    assert csc(1*pi/9) == csc(pi/9)
    assert csc(2*pi/9) == csc(2*pi/9)
    assert csc(3*pi/9) == sqrt(3)/2
    assert csc(4*pi/9) == cos(pi/18)
    assert csc(5*pi/9) == cos(pi/18)
    assert csc(6*pi/9) == sqrt(3)/2
    assert csc(7*pi/9) == csc(2*pi/9)
    assert csc(8*pi/9) == csc(pi/9)
    assert csc(9*pi/9) == 0
    assert csc(10*pi/9) == -csc(pi/9)
    assert csc(11*pi/9) == -csc(2*pi/9)
    assert csc(12*pi/9) == -sqrt(3)/2
    assert csc(13*pi/9) == -cos(pi/18)
    assert csc(14*pi/9) == -cos(pi/18)
    assert csc(15*pi/9) == -sqrt(3)/2
    assert csc(16*pi/9) == -csc(2*pi/9)
    assert csc(17*pi/9) == -csc(pi/9)
    assert csc(18*pi/9) == 0

    assert csc(-1*pi/9) == -csc(pi/9)
    assert csc(-2*pi/9) == -csc(2*pi/9)
    assert csc(-3*pi/9) == -sqrt(3)/2
    assert csc(-4*pi/9) == -cos(pi/18)
    assert csc(-5*pi/9) == -cos(pi/18)
    assert csc(-6*pi/9) == -sqrt(3)/2
    assert csc(-7*pi/9) == -csc(2*pi/9)
    assert csc(-8*pi/9) == -csc(pi/9)
    assert csc(-9*pi/9) == 0
    assert csc(-10*pi/9) == csc(pi/9)
    assert csc(-11*pi/9) == csc(2*pi/9)
    assert csc(-12*pi/9) == sqrt(3)/2
    assert csc(-13*pi/9) == cos(pi/18)
    assert csc(-14*pi/9) == cos(pi/18)
    assert csc(-15*pi/9) == sqrt(3)/2
    assert csc(-16*pi/9) == csc(2*pi/9)
    assert csc(-17*pi/9) == csc(pi/9)
    assert csc(-18*pi/9) == 0

# Case 4) a != 0 and b != 0
def test_csc_case4():
    x, y, n = symbols('x y n')
    r = Symbol('r', real=True)
    k = Symbol('k', integer=True)
    assert sin(-y + 2*pi) == -sin(y)
    assert sin(pi - y + 2*pi) == sin(y)
    assert sin(pi + y + 2*pi) == -sin(y)
    assert sin(2*pi - y + 2*pi) == -sin(y)
    assert sin(pi/2 + y + 2*pi) == cos(y)
    assert sin(pi/2 - y + 2*pi) == cos(y)

    assert sin(-y + 4*pi) == -sin(y)
    assert sin(pi - y + 4*pi) == sin(y)
    assert sin(pi + y + 4*pi) == -sin(y)
    assert sin(2*pi - y + 4*pi) == -sin(y)
    assert sin(pi/2 + y + 4*pi) == cos(y)
    assert sin(pi/2 - y + 4*pi) == cos(y)

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
    assert sin(pi - y) == sin(y)
    assert sin(pi + y) == -sin(y)
    assert sin(2*pi - y) == -sin(y)
    assert sin(pi/2 + y) == cos(y)
    assert sin(pi/2 - y) == cos(y)


# Case 1) a == 0 and b == 0
def test_cot_case1():
    assert cot(0) == zoo

# Case 2) a != 0 and b == 0
def test_cot_case2():
    assert cot(-y) == -cot(y)

# Case 4) a == 0 and b != 0
def test_cot_case3():
    assert cot(pi/6) == sqrt(3)
    assert cot(pi/4) == 1
    assert cot(pi/3) == 1/sqrt(3)
    assert cot(pi/2) == 0
    assert cot(pi/6 + pi) == sqrt(3)
    assert cot(pi/4 + pi) == 1
    assert cot(pi/3 + pi) == 1/sqrt(3)
    assert cot(pi/2 + pi) == 0
    assert cot(-16*pi/16) == zoo
    assert cot(-15*pi/16) == cot(pi/16)
    assert cot(-14*pi/16) == cot(pi/8)
    assert cot(-13*pi/16) == cot(3*pi/16)
    assert cot(-12*pi/16) == 1
    assert cot(-11*pi/16) == tan(3*pi/16)
    assert cot(-10*pi/16) == tan(pi/8)
    assert cot(-9*pi/16) == tan(pi/16)
    assert cot(-8*pi/16) == 0
    assert cot(-7*pi/16) == -tan(pi/16)
    assert cot(-6*pi/16) == -tan(pi/8)
    assert cot(-5*pi/16) == -tan(3*pi/16)
    assert cot(-4*pi/16) == -1
    assert cot(-3*pi/16) == -cot(3*pi/16)
    assert cot(-2*pi/16) == -cot(pi/8)
    assert cot(-1*pi/16) == -cot(pi/16)
    assert cot(0*pi/16) == zoo
    assert cot(1*pi/16) == cot(pi/16)
    assert cot(2*pi/16) == cot(pi/8)
    assert cot(3*pi/16) == cot(3*pi/16)
    assert cot(4*pi/16) == 1
    assert cot(5*pi/16) == tan(3*pi/16)
    assert cot(6*pi/16) == tan(pi/8)
    assert cot(7*pi/16) == tan(pi/16)
    assert cot(8*pi/16) == 0
    assert cot(9*pi/16) == -tan(pi/16)
    assert cot(10*pi/16) == -tan(pi/8)
    assert cot(11*pi/16) == -tan(3*pi/16)
    assert cot(12*pi/16) == -1
    assert cot(13*pi/16) == -cot(3*pi/16)
    assert cot(14*pi/16) == -cot(pi/8)
    assert cot(15*pi/16) == -cot(pi/16)
    assert cot(16*pi/16) == zoo

# Case 4) a != 0 and b != 0
def test_cot_case4():
    x, y = symbols('x y')
    assert cot(-y + pi) == -cot(y)
    assert cot(pi - y + pi) == -cot(y)
    assert cot(pi + y + pi) == cot(y)
    assert cot(2*pi - y + pi) == -cot(y)
    assert cot(pi/2 + y + pi) == -tan(y)
    assert cot(pi/2 - y + pi) == tan(y)
    assert cot(pi - y) == -cot(y)
    assert cot(pi + y) == cot(y)
    assert cot(2*pi - y) == -cot(y)
    assert cot(pi/2 + y) == -tan(y)
    assert cot(pi/2 - y) == tan(y)
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
    x, y = symbols('x y')
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



def test_attributes():
    x = Symbol('x')
    assert sin(x).args[:] == (x,)
    assert sin(x).args[0] != sin
    assert sin(x).args[0] == x

def test_sincos_rewrite():
    x = Symbol("x")
    y = Symbol("y")
    assert sin(pi/2-x) == cos(x)
    assert sin(pi-x) == sin(x)
    assert cos(pi/2-x) == sin(x)
    assert cos(pi-x) == -cos(x)

    assert sin(-x-y) == -sin(x+y)
    assert sin(-x-1) == -sin(x+1)
    assert cos(-x-y) == cos(x+y)
    assert cos(-x-1) == cos(x+1)
    assert sin(x-y) == sin(x-y)
    assert sin(y-x) == sin(y-x)
    assert cos(y-x) == cos(y-x)
    assert cos(x-y) == cos(x-y)
