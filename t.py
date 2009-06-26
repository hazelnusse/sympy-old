from sympy import *

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

stop







from sympy import sin, cos, pi, Symbol, S, Function, solve, latex
x = Symbol('x')
print 'sin(pi/2 - x) =',sin(pi/2 - x)
stop



from pydy import *
class GeneralizedCoordinate(Symbol):
    def __init__(self, name, depends_on=Symbol('t'), *args):
        self.dv = depends_on
    """
    def diff(self, var=Symbol('t')):
        if var == self:
            print 'var == self'
            return 1
        elif var == self.dv:
            print 'var != self'
            return GeneralizedCoordinate(self.name+'\'')
        else:
            print 'none'

        if n == 1:
            print 'n = ', n
            return GeneralizedCoordinate(self.name+'\'')
        else:
            print 'n = ', n
            d = GeneralizedCoordinate(self.name+'\'')
            return d.dt(n=n-1)
    """
    def _eval_derivative(self, s):
        if s == self:
            return S.One
        elif s == self.dv:
            return GeneralizedCoordinate(self.name+'\'')
        else:
            return S.Zero

    def _latex_(self):
        return '\dot{'+str(self.name)+'}'

t = Symbol('t')
x = GeneralizedCoordinate('x')
print 'x = ',x
print 'diff(x, t) = ', x.diff(t)
print 'dt(sin(x))', sin(x).diff(t)
e = sin(x).diff(t)
print 'type(sin(x).diff(t))', type(e)
print latex(x)
print print_pydy(sin(x))
stop


print 'diff(sin(x), t)', diff(sin(x), t)
stop



t = Symbol('t')
#x = SymbolF('x', iv=t)
x = Symbol('x')
print 'cos(0*pi/2 + x + 2)', cos(0*pi/S(2) + x + 2)
print 'cos(1*pi/2 + x + 2)', cos(1*pi/S(2) + x + 2)
print 'cos(2*pi/2 + x + 2)', cos(2*pi/S(2) + x + 2)
print 'cos(3*pi/2 + x + 2)', cos(3*pi/S(2) + x + 2)

print 'cos(4*pi/2 + x + 2)', cos(4*pi/S(2) + x + 2)
print 'cos(5*pi/2 + x + 2)', cos(5*pi/S(2) + x + 2)
print 'cos(6*pi/2 + x + 2)', cos(6*pi/S(2) + x + 2)
print 'cos(7*pi/2 + x + 2)', cos(7*pi/S(2) + x + 2)

stop
t = Symbol('t')
x = Symbol('x')
y = Function('y')(t)
a11,a12,a21,a22,b1,b2 = symbols('a11','a12','a21','a22','b1','b2')
soln = solve([a11*x + a12*y - b1, a21*x + a22*y - b2], x, y)
print 'soln', soln
print 'should equal', { y : (a11*b2 - a21*b1)/(a11*a22 - a12*a21), x :
    (a22*b1 - a12*b2)/(a11*a22 - a12*a21) }

stop

x,y = symbols('xy')
print sin(x+2*pi)
print tan(pi/S(2))
stop


print resmin
print 'x=1,', resmin.subs(x,1)
print 'x=2,', resmin.subs(x,2)
print sqrt(x*y)


e = sqrt(1-x**2).subs(x, 2)
print e
assert e == sqrt(-3)
e = (sqrt(1-x**2)**3).subs(x, 2)
print e
assert e == -3*sqrt(3)*I
e = (sqrt(1-x**2)).subs(x, 2)**3
print e
assert e == -3*sqrt(3)*I

print (-3)**(S(1)/2)
print (-3)**(S(3)/2)
print (-3)**(S(5)/2)
print (-3)**(S(7)/2)
