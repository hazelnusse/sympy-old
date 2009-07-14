from t import *
from sympy import S, pi, I, oo
x = symbols('x')
k = Symbol('k', integer=True, odd=True)
for n in range(-16, 17):
    print tan(n*pi/16)
stop

print sin(1*pi/9)
print sin(2*pi/9)
print sin(3*pi/9)
print sin(4*pi/9)
print sin(5*pi/9)
print sin(6*pi/9)
print sin(7*pi/9)
print sin(8*pi/9)
print sin(9*pi/9)
print sin(10*pi/9)
print sin(11*pi/9)
print sin(12*pi/9)
print sin(13*pi/9)
print sin(14*pi/9)
print sin(15*pi/9)
print sin(16*pi/9)
print sin(17*pi/9)
print sin(18*pi/9)

stop

print sin(0)
print sin(pi/12)
print sin(pi/6)
print sin(pi/4)
stop

print cos(0)
print cos(pi/12)
print cos(pi/6)
print cos(pi/4)

stop

print 'sin(x + 0) =', sin(x + 0)
print 'Octant 1'
print 'sin(x + pi/8) =', sin(x + pi/8)
print 'sin(x - 15*pi/8) =', sin(x - 15*pi/8)
print 'Octant 2'
print 'sin(x + 3*pi/8) =', sin(x + 3*pi/8)
print 'sin(x - 13*pi/8) =', sin(x - 13*pi/8)
print 'Octant 3'
print 'sin(x + 5*pi/8) =', sin(x + 5*pi/8)
print 'sin(x - 11*pi/8) =', sin(x - 11*pi/8)
print 'Octant 4'
print 'sin(x + 7*pi/8) =', sin(x + 7*pi/8)
print 'sin(x - 9*pi/8) =', sin(x - 9*pi/8)
print 'Octant 5'
print 'sin(x + 9*pi/8) =', sin(x + 9*pi/8)
print 'sin(x - 7*pi/8) =', sin(x - 7*pi/8)

print 'Octant 6'
print 'sin(x + 11*pi/8) =', sin(x + 11*pi/8)
print 'sin(x - 5*pi/8) =', sin(x - 5*pi/8)

print 'Octant 7'
print 'sin(x + 13*pi/8) =', sin(x + 13*pi/8)
print 'sin(x - 3*pi/8) =', sin(x - 3*pi/8)

print 'Octant 8'
print 'sin(x + 15*pi/8) =', sin(x + 15*pi/8)
print 'sin(x - pi/8) =', sin(x - pi/8)

print 'cos(x + 0) =', cos(x + 0)
print 'Octant 1'
print 'cos(x + pi/8) =', cos(x + pi/8)
print 'cos(x - 15*pi/8) =', cos(x - 15*pi/8)
print 'Octant 2'
print 'cos(x + 3*pi/8) =', cos(x + 3*pi/8)
print 'cos(x - 13*pi/8) =', cos(x - 13*pi/8)
print 'Octant 3'
print 'cos(x + 5*pi/8) =', cos(x + 5*pi/8)
print 'cos(x - 11*pi/8) =', cos(x - 11*pi/8)
print 'Octant 4'
print 'cos(x + 7*pi/8) =', cos(x + 7*pi/8)
print 'cos(x - 9*pi/8) =', cos(x - 9*pi/8)
print 'Octant 5'
print 'cos(x + 9*pi/8) =', cos(x + 9*pi/8)
print 'cos(x - 7*pi/8) =', cos(x - 7*pi/8)

print 'Octant 6'
print 'cos(x + 11*pi/8) =', cos(x + 11*pi/8)
print 'cos(x - 5*pi/8) =', cos(x - 5*pi/8)

print 'Octant 7'
print 'cos(x + 13*pi/8) =', cos(x + 13*pi/8)
print 'cos(x - 3*pi/8) =', cos(x - 3*pi/8)

print 'Octant 8'
print 'cos(x + 15*pi/8) =', cos(x + 15*pi/8)
print 'cos(x - pi/8) =', cos(x - pi/8)
stop

print cos(0)
print tan(0)
#print cot(0)

# Case 2: a + b*pi with a == 0 and b != 0
print '3*pi'
print sin(3*pi)
print cos(3*pi)
print tan(3*pi)
#print cot(3*pi)
print '3*pi/2'
print sin(3*pi/2)
print cos(3*pi/2)
print tan(3*pi/2)
#print cot(3*pi/2)
stop


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

#ex3 = tan(7*pi/18)+tan(5*pi/18)-sqrt(3)*tan(5*pi/18)*tan(7*pi/18)



"""
var("a b")
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
"""

