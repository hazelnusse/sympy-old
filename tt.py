from t import *

# Case 1:  a + b*pi with a == b == 0
print sin(0)
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
