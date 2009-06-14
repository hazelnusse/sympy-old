from sympy import SYMPY_DEBUG

from sympy.core import Basic, S, C, Add, Mul, Pow, Rational, Integer, \
        Derivative, Wild, Symbol, sympify

from sympy.core.numbers import igcd

from sympy.utilities import make_list, all
from sympy.functions import gamma, exp, sqrt

from sympy.simplify.cse_main import cse

from sympy.polys import Poly

import sympy.mpmath as mpmath

def fraction(expr, exact=False):
    """Returns a pair with expression's numerator and denominator.
       If the given expression is not a fraction then this function
       will assume that the denominator is equal to one.

       This function will not make any attempt to simplify nested
       fractions or to do any term rewriting at all.

       If only one of the numerator/denominator pair is needed then
       use numer(expr) or denom(expr) functions respectively.

       >>> from sympy import *
       >>> x, y = symbols('x', 'y')

       >>> fraction(x/y)
       (x, y)
       >>> fraction(x)
       (x, 1)

       >>> fraction(1/y**2)
       (1, y**2)

       >>> fraction(x*y/2)
       (x*y, 2)
       >>> fraction(Rational(1, 2))
       (1, 2)

       This function will also work fine with assumptions:

       >>> k = Symbol('k', negative=True)
       >>> fraction(x * y**k)
       (x, y**(-k))

       If we know nothing about sign of some exponent and 'exact'
       flag is unset, then structure this exponent's structure will
       be analyzed and pretty fraction will be returned:

       >>> fraction(2*x**(-y))
       (2, x**y)

       #>>> fraction(exp(-x))
       #(1, exp(x))

       >>> fraction(exp(-x), exact=True)
       (exp(-x), 1)

    """
    expr = sympify(expr)

    numer, denom = [], []

    for term in make_list(expr, Mul):
        if term.is_Pow:
            if term.exp.is_negative:
                if term.exp is S.NegativeOne:
                    denom.append(term.base)
                else:
                    denom.append(Pow(term.base, -term.exp))
            elif not exact and term.exp.is_Mul:
                coeff, tail = term.exp.args[0], Mul(*term.exp.args[1:])#term.exp.getab()

                if coeff.is_Rational and coeff.is_negative:
                    denom.append(Pow(term.base, -term.exp))
                else:
                    numer.append(term)
            else:
                numer.append(term)
        elif term.func is C.exp:
            if term.args[0].is_negative:
                denom.append(C.exp(-term.args[0]))
            elif not exact and term.args[0].is_Mul:
                coeff, tail = term.args[0], Mul(*term.args[1:])#term.args.getab()

                if coeff.is_Rational and coeff.is_negative:
                    denom.append(C.exp(-term.args[0]))
                else:
                    numer.append(term)
            else:
                numer.append(term)
        elif term.is_Rational:
            if term.is_integer:
                numer.append(term)
            else:
                numer.append(Rational(term.p))
                denom.append(Rational(term.q))
        else:
            numer.append(term)

    return Mul(*numer), Mul(*denom)

def numer(expr):
    return fraction(expr)[0]

def denom(expr):
    return fraction(expr)[1]

def fraction_expand(expr):
    a, b = fraction(expr)
    return a.expand() / b.expand()

def numer_expand(expr):
    a, b = fraction(expr)
    return a.expand() / b

def denom_expand(expr):
    a, b = fraction(expr)
    return a / b.expand()

def separate(expr, deep=False):
    """Rewrite or separate a power of product to a product of powers
       but without any expanding, ie. rewriting products to summations.

       >>> from sympy import *
       >>> x, y, z = symbols('x', 'y', 'z')

       >>> separate((x*y)**2)
       x**2*y**2

       >>> separate((x*(y*z)**3)**2)
       x**2*y**6*z**6

       >>> separate((x*sin(x))**y + (x*cos(x))**y)
       x**y*cos(x)**y + x**y*sin(x)**y

       #this does not work because of exp combining
       #>>> separate((exp(x)*exp(y))**x)
       #exp(x*y)*exp(x**2)

       >>> separate((sin(x)*cos(x))**y)
       cos(x)**y*sin(x)**y

       Notice that summations are left un touched. If this is not the
       requested behaviour, apply 'expand' to input expression before:

       >>> separate(((x+y)*z)**2)
       z**2*(x + y)**2

       >>> separate((x*y)**(1+z))
       x**(1 + z)*y**(1 + z)

    """
    expr = sympify(expr)

    if expr.is_Pow:
        terms, expo = [], separate(expr.exp, deep)

        if expr.base.is_Mul:
            t = [ separate(C.Pow(t,expo), deep) for t in expr.base.args ]
            return C.Mul(*t)
        elif expr.base.func is C.exp:
            if deep == True:
                return C.exp(separate(expr.base[0], deep)*expo)
            else:
                return C.exp(expr.base[0]*expo)
        else:
            return C.Pow(separate(expr.base, deep), expo)
    elif expr.is_Add or expr.is_Mul:
        return type(expr)(*[ separate(t, deep) for t in expr.args ])
    elif expr.is_Function and deep:
        return expr.func(*[ separate(t) for t in expr.args])
    else:
        return expr


def together(expr, deep=False):
    """Combine together and denest rational functions into a single
       fraction. By default the resulting expression is simplified
       to reduce the total order of both numerator and denominator
       and minimize the number of terms.

       Densting is done recursively on fractions level. However this
       function will not attempt to rewrite composite objects, like
       functions, interior unless 'deep' flag is set.

       By definition, 'together' is a complement to 'apart', so
       apart(together(expr)) should left expression unhanged.

       >>> from sympy import *
       >>> x, y, z = symbols('x', 'y', 'z')

       You can work with sums of fractions easily. The algorithm
       used here will, in an iterative style, collect numerators
       and denominator of all expressions involved and perform
       needed simplifications:

       #>>> together(1/x + 1/y)
       #(x + y)/(y*x)

       #>>> together(1/x + 1/y + 1/z)
       #(z*x + x*y + z*y)/(y*x*z)

       >>> together(1/(x*y) + 1/y**2)
       (x + y)/(x*y**2)

       Or you can just denest multi-level fractional expressions:

       >>> together(1/(1 + 1/x))
       x/(1 + x)

       It also perfect possible to work with symbolic powers or
       exponential functions or combinations of both:

       >>> together(1/x**y + 1/x**(y-1))
       x**(-y)*(1 + x)

       #>>> together(1/x**(2*y) + 1/x**(y-z))
       #x**(-2*y)*(1 + x**(y + z))

       #>>> together(1/exp(x) + 1/(x*exp(x)))
       #(1+x)/(x*exp(x))

       #>>> together(1/exp(2*x) + 1/(x*exp(3*x)))
       #(1+exp(x)*x)/(x*exp(3*x))

    """

    def _together(expr):

        from sympy.core.function import Function

        if expr.is_Add:
            items, coeffs, basis = [], [], {}

            for elem in expr.args:
                numer, q = fraction(_together(elem))

                denom = {}

                for term in make_list(q.expand(), Mul):
                    expo = S.One
                    coeff = S.One

                    if term.is_Pow:
                        if term.exp.is_Rational:
                            term, expo = term.base, term.exp
                        elif term.exp.is_Mul:
                            coeff, tail = term.exp.as_coeff_terms()
                            if coeff.is_Rational:
                                tail = C.Mul(*tail)
                                term, expo = Pow(term.base, tail), coeff
                        coeff = S.One
                    elif term.func is C.exp:
                        if term.args[0].is_Rational:
                            term, expo = S.Exp1, term.args[0]
                        elif term.args[0].is_Mul:
                            coeff, tail = term.args[0].as_coeff_terms()
                            if coeff.is_Rational:
                                tail = C.Mul(*tail)
                                term, expo = C.exp(tail), coeff
                        coeff = S.One
                    elif term.is_Rational:
                        coeff = Integer(term.q)
                        term = Integer(term.p)

                    if term in denom:
                        denom[term] += expo
                    else:
                        denom[term] = expo

                    if term in basis:
                        total, maxi = basis[term]

                        n_total = total + expo
                        n_maxi = max(maxi, expo)

                        basis[term] = (n_total, n_maxi)
                    else:
                        basis[term] = (expo, expo)

                    coeffs.append(coeff)
                items.append((numer, denom))

            numerator, denominator = [], []

            for (term, (total, maxi)) in basis.iteritems():
                basis[term] = (total, total-maxi)

                if term.func is C.exp:
                    denominator.append(C.exp(maxi*term.args[0]))
                else:
                    if maxi is S.One:
                        denominator.append(term)
                    else:
                        denominator.append(Pow(term, maxi))

            if all([ c.is_integer for c in coeffs ]):
                gcds = lambda x, y: igcd(int(x), int(y))
                common = Rational(reduce(gcds, coeffs))
            else:
                common = S.One

            product = reduce(lambda x, y: x*y, coeffs) / common

            for ((numer, denom), coeff) in zip(items, coeffs):

                expr, coeff = [], product / (coeff*common)

                for term in basis.iterkeys():
                    total, sub = basis[term]

                    if term in denom:
                        expo = total-denom[term]-sub
                    else:
                        expo = total-sub

                    if term.func is C.exp:
                        expr.append(C.exp(expo*term.args[0]))
                    else:
                        if expo is S.One:
                            expr.append(term)
                        else:
                            expr.append(Pow(term, expo))

                numerator.append(coeff*Mul(*([numer] + expr)))

            return Add(*numerator)/(product*Mul(*denominator))
        elif expr.is_Mul or expr.is_Pow:
            return type(expr)(*[ _together(t) for t in expr.args ])
        elif expr.is_Function and deep:
            return expr.func(*[ _together(t) for t in expr.args ])
        else:
            return expr

    return _together(separate(expr))

#apart -> partial fractions decomposition (will be here :)

def collect(expr, syms, evaluate=True, exact=False):
    """Collect additive terms with respect to a list of symbols up
       to powers with rational exponents. By the term symbol here
       are meant arbitrary expressions, which can contain powers,
       products, sums etc. In other words symbol is a pattern
       which will be searched for in the expression's terms.

       This function will not apply any redundant expanding to the
       input expression, so user is assumed to enter expression in
       final form. This makes 'collect' more predictable as there
       is no magic behind the scenes. However it is important to
       note, that powers of products are converted to products of
       powers using 'separate' function.

       There are two possible types of output. First, if 'evaluate'
       flag is set, this function will return a single expression
       or else it will return a dictionary with separated symbols
       up to rational powers as keys and collected sub-expressions
       as values respectively.

       >>> from sympy import *
       >>> x, y, z = symbols('x', 'y', 'z')
       >>> a, b, c = symbols('a', 'b', 'c')

       This function can collect symbolic coefficients in polynomial
       or rational expressions. It will manage to find all integer or
       rational powers of collection variable:

       >>> collect(a*x**2 + b*x**2 + a*x - b*x + c, x)
       c + x*(a - b) + x**2*(a + b)

       The same result can achieved in dictionary form:

       >>> d = collect(a*x**2 + b*x**2 + a*x - b*x + c, x, evaluate=False)
       >>> d[x**2]
       a + b
       >>> d[x]
       a - b
       >>> d[sympify(1)]
       c

       You can also work with multi-variate polynomials. However
       remember that this function is greedy so it will care only
       about a single symbol at time, in specification order:

       >>> collect(x**2 + y*x**2 + x*y + y + a*y, [x, y])
       x*y + y*(1 + a) + x**2*(1 + y)

       Also more complicated expressions can be used as patterns:

       >>> collect(a*sin(2*x) + b*sin(2*x), sin(2*x))
       (a + b)*sin(2*x)

       >>> collect(a*x*log(x) + b*(x*log(x)), x*log(x))
       x*(a + b)*log(x)

       It is also possible to work with symbolic powers, although
       it has more complicated behaviour, because in this case
       power's base and symbolic part of the exponent are treated
       as a single symbol:

       >>> collect(a*x**c + b*x**c, x)
       a*x**c + b*x**c

       >>> collect(a*x**c + b*x**c, x**c)
       x**c*(a + b)

       However if you incorporate rationals to the exponents, then
       you will get well known behaviour:

       >>> collect(a*x**(2*c) + b*x**(2*c), x**c)
       (x**2)**c*(a + b)

       Note also that all previously stated facts about 'collect'
       function apply to the exponential function, so you can get:

       >>> collect(a*exp(2*x) + b*exp(2*x), exp(x))
       (a + b)*exp(2*x)

       If you are interested only in collecting specific powers
       of some symbols then set 'exact' flag in arguments:

       >>> collect(a*x**7 + b*x**7, x, exact=True)
       a*x**7 + b*x**7

       >>> collect(a*x**7 + b*x**7, x**7, exact=True)
       x**7*(a + b)

       You can also apply this function to differential equations, where
       derivatives of arbitary order can be collected:

       >>> from sympy import Derivative as D
       >>> f = Function('f') (x)

       >>> collect(a*D(f,x) + b*D(f,x), D(f,x))
       (a + b)*D(f(x), x)

       >>> collect(a*D(D(f,x),x) + b*D(D(f,x),x), D(f,x))
       (a + b)*D(f(x), x, x)

       >>> collect(a*D(D(f,x),x) + b*D(D(f,x),x), D(f,x), exact=True)
       a*D(f(x), x, x) + b*D(f(x), x, x)

       >>> collect(a*D(D(f,x),x) + b*D(D(f,x),x) + a*D(f,x) + b*D(f,x), D(f,x))
       (a + b)*D(f(x), x) + (a + b)*D(f(x), x, x)

       Or you can even match both derivative order and exponent at time:

       >>> collect(a*D(D(f,x),x)**2 + b*D(D(f,x),x)**2, D(f,x)) == \
           (a + b)*D(D(f,x),x)**2
       True


    Notes
    =====
        - arguments are expected to be in expanded form, so you might have to call
           .expand prior to calling this function.
    """
    def make_expression(terms):
        product = []

        for term, rat, sym, deriv in terms:
            if deriv is not None:
                var, order = deriv

                while order > 0:
                    term, order = Derivative(term, var), order-1

            if sym is None:
                if rat is S.One:
                    product.append(term)
                else:
                    product.append(Pow(term, rat))
            else:
                product.append(Pow(term, rat*sym))

        return Mul(*product)

    def parse_derivative(deriv):
        # scan derivatives tower in the input expression and return
        # underlying function and maximal differentiation order
        expr, sym, order = deriv.expr, deriv.symbols[0], 1

        for s in deriv.symbols[1:]:
            if s == sym:
                order += 1
            else:
                raise NotImplementedError('Improve MV Derivative support in collect')

        while isinstance(expr, Derivative):
            s0 = expr.symbols[0]

            for s in expr.symbols:
                if s != s0:
                    raise NotImplementedError('Improve MV Derivative support in collect')

            if s0 == sym:
                expr, order = expr.expr, order+len(expr.symbols)
            else:
                break

        return expr, (sym, Rational(order))

    def parse_term(expr):
        """Parses expression expr and outputs tuple (sexpr, rat_expo, sym_expo, deriv)
        where:
         - sexpr is the base expression
         - rat_expo is the rational exponent that sexpr is raised to
         - sym_expo is the symbolic exponent that sexpr is raised to
         - deriv contains the derivatives the the expression

         for example, the output of x would be (x, 1, None, None)
         the output of 2**x would be (2, 1, x, None)
        """
        rat_expo, sym_expo = S.One, None
        sexpr, deriv = expr, None

        if expr.is_Pow:
            if isinstance(expr.base, Derivative):
                sexpr, deriv = parse_derivative(expr.base)
            else:
                sexpr = expr.base

            if expr.exp.is_Rational:
                rat_expo = expr.exp
            elif expr.exp.is_Mul:
                coeff, tail = expr.exp.as_coeff_terms()

                if coeff.is_Rational:
                    rat_expo, sym_expo = coeff, C.Mul(*tail)
                else:
                    sym_expo = expr.exp
            else:
                sym_expo = expr.exp
        elif expr.func is C.exp:
            if expr.args[0].is_Rational:
                sexpr, rat_expo = S.Exp1, expr.args[0]
            elif expr.args[0].is_Mul:
                coeff, tail = expr.args[0].as_coeff_terms()

                if coeff.is_Rational:
                    sexpr, rat_expo = C.exp(C.Mul(*tail)), coeff
        elif isinstance(expr, Derivative):
            sexpr, deriv = parse_derivative(expr)

        return sexpr, rat_expo, sym_expo, deriv

    def parse_expression(terms, pattern):
        """Parse terms searching for a pattern.
        terms is a list of tuples as returned by parse_terms
        pattern is an expression
        """
        pattern = make_list(pattern, Mul)

        if len(terms) < len(pattern):
            # pattern is longer than  matched product
            # so no chance for positive parsing result
            return None
        else:
            pattern = [ parse_term(elem) for elem in pattern ]

            elems, common_expo, has_deriv = [], None, False

            for elem, e_rat, e_sym, e_ord in pattern:
                if e_ord is not None:
                    # there is derivative in the pattern so
                    # there will by small performance penalty
                    has_deriv = True

                for j in range(len(terms)):
                    term, t_rat, t_sym, t_ord = terms[j]

                    if elem.is_Number:
                        # a constant is a match for everything
                        break

                    if elem == term and e_sym == t_sym:
                        if exact == False:
                            # we don't have to exact so find common exponent
                            # for both expression's term and pattern's element
                            expo = t_rat / e_rat

                            if common_expo is None:
                                # first time
                                common_expo = expo
                            else:
                                # common exponent was negotiated before so
                                # teher is no chance for pattern match unless
                                # common and current exponents are equal
                                if common_expo != expo:
                                    common_expo = 1
                        else:
                            # we ought to be exact so all fields of
                            # interest must match in very details
                            if e_rat != t_rat or e_ord != t_ord:
                                continue

                        # found common term so remove it from the expression
                        # and try to match next element in the pattern
                        elems.append(terms[j])
                        del terms[j]

                        break

                else:
                    # pattern element not found
                    return None
            return terms, elems, common_expo, has_deriv

    if evaluate:
        if expr.is_Mul:
            ret = 1
            for term in expr.args:
                ret *= collect(term, syms, True, exact)
            return ret
        elif expr.is_Pow:
            b = collect(expr.base, syms, True, exact)
            return C.Pow(b, expr.exp)

    summa = [ separate(i) for i in make_list(sympify(expr), Add) ]

    if isinstance(syms, list):
        syms = [ separate(s) for s in syms ]
    else:
        syms = [ separate(syms) ]

    collected, disliked = {}, S.Zero

    for product in summa:
        terms = [ parse_term(i) for i in make_list(product, Mul) ]

        for symbol in syms:
            if SYMPY_DEBUG:
                print "DEBUG: parsing of expression %s with symbol %s " % (str(terms), str(symbol))

            result = parse_expression(terms, symbol)

            if SYMPY_DEBUG:
                print "DEBUG: returned %s" %  str(result)

            if result is not None:
                terms, elems, common_expo, has_deriv = result

                # when there was derivative in current pattern we
                # will need to rebuild its expression from scratch
                if not has_deriv:
                    index = 1
                    for elem in elems:
                        index *= Pow(elem[0], elem[1])
                        if elem[2] is not None:
                            index **= elem[2]
                else:
                    index = make_expression(elems)

                terms = separate(make_expression(terms))
                index = separate(index)

                if index in collected.keys():
                    collected[index] += terms
                else:
                    collected[index] = terms

                break
        else:
            # none of the patterns matched
            disliked += product

    if disliked is not S.Zero:
        collected[S.One] = disliked

    if evaluate:
        return Add(*[ a*b for a, b in collected.iteritems() ])
    else:
        return collected

def ratsimp(expr):
    """
    Usage
    =====
        ratsimp(expr) -> joins two rational expressions and returns the simples form

    Notes
    =====
        Currently can simplify only simple expressions, for this to be really usefull
        multivariate polynomial algorithms are needed

    Examples
    ========
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> y = Symbol('y')
        >>> e = ratsimp(1/x + 1/y)
    """
    expr = sympify(expr)
    if expr.is_Pow:
        return Pow(ratsimp(expr.base), ratsimp(expr.exp))
    elif expr.is_Mul:
        res = []
        for x in expr.args:
            res.append( ratsimp(x) )
        return Mul(*res)
    elif expr.is_Function:
        return expr.func(*[ ratsimp(t) for t in expr.args ])

    #elif expr.is_Function:
    #    return type(expr)( ratsimp(expr[0]) )
    elif not expr.is_Add:
        return expr

    def get_num_denum(x):
        """Matches x = a/b and returns a/b."""
        a,b = map(Wild, 'ab')
        r = x.match(a/b)
        if r is not None and len(r) == 2:
            return r[a],r[b]
        return x, S.One

    x = expr.args[0]
    y = Add(*expr.args[1:])

    a,b = get_num_denum(ratsimp(x))
    c,d = get_num_denum(ratsimp(y))

    num = a*d+b*c
    denum = b*d

    # Check to see if the numerator actually expands to 0
    if num.expand() == 0:
        return 0

    return num/denum

def trigsimp(expr, deep=False, recursive=False):
    """
    Usage
    =====
        trigsimp(expr) -> reduces expression by using known trig identities

    Notes
    =====

    deep ........ apply trigsimp inside functions
    recursive ... use common subexpression elimination (cse()) and apply
                  trigsimp recursively (recursively==True is quite expensive
                  operation if the expression is large)

    Examples
    ========
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> y = Symbol('y')
        >>> e = 2*sin(x)**2 + 2*cos(x)**2
        >>> trigsimp(e)
        2
        >>> trigsimp(log(e))
        log(2*cos(x)**2 + 2*sin(x)**2)
        >>> trigsimp(log(e), deep=True)
        log(2)
    """
    # Implement the rule based algorithm presented in:
    # Automated and readable simplification of trigonometric expressions
    # by Hongguang Fu, Xiuqin Zhong, Zhenbing Zeng
    # Mathematical and Computer Modelling 44 (2006), 1169-1177

    a, b = map(Wild, 'ab')

    # A dictionary of rules.  Keys are 'TR1' - 'TR13', values are each a tuple
    # of length 2 tuples, the first entry being the on the left hand side of a
    # trig identity, the second entry being on the right hand side of the same
    # trig identity.
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
                          (cos(a) - cos(b), -2*sin((a + b)/2)*sin((a - b)/2)))
                  'TR10': 


    """
    from sympy.core.basic import S
    sin, cos, tan, cot = C.sin, C.cos, C.tan, C.cot
    if recursive:
        w, g = cse(expr)
        g = trigsimp_nonrecursive(g[0])
        for sub in reversed(w):
            g = g.subs(sub[0], sub[1])
            g = trigsimp_nonrecursive(g)
        result = g
    else:
        result = trigsimp_nonrecursive(expr, deep)

    # do some final simplifications like sin/cos -> tan:
    a,b,c = map(Wild, 'abc')
    matchers = (
            (a*sin(b)**c/cos(b)**c, a*tan(b)**c),
    )
    for pattern, simp in matchers:
        res = result.match(pattern)
        if res is not None:
            # if c is missing or zero, do nothing:
            if (not c in res) or res[c] == 0:
                continue
            # if "a" contains the argument of sin/cos "b", skip the
            # simplification:
            if res[a].has(res[b]):
                continue
            # simplify and finish:
            result = simp.subs(res)
            break

    return result
    """

def trigsimp_nonrecursive(expr, deep=False):
    """
    A nonrecursive trig simplifier, used from trigsimp.

    Usage
    =====
        trigsimp_nonrecursive(expr) -> reduces expression by using known trig
                                       identities

    Notes
    =====

    deep ........ apply trigsimp inside functions

    Examples
    ========
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> y = Symbol('y')
        >>> e = 2*sin(x)**2 + 2*cos(x)**2
        >>> trigsimp(e)
        2
        >>> trigsimp_nonrecursive(log(e))
        log(2*cos(x)**2 + 2*sin(x)**2)
        >>> trigsimp_nonrecursive(log(e), deep=True)
        log(2)
    """
    from sympy.core.basic import S
    sin, cos, tan, cot = C.sin, C.cos, C.tan, C.cot

    if expr.is_Function:
        if deep:
            return expr.func( trigsimp_nonrecursive(expr.args[0], deep) )
    elif expr.is_Mul:
        ret = S.One
        for x in expr.args:
            ret *= trigsimp_nonrecursive(x, deep)

        return ret
    elif expr.is_Pow:
        return Pow(trigsimp_nonrecursive(expr.base, deep),
                trigsimp_nonrecursive(expr.exp, deep))
    elif expr.is_Add:
        # TODO this needs to be faster

        # The types of trig functions we are looking for
        a,b,c = map(Wild, 'abc')
        matchers = (
            (a*sin(b)**2, a - a*cos(b)**2),
            (a*tan(b)**2, a*(1/cos(b))**2 - a),
            (a*cot(b)**2, a*(1/sin(b))**2 - a)
        )

        # Scan for the terms we need
        ret = S.Zero
        for term in expr.args:
            term = trigsimp_nonrecursive(term, deep)
            res = None
            for pattern, result in matchers:
                res = term.match(pattern)
                if res is not None:
                    ret += result.subs(res)
                    break
            if res is None:
                ret += term

        # Reduce any lingering artifacts, such as sin(x)**2 changing
        # to 1-cos(x)**2 when sin(x)**2 was "simpler"
        artifacts = (
            (a - a*cos(b)**2 + c, a*sin(b)**2 + c, cos),
            (a - a*(1/cos(b))**2 + c, -a*tan(b)**2 + c, cos),
            (a - a*(1/sin(b))**2 + c, -a*cot(b)**2 + c, sin)
        )

        expr = ret
        for pattern, result, ex in artifacts:
            # Substitute a new wild that excludes some function(s)
            # to help influence a better match. This is because
            # sometimes, for example, 'a' would match sec(x)**2
            a_t = Wild('a', exclude=[ex])
            pattern = pattern.subs(a, a_t)
            result = result.subs(a, a_t)
            if expr.is_number:
                continue

            m = expr.match(pattern)
            while m is not None:
                if m[a_t] == 0 or -m[a_t] in m[c].args or m[a_t] + m[c] == 0:
                    break
                expr = result.subs(m)
                m = expr.match(pattern)

        return expr
    return expr

def radsimp(expr):
    """
    Rationalize the denominator.

    Examples:
    =========
        >>> from sympy import *
        >>> radsimp(1/(2+sqrt(2)))
        1 - 2**(1/2)/2
        >>> x,y = map(Symbol, 'xy')
        >>> e = ( (2+2*sqrt(2))*x+(2+sqrt(8))*y )/( 2+sqrt(2) )
        >>> radsimp(e)
        x*2**(1/2) + y*2**(1/2)
    """
    n,d = fraction(expr)
    a,b,c = map(Wild, 'abc')
    r = d.match(a+b*sqrt(c))
    if r is not None:
        a = r[a]
        if r[b] == 0:
            b,c = 0,0
        else:
            b,c = r[b],r[c]

        syms = list(n.atoms(Symbol))
        n = collect( (n*(a-b*sqrt(c))).expand(), syms )
        d = a**2 - c*b**2

    return n/d

def powsimp(expr, deep=False):
    """
    Usage
    =====
        powsimp(expr, deep) -> reduces expression by combining powers with
            similar bases and exponents.

    Notes
    =====
        If deep is True then powsimp() will also simplify arguments of
        functions. By default deep is set to False.


    Examples
    ========
        >>> from sympy import *
        >>> x,n = map(Symbol, 'xn')
        >>> e = x**n * (x*n)**(-n) * n
        >>> powsimp(e)
        n**(1 - n)

        >>> powsimp(log(e))
        log(n*x**n*(n*x)**(-n))

        >>> powsimp(log(e), deep=True)
        log(n**(1 - n))
    """
    def _powsimp(expr):
        if expr.is_Pow:
            if deep:
                return C.Pow(powsimp(expr.base), powsimp(expr.exp))
            return expr
        elif expr.is_Function and deep:
            return expr.func(*[powsimp(t) for t in expr.args])
        elif expr.is_Add:
            return C.Add(*[powsimp(t) for t in expr.args])
        elif expr.is_Mul:
            # Collect base/exp data, while maintaining order in the
            # non-commutative parts of the product
            c_powers = {}
            nc_part = []
            for term in expr.args:
                if term.is_commutative:
                    b,e = term.as_base_exp()
                    c_powers[b] = c_powers.get(b, 0) + e
                else:
                    nc_part.append(term)

            # Pull out numerical coefficients from exponent
            for b,e in c_powers.items():
                exp_c, exp_t = e.as_coeff_terms()
                if not (exp_c is S.One) and exp_t:
                    del c_powers[b]
                    new_base = C.Pow(b, exp_c)
                    if new_base in c_powers:
                        c_powers[new_base] += C.Mul(*exp_t)
                    else:
                        c_powers[new_base] = C.Mul(*exp_t)

            # Combine bases whenever they have the same exponent which is
            # not numeric
            c_exp = {}
            for b, e in c_powers.items():
                if e in c_exp:
                    c_exp[e].append(b)
                else:
                    c_exp[e] = [b]

            # Merge back in the results of the above to form a new product
            for e in c_exp:
                bases = c_exp[e]
                if len(bases) > 1:
                    for b in bases:
                        del c_powers[b]
                    new_base = Mul(*bases)
                    if new_base in c_powers:
                        c_powers[new_base] += e
                    else:
                        c_powers[new_base] = e

            c_part = [ C.Pow(b,e) for b,e in c_powers.items() ]
            return C.Mul(*(c_part + nc_part))
        return expr

    return _powsimp(separate(expr, deep=deep))

def hypersimp(f, k):
    """Given combinatorial term f(k) simplify its consecutive term ratio
       ie. f(k+1)/f(k).  The input term can be composed of functions and
       integer sequences which have equivalent representation in terms
       of gamma special function.

       The algorithm performs three basic steps:

           (1) Rewrite all functions in terms of gamma, if possible.

           (2) Rewrite all occurrences of gamma in terms of products
               of gamma and rising factorial with integer,  absolute
               constant exponent.

           (3) Perform simplification of nested fractions, powers
               and if the resulting expression is a quotient of
               polynomials, reduce their total degree.

       If f(k) is hypergeometric then as result we arrive with a
       quotient of polynomials of minimal degree. Otherwise None
       is returned.

       For more information on the implemented algorithm refer to:

       [1] W. Koepf, Algorithms for m-fold Hypergeometric Summation,
           Journal of Symbolic Computation (1995) 20, 399-417
    """
    f = sympify(f)

    g = f.subs(k, k+1) / f

    g = g.rewrite(gamma)
    g = g.expand(func=True, basic=False)

    if g.is_rational_function(k):
        return Poly.cancel(g, k)
    else:
        return None

def hypersimilar(f, g, k):
    """Returns True if 'f' and 'g' are hyper-similar.

       Similarity in hypergeometric sens means that a quotient of
       f(k) and g(k) is a rational function in k.  This procedure
       is useful in solving recurrence relations.

       For more information see hypersimp().

    """
    f, g = map(sympify, (f, g))

    h = (f/g).rewrite(gamma)
    h = h.expand(func=True, basic=False)

    return h.is_rational_function(k)

def combsimp(expr):
    return expr

def simplify(expr):
    """Naively simplifies the given expression.

       Simplification is not a well defined term and the exact strategies
       this function tries can change in the future versions of SymPy. If
       your algorithm relies on "simplification" (whatever it is), try to
       determine what you need exactly  -  is it powsimp(), or radsimp(),
       or together(), or something else? And use this particular function
       directly, because those are well defined and thus your algorithm
       will be robust.

    """
    expr = Poly.cancel(powsimp(expr))
    return together(expr.expand())

def nsimplify(expr, constants=[], tolerance=None, full=False):
    """
    Numerical simplification -- tries to find a simple formula
    that numerically matches the given expression. The input should
    be possible to evalf to a precision of at least 30 digits.

    Optionally, a list of (rationally independent) constants to
    include in the formula may be given.

    A lower tolerance may be set to find less exact matches.

    With full=True, a more extensive search is performed
    (this is useful to find simpler numbers when the tolerance
    is set low).

    Examples:

        >>> from sympy import *
        >>> nsimplify(4/(1+sqrt(5)), [GoldenRatio])
        -2 + 2*GoldenRatio
        >>> nsimplify((1/(exp(3*pi*I/5)+1)))
        1/2 - I*(1/4 + 5**(1/2)/10)**(1/2)
        >>> nsimplify(I**I, [pi])
        exp(-pi/2)
        >>> nsimplify(pi, tolerance=0.01)
        22/7

    """
    expr = sympify(expr)

    prec = 30
    bprec = int(prec*3.33)

    constants_dict = {}
    for constant in constants:
        constant = sympify(constant)
        v = constant.evalf(prec)
        if not v.is_Real:
            raise ValueError("constants must be real-valued")
        constants_dict[str(constant)] = v._to_mpmath(bprec)

    exprval = expr.evalf(prec, chop=True)
    re, im = exprval.as_real_imag()

    # Must be numerical
    if not ((re.is_Real or re.is_Integer) and (im.is_Real or im.is_Integer)):
        return expr

    def nsimplify_real(x):
        orig = mpmath.mp.dps
        xv = x._to_mpmath(bprec)
        try:
            # We'll be happy with low precision if a simple fraction
            if not (tolerance or full):
                mpmath.mp.dps = 15
                rat = mpmath.findpoly(xv, 1)
                if rat is not None:
                    return Rational(-int(rat[1]), int(rat[0]))
            mpmath.mp.dps = prec
            newexpr = mpmath.identify(xv, constants=constants_dict,
                tol=tolerance, full=full)
            if not newexpr:
                raise ValueError
            if full:
                newexpr = newexpr[0]
            return sympify(newexpr)
        finally:
            mpmath.mp.dps = orig
    try:
        if re: re = nsimplify_real(re)
        if im: im = nsimplify_real(im)
    except ValueError:
        return expr

    return re + im*S.ImaginaryUnit


