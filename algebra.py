"""
Module accompanying the book `Algebra with Python` by Simon Thomas.

Copyright 2020.


algebra.py - the main module containing useful functions.

"""
import numpy as np
import matplotlib.pyplot as plt

from sympy import *
init_session()


# Plotting Settings
from distutils.spawn import find_executable
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# Use latex if installed
if find_executable("latex"):
    rc('text', usetex=True)


# ------------------------------------------------- #
# 	CHAPTER 1 - INTRODUCTION TO ALGEBRA
# ------------------------------------------------- #

"""

The functions for Chapter 1 should be run in the working Python shell
and not imported via this script. This is due to the way global variables
work are not shared across modules in a manner suitable for this book.

You can just copy and paste the following codes:

COPY BELOW HERE -----------------

#  inequality sign inverses
inv_sign = { Ge : Le, 
            Gt : Lt,
            Le : Ge,
            Lt : Gt }

def add(num):
    "
    Typically adds a number to the global Eq object varibale `eq`. 
    Can also add a number to both sides of an inequality.
    "
    global eq
    comps = []
    for arg in eq.args:
        comps.append(arg + num)
    if isinstance(eq, Equality):
        eq = Eq(comps[0],comps[1])
        return eq
    else:
        for f in [Ge, Gt, Le, Lt]:
            if isinstance(eq, f):
                eq = f(comps[0], comps[1])
                return eq


def sub(num):
    "
    Typically subtracts a number from the global Eq object varibale `eq`. 
    Can also subtracts a number from both sides of an inequality.Rules
    for inequality manipulation hardcoded.
    "
    global eq
    comps = []
    for arg in eq.args:
        comps.append(arg - num)
    if isinstance(eq, Equality):
        eq =  Eq(comps[0],comps[1])
        return eq
    else:
        for f in [Ge, Gt, Le, Lt]:
            if isinstance(eq, f):
                eq = f(comps[0], comps[1])
                return eq

def mul(num):
    "
    Typically multiples a each side of the global Eq object varibale `eq`
    by a number. Can also do the same for both sides of an inequality. Rules
    for inequality manipulation hardcoded.
    "
    global eq
    comps = []
    for arg in eq.args:
        comps.append(arg * num)
    if isinstance(eq, Equality):
        eq =  Eq(comps[0],comps[1])
        return eq
    else:
        for f in [Ge, Gt, Le, Lt]:
            if isinstance(eq, f):
                if num < 0:
                    eq = inv_sign[f](comps[0], comps[1])
                    return eq
                else:
                    eq = f(comps[0], comps[1])
                    return eq


def div(num):
    "
    Typically divides a each side of the global Eq object varibale `eq`
    by a number. Can also do the same for both sides of an inequality.
    Rules for inequality manipulation hardcoded.
    " 
    global eq
    comps = []
    for arg in eq.args:
        comps.append(arg / num)
    if isinstance(eq, Equality):
        eq =  Eq(comps[0],comps[1])
        return eq
    else:
        for f in [Ge, Gt, Le, Lt]:
            if isinstance(eq, f):
                if num < 0:
                    eq = inv_sign[f](comps[0], comps[1])
                    return eq
                else:
                    eq = f(comps[0], comps[1])
                    return eq

def swap():
    "
    Switches the sides of the equation
    "
    global eq
    if isinstance(eq, Equality):
        eq = Eq(eq.args[1], eq.args[0])
        return eq

END COPY -----------------------------------------------------

"""

# ------------------------------------ #
# 	CHAPTER 2 - LINEAR EQUATIONS
# ------------------------------------ #
def find_x_inter(eq):
    """
    Finds the x-intercept of an equation by substituting
    in 0 for y and solving for x.
    
    Input:
        eq - equation object containing variables x and y
        
    Output:
        x-intercept - float/int
    """
    return (solve(eq.subs(Symbol("y"),0),Symbol("x"))[0],0)

def find_y_inter(eq):
    """
    Finds the y-intercept of an equation by substituting
    in 0 for x and solving for y.
    
    Input:
        eq - equation object containing variables x and y
        
    Output:
        y-intercept - float/int
    
    """
    return (0,solve(eq.subs(Symbol("x"),0),Symbol("y"))[0])

def slope_from_points(point1, point2):
    """
    Finds the slope of a line passing
    through two points following the formual
    
    slope = (y2-y1) / (x2-x1)
    
    Input:
        point1 - tuple of length two corresponding to (x1, y1)
        point2 - tuple of length two corresponding to (x2, y2)
        
    Ouput:
        slope - the slope of the line
    """
    x1, y1 = point1[0], point1[1]
    x2, y2 = point2[0], point2[1]
    return (y2-y1) / (x2-x1)

# ------------------------------------ #
# 	CHAPTER 3 - FUNCTIONS
# ------------------------------------ #


def Func(sym, var, exp):
    """
    Create a custom function
    
    Input:
        sym -> Function name/symbol
        var -> function input variable
        exp -> expression to evalue when function is valled
        
    Out:
        CustomFunction object
    """
    class CustomFunction(Function):
        @classmethod
        def eval(cls, x):
            if x.is_Number or x.is_Float:
                return exp.subs(var, x)
        
        def __str__(self):
            pass
    CustomFunction.__name__ = sym
    CustomFunction.exp = exp
    
    return CustomFunction
    

# ------------------------------------ #
# 	CHAPTER 6 - RATIONAL EXPONENTS
# ------------------------------------ #


def n_root(S, n, error = 1e-5):
    """
    Computes the n-th root of S using Newton's method
    of approximation.
    
    Input:
        S - number to find root
        n - n-th root to find
        error - desired precision of root. default : 1e-5
    Outout:
        root - n-th root of S
    """
    # Step 1: Create known answer
    s, x, e = symbols("s x e")
    formula = Eq(s, (x+e)**n)
    # Expand out brackets
    formula = expand(formula)
    # Subtract last term
    formula = Eq(formula.args[0] - x**n,
                 formula.args[1] - x**n)
    formula = factor(formula, e)
    # Divide by right hand side to find e (-1 to cancel e)
    formula = formula.args[0] / (formula.args[1]/e - 1)
    # Cancel out e term
    formula = formula.subs(e, 1)
    
    # Step 2: Approximate root
    x_n = 1.5 # Initial estimate
    # Iterate
    count = 0
    while True:
        x_n_plus_1 = x_n + formula.subs([(s, S),(x, x_n)])
        if x_n_plus_1 < 0:
            x_n_plus_1 *= -1
        if abs(x_n_plus_1 - x_n) <= error:
            return x_n
        x_n = x_n_plus_1.n(10)
        count += 1
        # Raise exception
        if count > 1000:
            raise Exception("Loop limit reached")
    
# ------------------------------------ #
# 	CHAPTER 7 - POLYNOMIALS
# ------------------------------------ #

def quadsolve(A, B, C):
    """
    Solves a quadratic by using the quadratic equation
    
    x = -b +- \sqrt{b^2 - 4*a*c} / 2*a
    
    derived from  ax^2 + bx + c = 0
    
    Input:
    
        A - a coefficient
        B - b coefficient
        C - c coefficient
        
    Output:
    
        solutions to x
    """
    x, a, b, c = symbols("x a b c")
    quadForm = solve(a*x**2 + b*x + c, x)
    one = quadForm[0].subs({"a":A, "b":B, "c":C})
    two = quadForm[1].subs({"a":A, "b":B, "c":C})
    return [one, two]


def findVertex(A, B):
    """
    Use the shortcut to find the x-value of the vertex.
    """
    return Eq(x, -B /(2*A))

# ------------------------------------ #
# 	CHAPTER 9 - COMPLEX FUNCTIONS
# ------------------------------------ #

def getInverse(eq):
    """
    Returns the inverse function of the
    equation defined as y = f(x)
    """
    inv = solve(eq, x)[0]
    inv =inv.subs(y, x)
    return Eq(y, inv)

'''
These functions also operator on the `eq` object and
so also need to be run in the main script rather than 
imported.

def cancel_exp(num):
    """
    Cancels the n-th exponentiated term by taking the n-th root
    on the opposite side.
    """
    global eq
    comps = []
    for arg in eq.args:
        if isinstance(arg, Pow):
            comps.append(arg.args[0])
        else:
            comps.append(root(arg, num))
    eq =  Eq(comps[0],comps[1])
    return eq

def cancel_root(num):
    """
    Cancels the n-th root term by  exponentiating the opposited by n.
    """
    global eq
    comps = []
    for arg in eq.args:
        if isinstance(arg, Pow):
            comps.append(arg.args[0])
        else:
            comps.append(arg**num)
    eq =  Eq(comps[0],comps[1])
    return eq

'''
    
# ------------------------------------ #
# 	CHAPTER 10 - COMPLEX NUMBERS
# ------------------------------------ #

def cplot(exps):
    """
    Plots a list of complex numbers in the complex plane
    """
    xs = []
    ys = []
    for exp in exps:
        if len(exp.args) > 1:
            xs.append(exp.args[0])
            ys.append(exp.args[1].args[0])
        else:
            # No imaginary component
            xs.append(exp)
            ys.append(0)
            
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(xs, ys, marker='o', color=["r", "g", "b", "orange"])
    
    # Move axis to center
    ax.spines['left'].set_position('zero')
    ax.spines['bottom'].set_position('zero')
    
    # Eliminate upper and right axes
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    # Show ticks in the left and lower axes only
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    
    plt.show()

# ------------------------------------ #
# 	CHAPTER 11 - MORE POLYNOMIALS
# ------------------------------------ #

from sympy.polys.polytools import div as divide

def choose(n, k):
    """
    n choose k  = n! / (k!(n-k)!)
    """
    return factorial(n) / (factorial(k) * factorial(n - k))

def generate_polynomial_expression(n):
    """
    returns (a + b)^n using the binomial theorum
    
    (a + b)^n = \sum_{k=0}^n {n \choose k} a^{n-k} \cdot b^k
    """
    a, b = symbols("a b")
    exp = 0
    for k in range(n+1):
        exp += choose(n, k) * a**(n-k) * b**k
    return exp

def pascals_triangle(n):
    """
    Generates Pascal's triangle by finding the coefficients
    of the number of unique combinations of variables in
    polynomials.
    
    Input:
        n -  maximum of nth degree polynomial to find coefficients for
    
    Output:
        coeffs - list of coefficients in the first nth degree polynomials
    """
    coeffs = []
    for count in range(n+1):
        row = [] # need a row element to collect the row in
        for element in range(count + 1): 
            row.append(choose(count, element))
        coeffs.append(row)
    return coeffs


# ----------------------------------------------------------- #
# 	CHAPTER 12 - MORE RADICAL AND RATIONAL EXPRESSIONS
# ----------------------------------------------------------- #

def reciprocal(exp):
    """
    Returns the reciprocal of the given expression /equations.
    """
    if isinstance(exp, Eq):
        left_n, left_d = exp.args[0].as_numer_denom()
        right_n, right_d = exp.args[1].as_numer_denom()
        return Eq(left_d/ left_n, right_d / right_n)
    num, denom = exp.as_numer_denom()
    return denom / num

# ----------------------------------------------------------- #
# 	CHAPTER 14 - TRIGONOMETRIC FUNCTIONS
# ----------------------------------------------------------- #

def degrees2radians(degrees):
    """
    Converts degrees to radians
    """
    return (pi / 180)*degrees

def radians2degrees(radians):
    """
    Converts radians to degrees
    """
    return (180 / pi)*radians

def UnitCircle(objects=[], functions=[]):
    """
    Creates a base unit circle which can added to.
    """

    fig, ax = plt.subplots(figsize=(8, 8))
    # Create circle
    unit_circle = plt.Circle((0,0),radius=1, fill=False, lw=1)   
    ax.add_artist(unit_circle)
    x_axis = plt.Line2D((-1.5, 1.5),(0,0),lw=1, color="black")
    y_axis = plt.Line2D((0,0),(-1.5,1.5), lw=1, color="black")
    ax.add_artist(x_axis)
    ax.add_artist(y_axis)
    plt.grid(True)
    
    # Add user defined objects
    for item in objects:
        ax.add_artist(item)
    
    # Set limits
    ax.set_ylim([-1.1, 1.1])
    ax.set_xlim([-1.1, 1.1])
    
    for func in functions:
        eval(func)

from matplotlib.patches import ConnectionPatch, Arc, Polygon
from matplotlib.text import Text

def Line(start, end, **kargs):
    """
    Creates a line from start to end. Can take all the usual keyword
    arguments e.g. color, linewidth etc.

    Used in combination with the UnitCircle function.
    """
    return ConnectionPatch(start, end, "data", **kargs)

def plot_sinusoid(func=np.sin, a=None, b=None, c=None, d=None, limit=(0, 4*np.pi), step=0.2, **kwargs):
    """
    Plots a sinusoid function and shows the axis as multiplies of pi.
    
    Inputs:
        func - the function to plot. Default is np.sin
        a - the amplitude coefficient
        b - the period coefficient
        c - the starting position term
        d - the midline term
        limit - the range of the function
        step - the resolution of the plot. The small the step, the longer the time to compute
        **kwards - any other matplotlib keyword argument you want to use. e.g. color = "red"
    Ouput:
        None
    """
    # Create x and y mapping
    if not a:
        a = 1
    if not b:
        b = 1
    if not c:
        c = 0
    if not d:
        d = 0
    start = limit[0]
    end = limit[1]
    x = np.arange(start, end+step, step)
    y = a * func(b*x + c) + d

    #
    fig, ax = plt.subplots()
    ax.plot(x, y, **kwargs)

    # Set spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_position('center')
    ax.spines['left'].set_position('zero')

    # Set xticks as multiples of pi/4 (45 degrees)
    ticks = [i*np.pi/2 for i in range(10)]
    # Create labels as latex expressions 
    labels = [ "$" + latex(Rational(1, 4)*pi * i) + "$" for i in range(10) ]
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, size=16)

    plt.show()


# ----------------------------------------------------------- #
# 	CHAPTER 16 - INTRODUCTION TO VECTORS
# ----------------------------------------------------------- #

def center_align(ax):
    # Move left y-axis and bottim x-axis to centre, passing through (0,0)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')

    # Eliminate upper and right axes
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    
    ax.grid()
    return ax


Vector = Matrix # create a new pointer to Matrix object called `Vector`

def magnitude(vector):
    """
    Returns the magnitude of the input vector.
    """
    return root(np.sum(np.square(vector)), 2)

Vector.mag = magnitude # Assign function as class method

def unit(vector):
    return vector / vector.mag()

Vector.unit = unit # Assign function as class method

i = Vector([1, 0])
j = Vector([0, 1])

def VectorMagDir(magnitude, direction):
    """
    Returns the components of the vector with a 
    given magnitude and direction. Direction is in
    degrees (NOT radians).
    """
    y = sin(degrees2radians(direction)) * magnitude
    x = cos(degrees2radians(direction)) * magnitude
    u = Vector([x, y])
    return u







