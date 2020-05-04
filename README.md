# Algebra With Python
This repository contains the code accompanying the book Algebra With Python. Specifcally it contains the script
`algebra.py`. 

Check out the install guide at https://youtu.be/FmU1C9lS7VQ

If using an interactive IPython, Jupyter or QTconsole session, it should begin as follows:

```python
from algebra import *

#  inequality sign inverses
inv_sign = { Ge : Le, 
            Gt : Lt,
            Le : Ge,
            Lt : Gt }

def add(num):
    """
    Typically adds a number to the global Eq object varibale `eq`. 
    Can also add a number to both sides of an inequality.
    """
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
    """
    Typically subtracts a number from the global Eq object varibale `eq`. 
    Can also subtracts a number from both sides of an inequality.Rules
    for inequality manipulation hardcoded.
    """
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
    """
    Typically multiples a each side of the global Eq object varibale `eq`
    by a number. Can also do the same for both sides of an inequality. Rules
    for inequality manipulation hardcoded.
    """
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
    """
    Typically divides a each side of the global Eq object varibale `eq`
    by a number. Can also do the same for both sides of an inequality.
    Rules for inequality manipulation hardcoded.
    """ 
    global eq
    comps = []
    for arg in eq.args:
        comps.append(arg / num)
    if isinstance(eq, Equality):
        eq =  Eq(comps[0],comps[1])
        
        # Check for non-canceling terms
        
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
    """
    Swaps the sides of the equation
    """
    global eq
    if isinstance(eq, Equality):
        eq = Eq(eq.args[1], eq.args[0])
        return eq                   
                    
                                 

```

Throughout the book more functions will be introduced which will also need to be included. This is emphasised in each relevant chapter.

