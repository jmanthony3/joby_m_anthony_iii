from types import FunctionType
# import types
import numpy
import sympy

add_one = lambda x: x + 1
if isinstance(add_one,FunctionType):
    print(type(add_one))
print(add_one(2))

exp = 2+sympy.Symbol('x')
print(exp)
if isinstance(exp,(FunctionType, sympy.Expr)):
    print(type(exp))

def test(a, w=0): 
    if w == 0: w = 'yes'
    return a, w
print(test(1))

print(sum([1, 2]))