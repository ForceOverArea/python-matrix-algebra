from pymatlab import Matrix

# initialization is simple:
a = Matrix([   
    [1, 2, 3],   
    [4, 5, 6],   
])

b = Matrix([
    [ 7,  8],
    [ 9, 10],
    [11, 12],
])

c = a * b           # matrix multiplication is as easy as it should be

c *= 0.5            # scaling is easy too...
c.scale(2)          # ...and it's faster when you do it in-place

# ci = c.invert()     # inversion is also easy

d = 5 * c           # it is encouraged to multiply matrices like this...
e = c * 5           # ...but this is also valid

print(d == e)       # equality checking is easy

print(f"{c}")       # matrices have MATLAB-like string formatting...
print(f"{c:full}")  # ...and can be formatted as well