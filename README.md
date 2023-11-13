# python-matrix-algebra
A simple matrix algebra library for Python, written in both Python and C!

```python
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
```

This project is an attempt at making matrix operations computationally cheaper by using simpler
representations on the C-side of FFI. See how matrices are represented and indexed in the C source:
```c
typedef struct Matrix
{
    size_t rows;
    size_t cols;
    double data[];
} 
Matrix;

static inline double *getIndex(size_t i, size_t j, Matrix *m) 
{
    size_t index = indexMatrix(i, j, m);
    return &(m->data[index]);
}
```

This allows for less indirection when accessing values and less time `malloc`ing memory when initializing matrices.
