from ctypes import CDLL, c_double, c_size_t, c_void_p, c_bool
from os import path
from typing import Self

# DLL (aka DLL) file for C-written matrix library
PYMATLAB_DIR = path.dirname(__file__)
PYMATLAB_DLL = path.join(PYMATLAB_DIR, "pymatlablib.so")
DLL          = CDLL(PYMATLAB_DLL)

# Specify required function signature details
DLL.getSubset.argtypes     = [c_size_t, c_size_t, c_size_t, c_size_t, c_void_p]
DLL.getSubset.restype      = c_void_p

DLL.getValAtIndex.argtypes = [c_size_t, c_size_t, c_void_p]
DLL.getValAtIndex.restype  = c_double

DLL.invertMatrix.argtypes  = [c_void_p]
DLL.invertMatrix.restype   = c_void_p

DLL.killMatrix.argtypes    = [c_void_p]

DLL.matrixProduct.argtypes = [c_void_p, c_void_p]
DLL.matrixProduct.restype  = c_void_p

DLL.matrixEquality.argtypes= [c_void_p, c_void_p]
DLL.matrixEquality.restype = c_bool

DLL.newIdentity.argtypes   = [c_size_t]
DLL.newIdentity.restype    = c_void_p

DLL.newMatrix.restype      = c_void_p

DLL.scaleRow.argtypes      = [c_double, c_size_t, c_void_p]
DLL.scaleMatrix.argtypes   = [c_double, c_void_p]

DLL.setValAtIndex.argtypes = [c_size_t, c_size_t, c_double, c_void_p]

def is_valid_matrix_repr(data: list[list[float]]) -> bool:
    n = len(data[0])
    for row in data[1:]:
        if len(row) != n or False in [type(x) in [float, int] for x in row]:
            return False
    return True

def matrix_indices_are_valid(indices: tuple) -> bool:
    return type(indices) == tuple \
       and len(indices) == 2 \
       and type(indices[0]) == int and type(indices[1]) == int \
       and indices[0] >= 0 and indices[1] >= 0

class Matrix:

    def __init__(self, 
        data: list[list[float]] | None = None, 
        /, 
        rows: int | None = None, 
        cols: int | None = None,
        ptr:  c_void_p | None = None
    ) -> Self:
        """
        Initializes a new `Matrix` object.

        `data` - if present and correctly formatted, the matrix will be initialized with the 
                 given values in the 2-D list.
        `rows`, `cols` - if present in the absence of `data`, these keyword arguments will 
                 determine the shape of the matrix.
        `ptr`  - if present, allows the caller to give a specific raw pointer to an
                 already-initialized `Matrix`. If not, a new Matrix containing only 0
                 will be initialized.
        """
        if data and type(data) == list and is_valid_matrix_repr(data):
            self.rows = len(data)
            self.cols = len(data[0])
            self.raw_ptr: c_void_p = DLL.newMatrix(self.rows, self.cols)

            for i in range(self.rows):
                for j in range(self.cols):
                    self.__setitem__((i, j), c_double(data[i][j]))

        elif (not data) and rows and cols:
            self.rows = rows
            self.cols = cols
            if type(ptr) == c_void_p and ptr != 0:
                self.raw_ptr = ptr
            else:
                self.raw_ptr: c_void_p = DLL.newMatrix(rows, cols)

        else:
            raise Exception 

    def __rmul__(self, val: int | float) -> Self:
        res = c_void_p(DLL.getSubset(0, 0, self.rows-1, self.cols-1, self.raw_ptr))
        if res == 0:
            return None
        
        result = Matrix(
            rows = self.rows, 
            cols = self.cols, 
            ptr  = res
        )
        result.scale(val)

        return result

    def __mul__(self, val: Self | int | float) -> Self | None:
        if type(val) == Matrix:
            res = c_void_p(DLL.matrixProduct(self.raw_ptr, val.raw_ptr))
            if res == 0:
                return None

            return Matrix(
                rows = self.rows,
                cols = val.cols,
                ptr  = res
            )

        elif type(val) in [int, float]:
            return self.__rmul__(val)

        else:
            raise Exception

    def __getitem__(self, indices: tuple) -> float:
        if not matrix_indices_are_valid(indices):
            raise Exception
            
        return c_double(DLL.getValAtIndex(indices[0], indices[1], self.raw_ptr)).value

    def __setitem__(self, indices: tuple, new_value: float) -> None:
        if not matrix_indices_are_valid(indices):
            raise Exception

        DLL.setValAtIndex(indices[0], indices[1], new_value, self.raw_ptr)

    def __str__(self) -> str:
        res = "["
        lastcol = self.cols - 1

        for i in range(self.rows):
            for j in range(lastcol):
                res += f"{self[i, j]}, "
            res += f"{self[i, lastcol]}; "

        res = res[:-2] + "]"

        return res

    def __format__(self, __format_spec: str) -> str:
        if __format_spec != "full":
            return self.__str__()
        
        res = "\n"
        lastcol = self.cols - 1

        for i in range(self.rows):
            res += "|"
            for j in range(lastcol):
                res += f"{self[i, j]:10.3f}, "
            res += f"{self[i, lastcol]:10.3f}|\n"

        return res

    def __eq__(self, val: Self) -> bool:
        if type(val) != Matrix:
            return False
        
        return bool(DLL.matrixEquality(self.raw_ptr, val.raw_ptr))

    def __del__(self) -> None:
        DLL.killMatrix(self.raw_ptr)

    def scale(self, scalar: int | float) -> None:
        if type(scalar) not in [int, float]:
            raise Exception

        DLL.scaleMatrix(scalar, self.raw_ptr) # mutates the pointer in `self`

    def invert(self) -> Self | None:
        res = c_void_p(DLL.invertMatrix(self.raw_ptr))
        if res == 0:
            return None

        return Matrix(
            rows = self.rows, 
            cols = self.cols, 
            ptr  = res
        )        

class Identity(Matrix):

    def __init__(self, n: int) -> Self:
        """
        """
        self.rows = self.cols = n
        self.raw_ptr = c_void_p(DLL.newIdentity(n));