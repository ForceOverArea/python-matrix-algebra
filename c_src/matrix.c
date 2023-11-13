#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

/// @brief A contiguous chunk of memory representing a matrix. 
typedef struct Matrix
{
    size_t rows;
    size_t cols;
    double data[];
} 
Matrix;

/// @brief Calculates the index to access in the 1-dimensional array representation.
/// @param i the row to access, starting at index 0
/// @param j the column to access, starting at index 1
/// @param m the matrix pointer being accessed.
/// @return the 1-D index that should be accessed
static inline size_t indexMatrix(size_t i, size_t j, Matrix *m)
{
    return i * m->cols + j;
}

/// @brief Returns a pointer to a specific index of a `Matrix`
/// @param i the row to access, starting at 0
/// @param j the column to access, starting at 0
/// @param m the `Matrix` from which to get the value
/// @return a `double *` to the given `Matrix` index
static inline double *getIndex(size_t i, size_t j, Matrix *m) 
{
    size_t index = indexMatrix(i, j, m);
    return &(m->data[index]);
}

/// @brief Initializes a new pointer to a contiguous chunk of memory in the form of a matrix.
///        If this function fails to `malloc` a new pointer, it will return `NULL`.
/// @param rows the number of rows to include in the matrix
/// @param cols the number of columns to include in the matrix
/// @return `NULL` or a new `Matrix` pointer.
Matrix *newMatrix(size_t rows, size_t cols)
{
    size_t indices = rows * cols;
    Matrix *result = (Matrix *)malloc(sizeof(Matrix) + sizeof(double) * indices);
    if (!result)
        return NULL;

    result->rows = rows;
    result->cols = cols;
    for (size_t i = 0; i < indices; i++)
        result->data[i] = 0;

    return result;
}

/// @brief Serves the same purpose as `newMatrix`, but returns a pointer to an initialized 
///        identity matrix of size n
/// @param n the size of the identity matrix
/// @return a pointer to a new identity `Matrix`
Matrix *newIdentity(size_t n)
{
    size_t indices = sizeof(double) * (n * n);
    Matrix *result = (Matrix *)malloc(sizeof(Matrix) + indices);
    if (!result)
        return NULL;

    result->rows = n;
    result->cols = n;

    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            if (i == j)
                *getIndex(i, j, result) = 1;
            else
                *getIndex(i, j, result) = 0;

    return result;
}

/// @brief Adds a `double[]` to the given row of the given `Matrix *`,
///        returning a `bool` indicating whether the operation could be 
///        performed. If `vecToAdd` is not the same width as the matrix OR
///        if `row` is out-of-bounds, this operation will fail.
/// @param vecToAdd the `double[]` to add element-wise to the given `Matrix` row.
/// @param row the index of the row to add to (starts at 0)
/// @param m the `Matrix *` to add to
/// @return an indicator of success
bool addToRow(double vecToAdd[], size_t row, Matrix *m)
{
    if (row > m->cols)
        return false;

    for (size_t j = 0; j < m->cols; j++)
        *getIndex(row, j, m) += vecToAdd[j];

    return true;
}

/// @brief Scales a given row by a given scalar value, returning a `bool` indicating 
///        whether the operation could be performed. If `row` is out-of-bounds, this 
///        operation will fail.
/// @param scalar the scalar value to scale the row by.
/// @param row the row (starting at 0) to scale
/// @param m the `Matrix` to operate on. 
/// @return an indicator of success
bool scaleRow(double scalar, size_t row, Matrix *m)
{
    if (row > m->cols)
        return false;
    
    for (size_t j = 0; j < m->cols; j++)
        *getIndex(row, j, m) *= scalar;

    return true;
}

/// @brief 
/// @param a 
/// @param b 
/// @return 
bool matricesAreEqual(Matrix *a, Matrix *b)
{
    if (a->rows != b->rows || a->cols != b->cols)
        return false;

    for (size_t i = 0; i < a->rows; i++)
        for (size_t j = 0; j < a->rows; j++)
            if (*getIndex(i, j, a) != *getIndex(i, j, b))
                return false;

    return true;
}

/// @brief 
/// @param a 
/// @param b 
/// @return 
Matrix *augmentMatrix(Matrix *a, Matrix *b)
{
    if (a->rows != b->rows)
        return NULL;

    size_t rows = a->rows;
    size_t cols = a->cols + b->cols;
    Matrix *result = newMatrix(rows, cols);
    if (!result)
        return NULL;

    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            if (j < a->cols)
                result->data[indexMatrix(i, j, result)] 
                    = *getIndex(i, j,           a); // assign value from 'a'
            else
                result->data[indexMatrix(i, j, result)] 
                    = *getIndex(i, j - a->cols, b); // assign value from 'b'

    return result;
}

/// @brief Attempts to reduce a `Matrix` in place, mutating it's current values.
///        It is not recommended to call this function directly. Instead, consider 
///        whether `solveWith` or `invert` can be used. 
/// @param m the `Matrix` to reduce
/// @return a `bool` indicating whether the matrix could be successfully reduced.
bool reduceMatrix(Matrix *m)
{
    // only operate on the smaller number: # of rows OR # of cols
    size_t limit = (m->rows <= m->cols)? m->rows : m->cols;

    // upper right triangular
    for (size_t j = 0; j < limit; j++) // this is essentially just a counter here 
    {
        double a11 = *getIndex(j, j, m);
        for (size_t i = j + 1; i < limit; i++)
        {
            double denom = *getIndex(j, j, m);
            if (!denom)
                return false;
            double scalar  = a11 / denom;          // first item in matrix subsection
            double *rowPtr = &(m->data[indexMatrix(0, 0, m)]);   // address of contiguous 1st row vector
            scaleRow(-scalar, i, m);                            
            addToRow( rowPtr, i, m);
        }
    }

    // lower left triangular
    size_t n = limit - 1;
    for (size_t j = n; j >= 0; j--)
    {
        double ann = *getIndex(j, j, m);
        for (size_t i = j - 1; i >= 0; i--)
        {
            double denom = *getIndex(j, j, m);
            if (!denom)
                return false;
            double scalar  = ann / denom;          // first item in matrix subsection
            double *rowPtr = &(m->data[indexMatrix(0, n, m)]);  // address of contiguous last row vector
            scaleRow(-scalar, i, m);
            addToRow( rowPtr, i, m);
        }
    }

    // scale diagonal
    for (size_t j = 0; j < limit; j++)
    {
        double scalar = 1 / *getIndex(j, j, m); 
        scaleRow(scalar, j, m);
    }

    return true;
}

/// @brief Creates a new `Matrix` by copying values from a subset of another given
///        `Matrix`.
/// @param istart the first row to copy values from (inclusive) 
/// @param jstart the first column to copy values from (inclusive)
/// @param iend the last row to copy values from (inclusive)
/// @param jend the last column to copy values from (inclusive)
/// @param m the `Matrix` pointer to copy values from
/// @return a new `Matrix` containing the selected values
Matrix *subset(size_t istart, size_t jstart, size_t iend, size_t jend, Matrix *m)
{
    if (iend < istart || jend < jstart)
        return NULL;

    Matrix *result = newMatrix(iend + 1 - istart, jend + 1 - jstart);
    if (!result)
        return NULL;

    size_t currentIndex = 0;
    for (size_t i = istart; i <= iend; i++)
        for (size_t j = jstart; j <= jend; j++)
        {
            // if we did our math right, this is valid for all i, j.
            result->data[currentIndex] = *getIndex(i, j, m);             
            // printf("[%zu, %zu] (%zu) = %lf\n", i, j, currentIndex, result->data[currentIndex]);
            currentIndex++;
        }

    return result;
}

/// @brief Attempts to invert a matrix, preserving the original and returning a new
///        instance with the inverse value. This function returns `NULL` if the 
///        inversion process fails.
/// @param m the `Matrix` to invert
/// @return a new `Matrix` equal to the inverse of `m`
Matrix *invert(Matrix *m)
{
    if (m->rows != m->cols)
        return NULL;
    size_t n = m->rows;

    Matrix *ity = newIdentity(n);
    if (!ity)
        return NULL;

    // TODO: inline this part manually
    Matrix *aug = augmentMatrix(m, ity);
    free(ity); // we need to free this regardless of augment creation success
    if (!aug)
        return NULL;
    reduceMatrix(aug);

    size_t startRow = 0; 
    size_t startCol = n;
    size_t endRow   = n - 1;
    size_t endCol   = 2 * n - 1;

    Matrix *result = subset(startRow, startCol, endRow, endCol, aug);
    free(aug); // we need to free this regardless of de-augmentation success
    if (!result)
        return NULL;
    
    return result;
}

/// @brief Returns the product of two matrices
/// @param a the 'left' `Matrix`
/// @param b the 'right' `Matrix`
/// @return the product of `a` and `b` as a `Matrix`
Matrix *multiplyMatrices(Matrix *a, Matrix *b)
{
    if (a->cols != b->rows)
        return NULL;
    
    Matrix *result = newMatrix(a->rows, b->cols);
    if (!result)
        return NULL;

    for (size_t row = 0; row < a->rows; row++)  // for each row in a...
        for (size_t col = 0; col < b->cols; col++)  // for each col in b...
            for (size_t i   = 0; i   < a->cols;   i++)  // for each element in a's row and b's col...
                *getIndex(row, col, result) += (*getIndex(row,i,a)) * (*getIndex(i,col,b));

    return result;
}

/// @brief 
/// @param m 
/// @param solnVector 
/// @return 
Matrix *solveWith(Matrix *m, Matrix *solnVector)
{
    // TODO
    return NULL;
}