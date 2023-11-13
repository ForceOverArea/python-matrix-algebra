#ifndef MATRIX_H
#define MATRIX_H

#include <stdlib.h>
#include <stdbool.h>

typedef struct Matrix
{
    size_t rows;
    size_t cols;
    double data[];
} 
Matrix;

static inline size_t indexMatrix(size_t i, size_t j, Matrix *m)
{
    return i * m->cols + j;
}

static inline double *getIndex(size_t i, size_t j, Matrix *m) 
{
    size_t index = indexMatrix(i, j, m);
    return &(m->data[index]);
}

bool scaleRow(double scalar, size_t row, Matrix *m);

bool matricesAreEqual(Matrix *a, Matrix *b);

Matrix *subset(size_t istart, size_t jstart, size_t iend, size_t jend, Matrix *m);

Matrix *invert(Matrix *m);

Matrix *multiplyMatrices(Matrix *a, Matrix *b);

#endif