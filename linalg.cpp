#include "linalg.h"
#include "blas_wrapper.h"
#include "lapacke_wrapper.h"

Matrix Matrix::operator+(const scalar_t val) {
    Matrix p(_row, _col); std::fill_n(p.data(), _row * _col, val);
    Matrix ret = *this; 
    daxpy(1.0, p, ret);
    return ret;
}

Matrix Matrix::operator+(Matrix& p) {
    assert(_row == p.row());
    assert(_col == p.col());
    Matrix ret = *this; 
    daxpy(1.0, p, ret);
    return ret;
}

Matrix Matrix::operator+(Matrix&& p) {
    assert(_row == p.row());
    assert(_col == p.col());
    Matrix ret = *this;
    daxpy(1.0, p, ret);
    return ret;
}

Matrix Matrix::operator-(Matrix& p) {
    assert(_row == p.row());
    assert(_col == p.col());
    Matrix ret = *this;
    daxpy(-1.0, p, ret);
    return ret;
}

Matrix Matrix::operator-(scalar_t val) {
    Matrix p(_row, _col); std::fill_n(p.data(), _row * _col, val);
    Matrix ret = *this;
    daxpy(-1.0, p, ret);
}

Matrix Matrix::operator-(Matrix&& p) {
    Matrix ret = *this;
    daxpy(-1.0, p, ret);
    return ret;
}

Matrix Matrix::operator*(scalar_t val) {
    Matrix ret = *this;
    dscal(val, ret);
    return ret;
}

Matrix Matrix::operator*(Matrix& val) {
    Matrix ret = *this;
    dgemm(1.0, *this, val, 0., ret) ;
    return ret;
}


Matrix Matrix::operator*(Matrix&& val) {
    Matrix ret = *this;
    dgemm(1.0, *this, val, 0., ret) ;
    return ret;
}