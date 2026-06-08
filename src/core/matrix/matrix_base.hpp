// Copyright (c) 2015-2026, RAPtor Developer Team
// License: Simplified BSD, http://opensource.org/licenses/BSD-2-Clause
#ifndef RAPTOR_CORE_MATRIX_BASE_HPP
#define RAPTOR_CORE_MATRIX_BASE_HPP

#include "src/core/types.hpp"
#include "src/core/vector.hpp"

/**************************************************************
 *****   Matrix Base Class
 **************************************************************
 ***** This class constructs a sparse matrix, supporting simple linear
 ***** algebra operations.
 *****
 ***** Attributes
 ***** -------------
 ***** n_rows : int
 *****    Number of rows
 ***** n_cols : int
 *****    Number of columns
 ***** nnz : int
 *****    Number of nonzeros
 ***** Methods
 ***** -------
 ***** resize(int n_rows, int n_cols)
 *****    Resizes dimension of matrix to passed parameters
 *****
 ***** Virtual Methods
 ***** -------
 ***** format() 
 *****    Returns the format of the sparse matrix (COO, CSR, CSC)
 ***** sort()
 *****    Sorts the matrix by position.  Whether row-wise or 
 *****    column-wise depends on matrix format.
 ***** add_value(int row, int col, double val)
 *****     Adds val to position (row, col)
 *****     TODO -- make sure this is working for CSR/CSC
 **************************************************************/
namespace raptor
{
  class Matrix
  {
    public:
        Matrix(int _nrows, int _ncols, int _nnz = 0)
        {
            b_rows = 1;
            b_cols = 1;
            b_size = 1;

            n_rows = _nrows;
            n_cols = _ncols;
            nnz = _nnz;
            sorted = false;
            diag_first = false;
        }

        Matrix()
        {
            b_rows = 1;
            b_cols = 1;
            b_size = 1;

            n_rows = 0;
            n_cols = 0;
            nnz = 0;
            sorted = false;
            diag_first = false;
        }

        virtual ~Matrix(){}

        virtual format_t format() = 0;
        virtual void sort() = 0;
        virtual void move_diag() = 0;
        virtual void remove_duplicates() = 0;
        virtual void print() = 0;

        virtual void spmv(const double alpha, const double* x, 
                const double beta, double* b) = 0;
        virtual void spmv_T(const double alpha, const double* x,
                const double beta, double* b) = 0;
        void mult(const double* x, double* b);
        void mult_T(const double* x, double* b);
        void residual(const double* x, const double* b, double* r);

        void mult(Vector& x, Vector& b)
        {
            mult(x.data(), b.data());
        }
        void mult_T(Vector& x, Vector& b)
        {
            spmv_T(1.0, x.data(), 0.0, b.data());
        }


        int n_rows;
        int n_cols;
        int nnz;

        int b_rows;
        int b_cols;
        int b_size;

        bool sorted;
        bool diag_first;
  };

}

#endif
