#ifndef RAPTOR_CORE_BSR_MATRIX_HPP
#define RAPTOR_CORE_BSR_MATRIX_HPP

#include "src/core/matrix.hpp"

/**************************************************************
 *****   BSRMatrix Class (Inherits from CSRMatrix Base Class)
 **************************************************************
 ***** This class constructs a sparse matrix in BSR format.
 *****
 ***** Methods
 ***** -------
 ***** format() 
 *****    Returns the format of the sparse matrix (BSR)
 ***** sort()
 *****    Sorts the matrix.  Already in row-wise order, but sorts
 *****    the columns in each row.
 **************************************************************/
namespace raptor
{
  class CSRMatrix;
  class BSCMatrix;

  class BSRMatrix : public CSRMatrix
  {
    public:
        BSRMatrix(int _nrows, int _ncols, int _brows, 
                int _bcols, int _nnz = 0) 
        {
            n_rows = _nrows;
            n_cols = _ncols;
            b_rows = _brows;
            b_cols = _bcols;
            b_size = b_rows * b_cols;
            nnz = 0;
            sorted = false;
            diag_first = false;

            rowptr.resize(n_rows + 1);
            if (_nnz)
            {
                cols.reserve(_nnz);
                data.reserve(_nnz);
            }
        }

        BSRMatrix(int _nrows, int _ncols, int _brows, 
                int _bcols, double* _data) 
        {
            n_rows = _nrows;
            n_cols = _ncols;
            b_rows = _brows;
            b_cols = _bcols;
            b_size = b_rows * b_cols;
            sorted = false;
            diag_first = false;

            rowptr.resize(n_rows+1);

            nnz = 0;
            rowptr[0] = 0;
            for (int i = 0; i < n_rows; i++)
            {
                for (int j = 0; j < n_cols; j++)
                {
                    double sum = 0;
                    double* block = new double[b_size]();
                    for (int bi = 0; bi < b_rows; bi++)
                    {
                        int row = i*b_rows+bi;
                        for (int bj = 0; bj < b_cols; bj++)
                        {
                            int col = j*b_cols+bj;
                            double val = _data[row*n_cols*b_cols+col];
                            block[bi*b_cols+bj] = val;
                            sum += fabs(val);
                        }
                    }
                    if (sum > zero_tol)
                    {
                        cols.push_back(j);
                        data.push_back(block);
                        nnz++;
                    }
                    else
                    {
                        delete[] block;
                    }
                }
                rowptr[i+1] = nnz;
            }
        }

        BSRMatrix(CSRMatrix* A, int _brows, int _bcols)
        {
            b_rows = _brows;
            b_cols = _bcols;
            b_size = b_rows * b_cols;
            n_rows = A->n_rows / b_rows;
            n_cols = A->n_cols / b_cols;
            rowptr.resize(n_rows+1);

            // Map for each block row to accumulate block columns
            std::map<int, double*> block_map;

            rowptr[0] = 0;
            for (int i = 0; i < A->n_rows; i++)
            {
                int bi = i / b_rows; 
                int brow = i % b_rows;
                int row_start = A->rowptr[i];
                int row_end = A->rowptr[i+1];

                for (int j = row_start; j < row_end; j++)
                {
                    int col = A->cols[j];
                    int bj = col / b_cols;
                    int bcol = col % b_cols;

                    // map key `bj` is block column to be stored in BSR matrix
                    if (block_map.find(bj) == block_map.end())
                        block_map[bj] = new double[b_size]();

                    // Add value into appropriate block / positon
                    block_map[bj][brow*b_cols + bcol] = A->data[j];
                }

                // Last row in a block
                if ((i+1)%b_rows == 0)
                {
                    for (auto& [bj, block] : block_map)
                    {
                        // Push cols and data
                        cols.push_back(bj);
                        data.push_back(block);
                    }
                    rowptr[bi+1] = cols.size();
                    block_map.clear();
                }
            }
            nnz = cols.size();
        }


        ~BSRMatrix()
        {
            for (int i = 0; i < (int) data.size(); i++)
            {
                delete data[i];
            }
        }

        BSRMatrix* transpose()
        {
            BSCMatrix* T_bsc = new BSCMatrix(n_cols, n_rows, 
                    b_cols, b_rows, rowptr, cols, data);
            BSRMatrix* T = new BSRMatrix(T_bsc);
            delete T_bsc;
            return T;
        }

        format_t format()
        {
            return BSR;
        }


        void sort()
        {
            sort_helper(data);
        }

        void move_diag()
        {
            move_diag_helper(data);
        }

        void remove_duplicates()
        {
            remove_duplicates_helper(data);
        }

        void print()
        {
            int row_start, row_end; 
            for (int row = 0; row < n_rows; row++)
            {
                row_start = rowptr[row];
                row_end = rowptr[row+1];
                for (int j = row_start; j < row_end; j++)
                {
                    double* val = data[j];
                    for (int brow = 0; brow < b_rows; brow++)
                    {
                        for (int bcol = 0; bcol < b_cols; bcol++)
                        {
                            printf("A[%d][%d], Block[%d][%d] = %e\n", row, cols[j], brow, bcol, val[brow*b_cols+bcol]);
                        }
                    }
                }
            }
        }

        // Linear Algebra
        void spmv(const double alpha, const double* x,
                const double beta, double* b);
        void spmv_T(const double alpha, const double* x,
                const double beta, double* b);
        void add(BSRMatrix* B, BSRMatrix* C, double alpha, bool duplicates);
        void add(BSRMatrix* B, BSRMatrix* C, bool duplicates);
        void subtract(BSRMatrix* B, BSRMatrix* C);

        // Underlying storage will be contiguous
        std::vector<double*> data;
  };
}

#endif

