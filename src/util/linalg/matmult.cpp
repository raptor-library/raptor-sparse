#include "src/core/matrix.hpp"

using namespace raptor;


CSRMatrix* CSRMatrix::mult(const CSRMatrix* B, int* C_map)
{
    std::vector<int> next(B->n_cols, -1);
    std::vector<double> sums(B->n_cols, 0);

    CSRMatrix* C = new CSRMatrix(n_rows, B->n_cols);

    C->rowptr[0] = 0;
    for (int i = 0; i < n_rows; i++)
    {
        int head = -2;
        int length = 0;
        int row_start = rowptr[i];
        int row_end = rowptr[i+1];
        for (int j = row_start; j < row_end; j++)
        {
            int col = cols[j];
            double val = data[j];
            int row_start_B = B->rowptr[col];
            int row_end_B = B->rowptr[col+1];
            for (int k = row_start_B; k < row_end_B; k++)
            {
                int col_B = B->cols[k];
                sums[col_B] += val * B->data[k];
                if (next[col_B] == -1)
                {
                    next[col_B] = head;
                    head = col_B;
                    length++;
                }
            }
        }
        for (int j = 0; j < length; j++)
        {
            int col = head;
            head = next[col];
            next[col] = -1;
            double val = sums[col];
            sums[col] = 0;

            if (fabs(val) > zero_tol)
            {
                if (C_map)
                {
                    col = C_map[col];
                }
                C->cols.push_back(col); 
                C->data.push_back(val);       
            }
        }
        C->rowptr[i+1] = C->cols.size();
    }
    C->nnz = C->cols.size();

    return C;
}


CSRMatrix* CSRMatrix::mult_T(const CSRMatrix* B, int* C_map)
{
    CSCMatrix* AT = new CSCMatrix(this);

    std::vector<int> next(B->n_cols, -1);
    std::vector<double> sums(B->n_cols, 0);

    CSRMatrix* C = new CSRMatrix(n_cols, B->n_cols);

    C->rowptr[0] = 0;
    for (int i = 0; i < AT->n_cols; i++)
    {
        int head = -2;
        int length = 0;
        int col_start = AT->colptr[i];
        int col_end = AT->colptr[i+1];
        for (int j = col_start; j < col_end; j++)
        {
            int row = AT->rows[j];
            double val = AT->data[j];
            int row_start_B = B->rowptr[row];
            int row_end_B = B->rowptr[row+1];
            for (int k = row_start_B; k < row_end_B; k++)
            {
                int col = B->cols[k];
                sums[col] += val * B->data[k];
                if (next[col] == -1)
                {
                    next[col] = head;
                    head = col;
                    length++;
                }
            }
        }
        for (int j = 0; j < length; j++)
        {
            int col = head;
            double val = sums[col];
            head = next[col];
            next[col] = -1;
            sums[col] = 0;
            if (fabs(val) > zero_tol)
            {
                if (C_map)
                {
                    col = C_map[col];
                }
                C->cols.push_back(col);
                C->data.push_back(val);
            }
        }
        C->rowptr[i+1] = C->cols.size();
    }
    C->nnz = C->cols.size();

    delete AT;
    return C;
}




BSRMatrix* BSRMatrix::mult(const BSRMatrix* B, int* C_map)
{
    std::vector<int> next(B->n_cols, -1);
    std::vector<std::array<double, b_size>> sums(B->n_cols);
    for (std::array<double, b_size>& s : sums)
    {
        s.fill(0.0);
    }

    BSRMatrix<b_rows, b_cols>* C = new BSRMatrix<b_rows, b_cols>(n_rows, B->n_cols);

    C->rowptr[0] = 0;
    for (int i = 0; i < n_rows; i++)
    {
        int head = -2;
        int length = 0;
        int row_start = rowptr[i];
        int row_end = rowptr[i+1];
        for (int j = row_start; j < row_end; j++)
        {
            int col = cols[j];
            int row_start_B = B->rowptr[col];
            int row_end_B = B->rowptr[col+1];
            for (int k = row_start_B; k < row_end_B; k++)
            {
                int col_B = B->cols[k];
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, b_rows, b_cols, 
                        b_cols, 1.0, data[j].data(), b_cols, B->data[k].data(), b_cols,
                        1.0, sums[col_B].data(), b_cols);
                if (next[col_B] == -1)
                {
                    next[col_B] = head;
                    head = col_B;
                    length++;
                }
            }
        }
        for (int j = 0; j < length; j++)
        {
            int col = head;
            head = next[col];
            next[col] = -1;
            std::array<double, b_size> block = sums[col];

            if (abs_val(sums[col]) > zero_tol)
            {
                if (C_map)
                {
                    C->cols.push_back(C_map[col]);
                }
                else
                {
                    C->cols.push_back(col); 
                }
                C->data.push_back(block);       
            }
            sums[col].fill(0.0);
        }
        C->rowptr[i+1] = C->cols.size();
    }
    C->nnz = C->cols.size();

    return C;
}


BSRMatrix* BSRMatrix::mult_T(const BSRMatrix* B, int* C_map)
{
    BSCMatrix<b_rows, b_cols>* AT = new BSCMatrix<b_rows, b_cols>(this);

    std::vector<int> next(B->n_cols, -1);
    std::vector<std::array<double, b_size>> sums(B->n_cols);
    for (std::array<double, b_size>& s : sums)
    {
        s.fill(0.0);
    }

    BSRMatrix<b_cols, b_cols>* C = new BSRMatrix<b_cols, b_cols>(n_cols, B->n_cols);

    C->rowptr[0] = 0;
    for (int i = 0; i < AT->n_cols; i++)
    {
        int head = -2;
        int length = 0;
        int col_start = AT->colptr[i];
        int col_end = AT->colptr[i+1];
        for (int j = col_start; j < col_end; j++)
        {
            int row = AT->rows[j];
            int row_start_B = B->rowptr[row];
            int row_end_B = B->rowptr[row+1];
            for (int k = row_start_B; k < row_end_B; k++)
            {
                int col = B->cols[k];
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, b_cols, b_cols, 
                        b_rows, 1.0, AT->data[j].data(), b_rows, B->data[k].data(), b_cols,
                        1.0, sums[col].data(), b_cols);
                if (next[col] == -1)
                {
                    next[col] = head;
                    head = col;
                    length++;
                }
            }
        }
        for (int j = 0; j < length; j++)
        {
            int col = head;
            head = next[col];
            next[col] = -1;
            std::array<double, b_size> block = sums[col];

            if (abs_val(block) > zero_tol)
            {
                if (C_map)
                {
                    C->cols.push_back(C_map[col]);
                }
                else
                {
                    C->cols.push_back(col); 
                }
                C->data.push_back(block);       
            }
            sums[col].fill(0.0);
        }
        C->rowptr[i+1] = C->cols.size();
    }
    C->nnz = C->cols.size();

    delete AT;
    return C;
}

