// Copyright (c) 2015-2017, RAPtor Developer Team
// License: Simplified BSD, http://opensource.org/licenses/BSD-2-Clause
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>

#include "raptor-sparse.hpp"

#ifdef USING_MFEM
#include "raptor/external/mfem_wrapper.hpp"
#endif

using namespace raptor;


void mult(CSRMatrix* A, double* x, double* b, double alpha, double beta)
{
    int start, end;
    double val;
    for (int i = 0; i < A->n_rows; i++)
    {
        start = A->idx1[i];
        end = A->idx1[i+1];
        val = 0;
        for (int j = start; j < end; j++)
            val += A->vals[j] * x[A->idx2[j]];
        b[i] = alpha * val + beta * b[i];
    }
}

void mult_append(CSCMatrix* A, double* x, double* b)
{
    int start, end, row;

    for (int col = 0; col < A->n_cols; col++)
    {
        start = A->idx1[col];
        end = A->idx1[col+1];
        for (int j = start; j < end; j++)
        {
            row = A->idx2[j];
            b[row] += A->vals[j] * x[col];
        }
    }
}


void init_sends(ParCSRMatrix* A, double* values)
{
    int proc, start, end, idx;

    // Initialize Sends
    for (int i = 0; i < A->comm->send_data->num_msgs; i++)
    {
        proc = A->comm->send_data->procs[i];
        start = A->comm->send_data->indptr[i];
        end = A->comm->send_data->indptr[i+1];
        for (int j = start; j < end; j++)
        {
            idx = A->comm->send_data->indices[j];
            A->comm->send_data->buffer[j] = values[idx];
        }
        MPI_Isend(&(A->comm->send_data->buffer[start]), (end - start), MPI_DOUBLE, 
                proc, A->comm->key, A->comm->mpi_comm, &(A->comm->send_data->requests[i]));
   }
}

void init_recvs(ParCSRMatrix* A)
{
    int proc, start, end;

    for (int i = 0; i < A->comm->recv_data->num_msgs; i++)
    {
        proc = A->comm->recv_data->procs[i];
        start = A->comm->recv_data->indptr[i];
        end = A->comm->recv_data->indptr[i+1];
        MPI_Irecv(&(A->comm->recv_data->buffer[start]), (end - start), MPI_DOUBLE, 
            proc, A->comm->key, A->comm->mpi_comm, &(A->comm->recv_data->requests[i]));
    }
}

void wait_sends(ParCSRMatrix* A)
{
    if (A->comm->send_data->num_msgs)
        MPI_Waitall(A->comm->send_data->num_msgs, 
                A->comm->send_data->requests.data(),
                MPI_STATUSES_IGNORE);
}

void wait_recvs(ParCSRMatrix* A)
{
    if (A->comm->recv_data->num_msgs)
        MPI_Waitall(A->comm->recv_data->num_msgs,
                A->comm->recv_data->requests.data(),
                MPI_STATUSES_IGNORE);
}

void standard_spmv(ParCSRMatrix* A, ParVector& x, ParVector& b)
{
    CSRMatrix* A_on = (CSRMatrix*)A->on_proc;
    CSRMatrix* A_off = (CSRMatrix*)A->off_proc;

    init_recvs(A);
    init_sends(A, x.local.data());
    mult(A_on, x.local.data(), b.local.data(), 1.0, 0.0);
    wait_recvs(A);
    mult(A_off, A->comm->recv_data->buffer.data(), b.local.data(), 1.0, 1.0);
    wait_sends(A);
}

void early_comp_spmv(ParCSRMatrix* A, ParVector& x, ParVector& b)
{
    CSRMatrix* A_on = (CSRMatrix*)A->on_proc;
    CSRMatrix* A_off = (CSRMatrix*)A->off_proc;

    mult(A_on, x.local.data(), b.local.data(), 1.0, 0.0);
    init_recvs(A);
    init_sends(A, x.local.data());
    wait_recvs(A);
    mult(A_off, A->comm->recv_data->buffer.data(), b.local.data(), 1.0, 1.0);
    wait_sends(A);
}

void CSC_spmv(ParCSRMatrix* A, ParVector& x, ParVector& b, CSCMatrix* A_off_csc)
{
    CSRMatrix* A_on = (CSRMatrix*)A->on_proc;

    mult(A_on, x.local.data(), b.local.data(), 1.0, 0.0);
    init_recvs(A);
    init_sends(A, x.local.data());
    wait_recvs(A);
    mult_append(A_off_csc, A->comm->recv_data->buffer.data(), b.local.data());
    wait_sends(A);
}

void CSC_early_recv_spmv(ParCSRMatrix* A, ParVector& x, ParVector& b, CSCMatrix* A_off_csc)
{
    CSRMatrix* A_on = (CSRMatrix*)A->on_proc;

    int proc, start, end, idx;
    int col_start, col_end, row;
    MPI_Status status;

    mult(A_on, x.local.data(), b.local.data(), 1.0, 0.0);

    init_sends(A, x.local.data());

    for (int i = 0; i < A->comm->recv_data->num_msgs; i++)
    {
        /***************************
        **** Recv a Single Message
        ***************************/
        MPI_Probe(MPI_ANY_SOURCE, A->comm->key, A->comm->mpi_comm, &status);
        proc = status.MPI_SOURCE;

        // Find index of proc in recv comm
        std::vector<int>::iterator it = std::find(A->comm->recv_data->procs.begin(),
                A->comm->recv_data->procs.end(), proc);
        idx = it - A->comm->recv_data->procs.begin();

        // Get start,end of recv
        start = A->comm->recv_data->indptr[idx];
        end = A->comm->recv_data->indptr[idx+1];

        MPI_Recv(&(A->comm->recv_data->buffer[start]), (end - start), MPI_DOUBLE,
                proc, A->comm->key, A->comm->mpi_comm, &status);

        /***************************
        **** Multiply this message
        ***************************/
        for (int col = start; col < end; col++)
        {
            col_start = A_off_csc->idx1[col];
            col_end = A_off_csc->idx1[col+1];
            for (int j = col_start; j < col_end; j++)
            {
                row = A_off_csc->idx2[j];
                b.local[row] += A_off_csc->vals[j] * A->comm->recv_data->buffer[col];
            }
        }
    }

    wait_sends(A);
}





using namespace raptor;
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int dim=0;
    int n = 5;
    int system = 0;

    if (argc > 1)
    {
        system = atoi(argv[1]);
    }

    ParCSRMatrix* A=NULL;
    ParVector x;
    ParVector b;

    double t0, tfinal;

    int cache_len = 10000;

    std::vector<double> cache_array(cache_len);

    if (system < 2)
    {
        double* stencil = NULL;
        std::vector<int> grid;
        if (argc > 2)
        {
            n = atoi(argv[2]);
        }

        if (system == 0)
        {
            dim = 3;
            grid.resize(dim, n);
            stencil = laplace_stencil_27pt();
        }
        else if (system == 1)
        {
            dim = 2;
            grid.resize(dim, n);
            double eps = 0.001;
            double theta = M_PI/8.0;
            if (argc > 3)
            {
                eps = atof(argv[3]);
                if (argc > 4)
                {
                    theta = atof(argv[4]);
                }
            }
            stencil = diffusion_stencil_2d(eps, theta);
        }
        A = par_stencil_grid(stencil, grid.data(), dim);
        delete[] stencil;
    }
#ifdef USING_MFEM
    if (system < 2)
    {
        double* stencil = NULL;
        std::vector<int> grid;
        if (argc > 2)
        {
            n = atoi(argv[2]);
        }

        if (system == 0)
        {
            dim = 3;
            grid.resize(dim, n);
            stencil = laplace_stencil_27pt();
        }
        else if (system == 1)
        {
            dim = 2;
            grid.resize(dim, n);
            double eps = 0.001;
            double theta = M_PI/4.0;
            if (argc > 3)
            {
                eps = atof(argv[3]);
                if (argc > 4)
                {
                    theta = atof(argv[4]);
                }
            }
            stencil = diffusion_stencil_2d(eps, theta);
        }
        A = par_stencil_grid(stencil, grid.data(), dim);
        delete[] stencil;
    }
#endif
    else if (system == 3)
    {
        const char* file = "../../test_data/rss_A0.pm";
        if (argc > 2)
        {
            file = argv[2];
        }
        A = readParMatrix(file);
    }

    if (system != 2)
    {
        x = ParVector(A->global_num_cols, A->on_proc_num_cols);
        b = ParVector(A->global_num_rows, A->local_num_rows);
        x.set_const_value(1.0);
    }

    ParVector b_new = ParVector(A->global_num_rows, A->local_num_rows);
    CSCMatrix* A_off_csc = A->off_proc->to_CSC();

    /************************
    **** Correctness Checks
    ************************/
    A->mult(x, b); // Standard RAPtor SpMV

    standard_spmv(A, x, b_new);
    for (int i = 0; i < A->local_num_rows; i++)
        if (fabs(b[i] - b_new[i]) > 1e-06)
            printf("Difference in RESULTS!\n");

    early_comp_spmv(A, x, b_new);
    for (int i = 0; i < A->local_num_rows; i++)
        if (fabs(b[i] - b_new[i]) > 1e-06)
            printf("Difference in RESULTS!\n");
    
    CSC_spmv(A, x, b_new, A_off_csc);
    for (int i = 0; i < A->local_num_rows; i++)
        if (fabs(b[i] - b_new[i]) > 1e-06)
            printf("Difference in RESULTS!\n");

    CSC_early_recv_spmv(A, x, b_new, A_off_csc);
    for (int i = 0; i < A->local_num_rows; i++)
        if (fabs(b[i] - b_new[i]) > 1e-06)
            printf("Difference in RESULTS!\n");

    int n_iters;
    MPI_Barrier(MPI_COMM_WORLD);

    /************************
    **** Time RAPtor SpMV
    ************************/
    // Warm-up
    A->mult(x, b);

    // Calculate iterations
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
    A->mult(x, b);
    tfinal = MPI_Wtime() - t0;
    MPI_Allreduce(&tfinal, &t0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    n_iters = 1.0 / t0 + 1;

    for (int test = 0; test < 5; test++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        t0 = MPI_Wtime();
        for (int i = 0; i < n_iters; i++)
        {
            A->mult(x, b);
        }
        tfinal = (MPI_Wtime() - t0) / n_iters;
        MPI_Reduce(&tfinal, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0) printf("RAPtor SpMV Time: %e\n", t0);
    }

    /************************
    **** Time Standard SpMV
    ************************/
    // Warm-up
    standard_spmv(A, x, b_new);

    // Calculate iterations
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
    standard_spmv(A, x, b_new);
    tfinal = MPI_Wtime() - t0;
    MPI_Allreduce(&tfinal, &t0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    n_iters = 1.0 / t0 + 1;

    for (int test = 0; test < 5; test++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        t0 = MPI_Wtime();
        for (int i = 0; i < n_iters; i++)
        {
            standard_spmv(A, x, b_new);
        }
        tfinal = (MPI_Wtime() - t0) / n_iters;
        MPI_Reduce(&tfinal, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0) printf("Standard SpMV Time: %e\n", t0);
    }

    /************************
    **** Time Early Computation SpMV
    ************************/
    // Warm-up
    early_comp_spmv(A, x, b_new);

    // Calculate iterations
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
    early_comp_spmv(A, x, b_new);
    tfinal = MPI_Wtime() - t0;
    MPI_Allreduce(&tfinal, &t0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    n_iters = 1.0 / t0 + 1;

    for (int test = 0; test < 5; test++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        t0 = MPI_Wtime();
        for (int i = 0; i < n_iters; i++)
        {
            early_comp_spmv(A, x, b_new);
        }
        tfinal = (MPI_Wtime() - t0) / n_iters;
        MPI_Reduce(&tfinal, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0) printf("Early Comp SpMV Time: %e\n", t0);
    }

    /************************
    **** Time CSC SpMV
    ************************/
    // Warm-up
    CSC_spmv(A, x, b_new, A_off_csc);

    // Calculate iterations
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
    CSC_spmv(A, x, b_new, A_off_csc);
    tfinal = MPI_Wtime() - t0;
    MPI_Allreduce(&tfinal, &t0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    n_iters = 1.0 / t0 + 1;

    for (int test = 0; test < 5; test++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        t0 = MPI_Wtime();
        for (int i = 0; i < n_iters; i++)
        {
            CSC_spmv(A, x, b_new, A_off_csc);
        }
        tfinal = (MPI_Wtime() - t0) / n_iters;
        MPI_Reduce(&tfinal, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0) printf("CSC SpMV Time: %e\n", t0);
    }


    /************************
    **** Time Early Recv CSC SpMV
    ************************/
    // Warm-up
    CSC_early_recv_spmv(A, x, b_new, A_off_csc);

    // Calculate iterations
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
    CSC_early_recv_spmv(A, x, b_new, A_off_csc);
    tfinal = MPI_Wtime() - t0;
    MPI_Allreduce(&tfinal, &t0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    n_iters = 1.0 / t0 + 1;

    for (int test = 0; test < 5; test++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        t0 = MPI_Wtime();
        for (int i = 0; i < n_iters; i++)
        {
            CSC_early_recv_spmv(A, x, b_new, A_off_csc);
        }
        tfinal = (MPI_Wtime() - t0) / n_iters;
        MPI_Reduce(&tfinal, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0) printf("Early Recv CSC SpMV Time: %e\n", t0);
    }

    delete A_off_csc;
    delete A;

    MPI_Finalize();

    return 0;
}

