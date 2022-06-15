#include "matrix_funs_intel_mkl.h"


/* build transpose of matrix : Mt = M^T */
void matrix_build_transpose(mat *Mt, mat *M){
    int i,j;
    for(i=0; i<(M->nrows); i++){
        for(j=0; j<(M->ncols); j++){
            matrix_set_element(Mt,j,i,matrix_get_element(M,i,j)); 
        }
    }
}

void linear_solve_Uxb(mat *A, mat *b) {
    LAPACKE_dtrtrs(LAPACK_ROW_MAJOR, 'U', 'N', 'N',  //unclear
        b->nrows,
        b->ncols, 
        A->d,
        A->ncols,
        b->d,
        b->ncols
    );
}

/*AX = B*/
void linear_solver(mat* A, mat* B, mat* X)
{
    mat *Q = matrix_new(A->nrows, A->ncols);
    mat *R = matrix_new(A->ncols, A->ncols);
    compact_QR_factorization(A, Q, R);
    matrix_transpose_matrix_mult_row(Q, B, X);
    linear_solve_Uxb(R, X); //unclear
    matrix_delete(Q);
    matrix_delete(R);
}

//input A S1 S2 S3 S4
//out Xt = A^T*S1, Y=A*S2, Z = S3^T*A*S4
void matrix_sketch_disk(FILE *A, mat *S1, mat *S2, mat* S3, mat* S4, mat* Xt, mat* Y, mat* Z, int row, int col, int l){
    int row_size = l;
    int read_row_size = row_size;

    double alpha, beta1, beta0;
    int i, j; // count
    int m=row, n=col, k=S1->ncols, s=S4->ncols;
    float *M_f = (float*)malloc(read_row_size*n*sizeof(float));
    double *M = (double*)malloc(read_row_size*n*sizeof(double));;
    // printf("matrix_transpose_matrix_mult_disk is running\n");
    alpha = 1.0;
    beta1 = 1.0;
    beta0 = 0.0;
    // initial Xt=0
    #pragma omp parallel shared(Xt) private(i) 
        {
            #pragma omp parallel for
            for(i=0; i < (Xt->nrows*Xt->ncols); i++){
                Xt->d[i] = 0.0;          
            }
        }
    mat* Z1 = matrix_new(row_size, s);
    mat* Zt = matrix_new(s, s);

    struct timeval start_timeval_1, end_timeval_1;
    struct timeval start_timeval_2, end_timeval_2;
    double sum = 0;
    double time_1, time_2;
    gettimeofday(&start_timeval_1, NULL);

    for (i = 0; i < m; i += row_size){
        if (row_size > (m-i) )
            read_row_size = m-i;
        gettimeofday(&start_timeval_2, NULL);   //time_2
        fread(M_f, sizeof(float), n*read_row_size, A);
        // cblas_dcopy(k, B->d+i*k, 1, g_row, 1); //g_row = g[i];

        #pragma omp parallel shared(M,M_f,n,read_row_size) private(j) 
        {
            #pragma omp parallel for
            for(j=0; j < n * read_row_size; j++){
                M[j] = M_f[j];          
            }
        }
        gettimeofday(&end_timeval_2, NULL);
        sum += get_seconds_frac(start_timeval_2 ,end_timeval_2);
        
        //Xt = A^T*S1
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, k, read_row_size, alpha, M, n, S1->d+i*k, k, beta1, Xt->d, k);
        
        //Y = A*S2
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, read_row_size, k, n, alpha, M, n, S2->d, k, beta0, Y->d+i*k, k);
        
        //Z1 = A*S4
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, read_row_size, s, n, alpha, M, n, S4->d, s, beta0, Z1->d, s);
        
        //Zt = Z1^T*S3
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, s, s, read_row_size, alpha, Z1->d, s, S3->d+i*s, s, beta1, Zt->d, s);

    }
    
    matrix_build_transpose(Z, Zt);
    
    gettimeofday(&end_timeval_1, NULL);
    time_1 = get_seconds_frac(start_timeval_1 ,end_timeval_1);
    
    time_2 = sum;
    printf("Time for reading data file_(fread-time2): %g second\n",time_2);
    printf("Time for matrix_transpose_matrix_mult: %g second\n", time_1);

    free(M_f);
    free(M);
    matrix_delete(Zt);
    matrix_delete(Z1);
}

void SketchySVD_file(char* filename, int m, int n, int r, mat** U, mat** S, mat **V)
{
    int k = 4*r + 1;
    int s = 2*k + 1;
    mat* S1 = matrix_new(m, k);
    mat* S2 = matrix_new(n, k);
    mat* S3 = matrix_new(m, s);
    mat* S4 = matrix_new(n, s);
    initialize_random_matrix(S1);
    initialize_random_matrix(S2);
    initialize_random_matrix(S3);
    initialize_random_matrix(S4);
    
    mat* Xt = matrix_new(n, k);
    mat* Y = matrix_new(m, k);
    mat* Z = matrix_new(s, s);
    FILE* fid = fopen(filename, "rb");
    matrix_sketch_disk(fid, S1, S2, S3, S4, Xt, Y, Z, m, n, r);
    fclose(fid);
    
    QR_factorization_getQ_inplace(Y);
    QR_factorization_getQ_inplace(Xt);
    
    mat* T1 = matrix_new(s, k);
    matrix_transpose_matrix_mult_row(S3, Y, T1);
    mat* T2 = matrix_new(k, s);
    linear_solver(T1, Z, T2);
    mat* T2T = matrix_new(s, k);
    matrix_build_transpose(T2T, T2);
    matrix_delete(T2);
    mat* T3 = matrix_new(s, k);
    matrix_transpose_matrix_mult_row(S4, Xt, T3);
    mat* CT = matrix_new(k, k);
    linear_solver(T3, T2T, CT);
    
    mat* Uk = matrix_new(k, k);
    mat* Sk = matrix_new(k, 1);
    mat* Vk = matrix_new(k, k);
    
    svd_row(CT, Vk, Sk, Uk);
    mat* Ukr = matrix_new(k, r);
    mat* Vkr = matrix_new(k, r);
    int inds[r];
    int i;
    *S = matrix_new(r,1);
    for(i = 0;i<r;i++)
    {
        inds[i] = i;
        (*S)->d[i] = Sk->d[i];
    }
    matrix_get_selected_columns_new(Uk, inds, Ukr);
    *U = matrix_new(m, r);
    matrix_matrix_mult_row(Y, Ukr, *U);
    matrix_get_selected_columns_new(Vk, inds, Vkr);
    *V = matrix_new(n, r);
    matrix_matrix_mult_row(Xt, Vkr, *V);
    matrix_delete(S1);
    matrix_delete(S2);
    matrix_delete(S3);
    matrix_delete(S4);
    matrix_delete(Y);
    matrix_delete(Xt);
    matrix_delete(Z);
    matrix_delete(T1);
    matrix_delete(T2T);
    matrix_delete(T3);
    matrix_delete(CT);
    matrix_delete(Uk);
    matrix_delete(Sk);
    matrix_delete(Vk);
    matrix_delete(Ukr);
    matrix_delete(Vkr);   
}
