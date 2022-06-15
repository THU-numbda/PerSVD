#include "matrix_funs_intel_mkl.h"

/*[U, S, V] = eigSVD(A)*/
void eigSVD(mat* A, mat *U, mat *S, mat *V)
{   
    matrix_transpose_matrix_mult_row(A, A, V);
    LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', V->ncols, V->d, V->ncols, S->d);
    mat *V1 = matrix_new(V->ncols, V->ncols);
    matrix_copy(V1, V);
    MKL_INT i, j;
    #pragma omp parallel shared(V1,S) private(i,j) 
    {
    #pragma omp for 
        for(i=0; i<V1->nrows; i++)
        {
            S->d[i] = sqrt(S->d[i]);
            for(j=0; j<V1->ncols;j++)
            {           
                V1->d[j*V1->ncols+i] /= S->d[i];
            }
        }
    }
    mat *Uc = matrix_new(U->nrows, U->ncols);
    matrix_matrix_mult_row(A, V1, Uc);
    matrix_copy(U, Uc);
    matrix_delete(Uc);
    matrix_delete(V1);
}

void basic_rSVD(char *filename, int m, int n, int k, int l, int p, mat **U, mat **S, mat **V)
{
    int i, j;
    FILE *fid;
    mat *Qt = matrix_new(m, l);
    mat *Q = matrix_new(n, l);
    
    initialize_random_matrix(Q);
    fid = fopen(filename, "rb");
    matrix_matrix_mult_disk(fid, Q, Qt, m, n, k);
    fclose(fid);
    QR_factorization_getQ_inplace(Qt);
    for(i=1;i<=p;i++)
    {
        fid = fopen(filename, "rb");
        matrix_transpose_matrix_mult_disk(fid, Qt, Q, m, n, k);
        fclose(fid);
        fid = fopen(filename, "rb");
        matrix_matrix_mult_disk(fid, Q, Qt, m, n, k);
        fclose(fid);
        QR_factorization_getQ_inplace(Qt);
    }
    fid = fopen(filename, "rb");
    matrix_transpose_matrix_mult_disk(fid, Qt, Q, m, n, k);
    fclose(fid);
    
    mat* UU = matrix_new(l, l);
    mat* VV = matrix_new(n, l);
    mat* SS = matrix_new(l, 1);
    svd_row(Q, VV, SS, UU);
    //puts("1");
    mat* UUk = matrix_new(l, k);
    int inds[k];
    *S = matrix_new(k, 1);
    for(i=0;i<k;i++)
    {
        inds[i] = i;
        (*S)->d[i] = SS->d[i];
    }
    *U = matrix_new(m, k);
    matrix_get_selected_columns_new(UU, inds, UUk);
    matrix_matrix_mult_row(Qt, UUk, *U);
    *V = matrix_new(n, k);
    matrix_get_selected_columns_new(VV, inds, *V);
    matrix_delete(Q);
    matrix_delete(Qt);
    matrix_delete(UU);
    matrix_delete(VV);
    matrix_delete(SS);
    matrix_delete(UUk);
}

void PerSVD(char *filename, int m, int n, int k, int l, int p, mat **U, mat **S, mat **V)
{
    int i, j;
    FILE *fid;
    mat *Qt = matrix_new(m, l);
    mat *Q = matrix_new(n, l);
    
    mat* VV = matrix_new(l, l);
    mat* SS = matrix_new(l, 1);
    
    initialize_random_matrix(Q);
    eigSVD(Q, Q, SS, VV);
    
    mat* D1 = matrix_new(l, l);
    mat* D2 = matrix_new(l, l);
    mat* st = matrix_new(l, 1);
    
    int niter = p+1;
    
    double alpha = 0;
    mat* QQ;
    
    for(i=1;i<=niter;i++)
    {
        QQ = matrix_new(n, l);
        
        
        fid = fopen(filename, "rb");
        matrix_union_matrix_mult_disk_mem(fid, Q, Qt, QQ, m, n, k);
        fclose(fid);
        
        if (i == niter) break;
        
        matrix_transpose_matrix_mult_row(QQ, QQ, D1);
        matrix_transpose_matrix_mult_row(Qt, Qt, D2);
        double alpha1 = alpha;
        double alpha2 = 0;
        double alpha_tol;
        for(j=0;j<100;j++)
        {
            matrix_copy(VV, D1);
            matrix_sub_d(VV, D2, 2.0*alpha1);
            LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', VV->ncols, VV->d, VV->ncols, st->d);
            double stl = sqrt(st->d[0]+alpha1*alpha1);
            if(stl < 1e-10) break;
            if(stl < alpha1) break;
            alpha2 = alpha1;
            alpha1 = (alpha1+stl)/2;
            alpha_tol = (alpha1-alpha2)/alpha1;
            if(alpha_tol < 1e-2) break;
        }
        alpha = alpha1;

        matrix_sub_d(QQ, Q, alpha);
        svd_row(QQ, Q, SS, VV);
                
        SS->d[0] = SS->d[l-1];
        
        if(alpha<SS->d[0]) alpha = (alpha+SS->d[0])/2;
        
        matrix_delete(QQ);
    }
    matrix_copy(Q, QQ);
    matrix_delete(QQ);
    eigSVD(Qt, Qt, SS, VV);
    mat* Sl = matrix_new(l, l);
    for(i=0;i<l;i++)
    {
        Sl->d[i*l+i] = 1.0/SS->d[i];
    }
    mat* SiV = matrix_new(l, l);
    matrix_matrix_mult_row(VV, Sl, SiV);
    mat* B = matrix_new(n, l);
    matrix_matrix_mult_row(Q, SiV, B);
    matrix_delete(SiV);
    eigSVD(B, Q, SS, VV);
    matrix_delete(B);
    int inds[k];
    int s = l - k;
    *S = matrix_new(k, 1);
    for(i=s;i<l;i++)
    {
        inds[i-s] = i;
        (*S)->d[i-s] = SS->d[i]; 
    }
    *U = matrix_new(m, k);
    mat *VV2 = matrix_new(k+s, k);
    matrix_get_selected_columns_new(VV, inds, VV2);
    matrix_delete(VV);
    matrix_matrix_mult_row(Qt, VV2, *U);
    matrix_delete(Qt);
    *V = matrix_new(n, k);
    matrix_get_selected_columns_new(Q, inds, *V);
    matrix_delete(Q);
    matrix_delete(SS);
    matrix_delete(VV2);
}
