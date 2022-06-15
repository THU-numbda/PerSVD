#include "matrix_funs_intel_mkl.h"

#include "matrix_funs_intel_mkl.h"

/* parameter setting region */
//input data file

char filename[] = "Dense1.dat"; 
int m = 1000;
int n = 1000;   				// column number of matrix          
const int k = 50;
int l = k + k/2;

//Other settings
BOOL outputflag= true;			//output U, S, V to disk file?
/* end parameter setting */

char expname[50];

void PerSVD_test()
{
    puts("PerSVD begin:\n");
    mat* U1;
    mat* S1;
    mat* V1;
    struct timeval start_timeval, end_timeval;
    
    
    gettimeofday(&start_timeval, NULL);
    PerSVD(filename, m, n, k, l, 2, &U1, &S1, &V1);
    gettimeofday(&end_timeval, NULL);
    
    printf("\nTotal runtime of PerSVD: %f second\n\n", get_seconds_frac(start_timeval,end_timeval));    
	    
    if (outputflag)
            {
                char U_buffer[100];
                char V_buffer[100];
                char S_buffer[100];
                sprintf(U_buffer,"%s_persvd_k=%d_U.dat", expname,k);
                sprintf(V_buffer,"%s_persvd_k=%d_V.dat", expname,k);
                sprintf(S_buffer,"%s_persvd_k=%d_S.dat", expname,k);
                // printf("%s",U_buffer);
                FILE* fid_U = fopen(U_buffer,"w");//mark
                FILE* fid_V = fopen(V_buffer,"w");//mark
                FILE* fid_S = fopen(S_buffer,"w");//mark
                
                fwrite(U1->d, sizeof(double), m*k, fid_U);
                fwrite(V1->d, sizeof(double), n*k, fid_V);
                fwrite(S1->d, sizeof(double), k, fid_S);
                fclose(fid_U);
                fclose(fid_V);
                fclose(fid_S);
            }
                
}

void basic_test()
{
    puts("Basic randomized SVD begin:\n");
    mat* U1;
    mat* S1;
    mat* V1;
    struct timeval start_timeval, end_timeval;
    
    
    gettimeofday(&start_timeval, NULL);
    basic_rSVD(filename, m, n, k, l, 2, &U1, &S1, &V1);
    gettimeofday(&end_timeval, NULL);
    
    printf("\nTotal runtime of basic randomized SVD: %f second\n\n", get_seconds_frac(start_timeval,end_timeval));    
	    
    if (outputflag)
            {
                char U_buffer[100];
                char V_buffer[100];
                char S_buffer[100];
                sprintf(U_buffer,"%s_basic_k=%d_U.dat", expname,k);
                sprintf(V_buffer,"%s_basic_k=%d_V.dat", expname,k);
                sprintf(S_buffer,"%s_basic_k=%d_S.dat", expname,k);
                // printf("%s",U_buffer);
                FILE* fid_U = fopen(U_buffer,"w");//mark
                FILE* fid_V = fopen(V_buffer,"w");//mark
                FILE* fid_S = fopen(S_buffer,"w");//mark
                
                fwrite(U1->d, sizeof(double), m*k, fid_U);
                fwrite(V1->d, sizeof(double), n*k, fid_V);
                fwrite(S1->d, sizeof(double), k, fid_S);
                fclose(fid_U);
                fclose(fid_V);
                fclose(fid_S);
            }
                
}

void sketchy_test()
{
    puts("Tropp's algorithm begin:\n");
    mat* U1;
    mat* S1;
    mat* V1;
    struct timeval start_timeval, end_timeval;
    
    
    gettimeofday(&start_timeval, NULL);
    SketchySVD_file(filename, m, n, k, &U1, &S1, &V1);
    gettimeofday(&end_timeval, NULL);
    
    printf("\nTotal runtime of Tropp's alg.: %f second\n\n", get_seconds_frac(start_timeval,end_timeval));    
	    
    if (outputflag)
            {
                char U_buffer[100];
                char V_buffer[100];
                char S_buffer[100];
                sprintf(U_buffer,"%s_tropp_k=%d_U.dat", expname,k);
                sprintf(V_buffer,"%s_tropp_k=%d_V.dat", expname,k);
                sprintf(S_buffer,"%s_tropp_k=%d_S.dat", expname,k);
                // printf("%s",U_buffer);
                FILE* fid_U = fopen(U_buffer,"w");//mark
                FILE* fid_V = fopen(V_buffer,"w");//mark
                FILE* fid_S = fopen(S_buffer,"w");//mark
                
                fwrite(U1->d, sizeof(double), m*k, fid_U);
                fwrite(V1->d, sizeof(double), n*k, fid_V);
                fwrite(S1->d, sizeof(double), k, fid_S);
                fclose(fid_U);
                fclose(fid_V);
                fclose(fid_S);
            }
                
}

int main() {
	char * ptr1, *ptr2;
	int len, i;
	ptr1= strrchr(filename, '/');
	ptr2= strchr(filename, '.');
	if(ptr1== NULL){
	 	i=0;
		while(filename[i]!='.'){
			expname[i]= filename[i];
			i++;
		}
		expname[i]='\0';
	}
	else{
		len= ptr2-ptr1-1;
		for(i=0; i<len; i++){
			expname[i]= ptr1[i+1];
		}
		expname[i]='\0';
	}

	printf("***************************************************************\n");
	printf("Test PeSVD and other counterparts.\n");
	printf("***************************************************************\n\n");
	printf("Input matrix file: %s\n\n", filename);

    PerSVD_test();
    basic_test();
    sketchy_test();
    return 0;
    
}
