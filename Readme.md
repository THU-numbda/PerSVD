Programs of PerSVD Algorithm
---
### 1.Main Algorithms

1.matlab/rSVD_fp.m ---- randomized SVD algorithm with fewer passes (Alg. 3 in our paper).

2.matlab/PerSVD_once.m ---- PerSVD with shifted power iteration with fixed shift value (Alg. 4 in our paper).

2.matlab/PerSVD_update.m ---- PerSVD with shifted power iteration, dynamic scheme to update the shift value and pass-efficient scheme (Alg. 5 in our paper).

3.mkl/rsvd.c ---- PerSVD algorithm in our paper which is implemented in C with MKL and OpenMP (Alg. 5 in our paper).

### 2.Experiments for Testing

(1)The program for testing PerSVD with MKL is in "mkl/". The MKL library needs the support of Intel MKL [1]. When all the libraries have been prepared, firstly modified the path of MKL in makefile, and secondly use "make" to produce the executable program "persvd_test". The result of program is an example of testing the dataset Dense1 in size 1000 x 1000.

(2)matlab/PerSVD_test.m is used to test the effectiveness of shifted power iteration. The comparison is between Alg. 1 (basic_rSVD.m)  [2], Alg. 3  (rSVD_fp.m), Alg. 4  (PerSVD_once.m)  and Alg. 5 (PerSVD.m) on Dense1 in size 1000 x 1000.

(3)It should be mentioned that the singular values computed by eigSVD are in ascending order.

### Reference

[1]  Intel oneAPI Math Kernel Library. https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html, 2021. 

[2] N Halko, P. G Martinsson, and J. A Tropp. Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions. Siam Review, 53(2):217-288, 2011. 