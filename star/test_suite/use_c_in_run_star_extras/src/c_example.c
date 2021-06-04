
#include <stdio.h>

int example_c(
      int *iopt, double *res, int *n, double A[], int Ai[])
{
    int i, sz;   
    i = *iopt;
    sz = *n;
    *res = i*i;
    if ( i >= 0 && i < sz) {       
	    printf("i=%d passed to example_c: A[i]=%f Ai[i]=%d\n", i, A[i], Ai[i]);
    } else {
	    printf("Invalid i=%d passed to example_c\n",i);
	    return -1;
    }
    return 0;
}
