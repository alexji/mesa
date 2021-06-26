//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "BlockMat.hpp"
#include "BRelax.hpp"
#include "BILUK.hpp"
#include "BTIF.hpp"
#include "BPKIT.hpp"
#include "HBTMat.hpp"
#include "FGMRES.hpp"
#include "BpResource.hpp"

#if 0
void (*_new_handler)();

static void freeStoreException()
{
    std::cerr << "bp: memory exhausted" << std::endl;
    std::cout << std::endl; // flush stdout
    exit(1);
}
#endif

void condest(int n, GlobalPrecon& M, BlockMat& B)
{
    // A bound of the infinity-norm condition number of the preconditioner
    // equal to infinity-norm of M^-1 e, where e is the vector of all ones.

    double *u = new double[n];
    double *v = new double[n];
    int i;

    for (i=0; i<n; i++)
       u[i] = 1.0;

    M.apply(n, 1, u, n, v, n);

    double m = 0.0;
    for (i=0; i<n; i++)
        if (ABS(v[i]) > m) 
            m = ABS(v[i]);

    // std::cout << "condest: " << m << std::endl;
    printf("condest: %.2e\n", m);

#if 0
    // measure accuracy
    B.mult(n,1,v, n, u, n);
    m = 0.0;
    for (i=0; i<n; i++)
        m += u[i]*u[i];

    m = 1.0/sqrt(m);
    double temp = 1.0/sqrt(n);
    double temp1;

    double inacc = 0.0;
    for (i=0; i<n; i++)
    {
        temp1 = u[i]*m - temp;
        inacc += temp1*temp1;
    }
    inacc = sqrt(inacc);
    printf("inacc: %.2e\n", inacc);
#endif

    delete [] v;
    delete [] u;
}

void set_localprecon(GlobalPrecon *M, char *argv[])
{
    if (!strcmp(argv[0],      "LU"))
             M->localprecon(LP_LU);

    else if (!strcmp(argv[0], "INVERSE"))
             M->localprecon(LP_INVERSE);

    else if (!strcmp(argv[0], "SVD"))
             M->localprecon(LP_SVD,    atof(argv[1]), atof(argv[2]));

    else if (!strcmp(argv[0], "RILUK"))
             M->localprecon(LP_RILUK,  atoi(argv[1]), atof(argv[2]));

    else if (!strcmp(argv[0], "ILUT"))
             M->localprecon(LP_ILUT,   atoi(argv[1]), atof(argv[2]));

    else if (!strcmp(argv[0], "APINV_TRUNC"))
             M->localprecon(LP_APINV_TRUNC, atoi(argv[1]));

    else if (!strcmp(argv[0], "APINV_BANDED"))
             M->localprecon(LP_APINV_BANDED, atoi(argv[1]));

    else if (!strcmp(argv[0], "APINV0"))
             M->localprecon(LP_APINV0);

    else if (!strcmp(argv[0], "APINVS"))
             M->localprecon(LP_APINVS, atoi(argv[1]), atoi(argv[2]));

    else if (!strcmp(argv[0], "DIAG"))
             M->localprecon(LP_DIAG);

    else if (!strcmp(argv[0], "TRIDIAG"))
             M->localprecon(LP_TRIDIAG);

    else if (!strcmp(argv[0], "SOR"))
             M->localprecon(LP_SOR,    atof(argv[1]), atoi(argv[2]));

    else if (!strcmp(argv[0], "SSOR"))
             M->localprecon(LP_SSOR,   atof(argv[1]), atoi(argv[2]));

    else if (!strcmp(argv[0], "GMRES"))
             M->localprecon(LP_GMRES,  atoi(argv[1]), atof(argv[2]));

    else if (!strcmp(argv[0], "DIAGDOM"))
             M->localprecon(LP_DIAGDOM,atof(argv[1]));

    else if (!strcmp(argv[0], "GERSH"))
             M->localprecon(LP_GERSH,  atof(argv[1]));

    else
    {
      std::cerr << "bp: local preconditioner not defined: " << argv[0] << std::endl;
        exit(1);
    }
}

int main(int argc, char *argv[])
{
#if 0
    _new_handler = freeStoreException;
#endif

    if (argc == 1)
    {
   std::cout << "Usage: bp HBfile bsizef global_precon gparam local_precon lparam\n";
       std::cout << "                                             \n";
       std::cout << "global_precon:                               \n";
       std::cout << "       None                                  \n";
       std::cout << "       BJacobi                               \n";
       std::cout << "       BSOR         omega steps              \n";
       std::cout << "       BSSOR        omega steps              \n";
       std::cout << "       BILUK        level                    \n";
       std::cout << "       BTIF                                  \n";
       std::cout << "                                             \n";
       std::cout << "local_precon:                                \n";
       std::cout << "       LU                                    \n";
       std::cout << "       INVERSE                               \n";
       std::cout << "       SVD          rthresh athresh          \n";
       std::cout << "       RILUK        level omega              \n";
       std::cout << "       ILUT         lfil threshold           \n";
       std::cout << "       APINV_TRUNC  semibw                   \n";
       std::cout << "       APINV_BANDED semibw                   \n";
       std::cout << "       APINV0                                \n";
       std::cout << "       APINVS       lfil self_precon         \n";
       std::cout << "       DIAG                                  \n";
       std::cout << "       TRIDIAG                               \n";
       std::cout << "       SOR          omega iterations         \n";
       std::cout << "       SSOR         omega iterations         \n";
       std::cout << "       GMRES        iterations tol           \n";
//     std::cout << "       DIAGDOM                               \n";
//     std::cout << "       GERSH        thresh                   \n";
       exit(0);
    }

    HBTMat H(argv[1]);
    double *rhs = H.get_rhs();
    double *guess = H.get_guess();
    double *exact = H.get_exact();

    int n = H.dimrow();
    int dense_size = BpGetInt("bp.dense_size", 16);
    int blocksize;
    BlockType blocktype;
    BlockMat *B;
    GlobalPrecon *M;

    //
    // Set Block Matrix
    //

    if (isdigit(argv[2][0]))
    {
        blocksize = atoi(argv[2]);
        if (blocksize == 0)
            blocksize = n;

        if (blocksize <= dense_size)
            blocktype = DENSE;
        else
            blocktype = CSR;

	B = new BlockMat(H, blocksize, blocktype);
        H.FreeMat();
    }
    else
    {
        std::cout << "Using block partitioning from file: " << argv[2] << std::endl;
        FILE *fp = fopen(argv[2], "r");
        if (fp == NULL)
	{
            std::cerr << "bp: could not open file: " << argv[2] << std::endl;
	    exit(1);
	}

	int *partit = new int[n+1];
	int i = 0;

        while (fscanf(fp, "%d", &partit[i]) != EOF)
	{
	    // std::cout << i << partit[i] << std::endl;
	    i++;
	}

	if (partit[0] != 0 || partit[i-1] != n)
	{
            std::cerr << "bp: bad partitioning file. " << partit[i-1] << std::endl;
            exit(1);
	}

	// use the first blocksize to determine block type
        if (partit[1] <= dense_size)
            blocktype = DENSE;
        else
            blocktype = CSR;

        // std::cout << "Number of block equations: " << i-1 << std::endl;
	B = new BlockMat(H, i-1, partit, blocktype);
	delete partit;
        H.FreeMat();
	fclose(fp);
    }

    // std::cout << *B << std::endl;

    //
    // Set Global Preconditioner
    //

    if (!strcmp(argv[3], "None"))
    {
        M = new None;
    }
    else if (!strcmp(argv[3], "BJacobi"))
    {
        M = new BJacobi;
        set_localprecon(M, &argv[4]);
        ((BJacobi *)M)->setup(*B);
    }
    else if (!strcmp(argv[3], "BSOR"))
    {
        M = new BSOR;
        set_localprecon(M, &argv[6]);
        ((BSOR *)M)->setup(*B, atof(argv[4]), atoi(argv[5]));
    }
    else if (!strcmp(argv[3], "BSSOR"))
    {
        M = new BSSOR;
        set_localprecon(M, &argv[6]);
        ((BSSOR *)M)->setup(*B, atof(argv[4]), atoi(argv[5]));
    }
    else if (!strcmp(argv[3], "BILUK"))
    {
        M = new BILUK;
        set_localprecon(M, &argv[5]);
        ((BILUK *)M)->setup(*B, atoi(argv[4]));
    }
    else if (!strcmp(argv[3], "BTIF"))
    {
        M = new BTIF;
        set_localprecon(M, &argv[4]);
        ((BTIF *)M)->setup(*B);
    }
    else
    {
        std::cerr << "bp: global preconditioner not defined: " << argv[3] << std::endl;
        exit(1);
    }

    //
    // Iterative Solution
    //

    if (BpGetInt("bp.condest", 0))
        condest(n, *M, *B);

    fgmres f;
    f.solve(*B, guess, rhs, *M);
    f.status(std::cout);

    delete M;
    delete B;

    return 0;
}
