//
//                BPKIT 2.0 Block Preconditioner Toolkit
//                  Authors:  E. Chow and M. A. Heroux
//    Copyright (c) 1995-1996  The Regents of the University of Minnesota

#include "DenseMat.hpp"

DenseMat::DenseMat(const DenseMat& A)
{
    nrow = A.nrow; 
    ncol = A.ncol;
    register double *p = a = new double[nrow*ncol];
    register double *q = A.a;
    for (int i=0; i<nrow*ncol; i++)
        *p++ = *q++;
}

void DenseMat::Print(std::ostream& os) const
{
    // check not an implicit inverse
    assert (a != NULL);

    const double *p = a;
    for (int j=0; j<numcol(); j++)
    for (int i=0; i<numrow(); i++)
        os << i+1 << "  " << j+1 << "  " << *p++ << std::endl;
}

// non-member functions

std::ostream& operator << (std::ostream& os, const DenseMat& mat)
{
    // should check not an implicit inverse

    const double *a = &mat(0,0);
    for (int j=0; j<mat.numcol(); j++)
    for (int i=0; i<mat.numrow(); i++)
        os << i+1 << "  " << j+1 << "  " << *a++ << std::endl;

    return os;
}
