
#include "paul.h"

double get_moment_arm( double , double );
double get_dV( double , double );

void output( struct domain * theDomain , char * filestart ){

   struct cell * theCells = theDomain->theCells;
   int nz = theDomain->nz;

   char filename[256];
   sprintf(filename,"%s.dat",filestart);

   FILE * pFile = fopen( filename , "w" );
   fprintf(pFile,"#r           dr           Density      Pressure     Velocity     X            Alpha\n");

   int i,q;
   for( i=0 ; i<nz ; ++i ){
      struct cell * c = theCells+i;
      double rp = c->riph;
      double dr = c->dr; 
      double rm = rp-dr;
      double r  = .5*(rp+rm);//get_moment_arm( rp , rm );
      fprintf(pFile,"%e %e ",r,dr);
      for( q=0 ; q<NUM_Q ; ++q ){
         fprintf(pFile,"%e ",c->prim[q]);
      }
      fprintf(pFile,"%e ",c->miph);
      fprintf(pFile,"\n");
   }
   fclose( pFile );

}
