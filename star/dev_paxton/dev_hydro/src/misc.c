#include "paul.h"
#include <string.h>

double get_dA( double );
double get_dV( double , double );
double get_g( struct cell * );
double get_GMr( struct cell * );
double get_moment_arm( double , double );

double mindt( double * , double , double , double , double );

double getmindt( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int nz = theDomain->nz;
   int gE = theDomain->theParList.grav_e_mode;
   int gFl = theDomain->theParList.grav_flag;

   double dt = 1e100;
   int imin = 0;
   int imax = nz;

   int i;
   for( i=imin ; i<imax ; ++i ){
      int im = i-1;
      if( i==0 ) im = 0;
      struct cell * c = theCells+i;
      double dr = c->dr;
      double r = c->riph-.5*dr;
      double wm = theCells[im].wiph;
      double wp = c->wiph;
      double w = .5*(wm+wp);
      double g = 0.0;
      if( gFl && gE == 1 ) g = get_g( c );
      double dt_temp = mindt( c->prim , w , r , g , dr );
      if( dt > dt_temp ) dt = dt_temp;
   }
   dt *= theDomain->theParList.CFL; 
   //MPI_Allreduce( MPI_IN_PLACE , &dt , 1 , MPI_DOUBLE , MPI_MIN , MPI_COMM_WORLD );

   return( dt );
}

void initial( double * , double * );
void cons2prim( double * , double * , double , double );
double get_vr( double * );

void set_wcell( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int mesh_motion = theDomain->theParList.Mesh_Motion;
   int nz = theDomain->nz;
   int bufferzone = 1;
   double Rmax = theCells[nz-1].riph;
   //MPI_Allreduce( MPI_IN_PLACE , &Rmax , 1 , MPI_DOUBLE , MPI_MAX , MPI_COMM_WORLD );
   double Rbuf = 0.8*Rmax;

   int i;
   for( i=0 ; i<nz ; ++i ){
      struct cell * cL = theCells+i;  
      double w = 0.0;
      if( mesh_motion && i<nz-1 ){
         struct cell * cR = theCells+i+1;
         double wL = get_vr( cL->prim );
         double wR = get_vr( cR->prim );
         w = .5*(wL + wR); 
         if( i==0 && 0==0 ) w = 0.5*wR*(cL->riph)/(cR->riph-.5*cR->dr);//*(cR->riph - .5*cR->dr)/(cL->riph);//0.0;//2./3.*wR;
         if( bufferzone && cL->riph > Rbuf ) w *= (Rmax-cL->riph)/(Rmax-Rbuf);//w = 0.0;
         //if( i==0 && 0==0 ) w = 0.0;
      }
      cL->wiph = w;
   }
}

void adjust_RK_cons( struct domain * theDomain , double RK ){

   struct cell * theCells = theDomain->theCells;
   int nz = theDomain->nz;

   int i,q;
   for( i=0 ; i<nz ; ++i ){
      struct cell * c = theCells+i;
      for( q=0 ; q<NUM_Q ; ++q ){
         c->cons[q] = (1.-RK)*c->cons[q] + RK*c->RKcons[q];
      }
   }
}

void move_cells( struct domain * theDomain , double RK , double dt){

   struct cell * theCells = theDomain->theCells;
   int nz = theDomain->nz;
   int i;
   for( i=0 ; i<nz ; ++i ){
      struct cell * c = theCells+i;
      c->riph += c->wiph*dt;
   }

}

void calc_dr( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int nz = theDomain->nz;

   int i;
   for( i=1 ; i<nz ; ++i ){
      int im = i-1;
      double rm = theCells[im].riph;
      double rp = theCells[i ].riph;
      double dr = rp-rm;
      theCells[i].dr = dr;
   }
   theCells[0].dr = theCells[0].riph;

}

void calculate_mass( struct domain * );
void calculate_pot( struct domain * );

void calc_prim( struct domain * theDomain ){


   struct cell * theCells = theDomain->theCells;
   int gE = theDomain->theParList.grav_e_mode;
   if( gE ) calculate_mass( theDomain );
   if( gE == 2 ) calculate_pot( theDomain );

   int nz = theDomain->nz;

   int i;
   #pragma omp for private(i)
   for( i=0 ; i<nz ; ++i ){
      struct cell * c = theCells+i;
      double rp = c->riph;
      double rm = rp-c->dr;
      double dV = get_dV( rp , rm );
//      double g = 0.0;
//      if( gE == 1 ) g = get_g( c );
//      double pot = 0.0;
//      if( gE == 2 ) pot = c->pot;
      double GMr = 0.0;
      if( gE == 3 ) GMr = get_GMr( c );
      cons2prim( c->cons , c->prim , GMr , dV );
   }

}

void plm( struct domain *);
void riemann_flux( struct cell * , struct cell * , double );

void radial_flux( struct domain * theDomain , double dt ){

   int gE = theDomain->theParList.grav_e_mode;
   if( gE ) calculate_mass(theDomain);
   if( gE == 2 ){
      calculate_pot( theDomain );
   }
   struct cell * theCells = theDomain->theCells;
   int nz = theDomain->nz;
   int i;
   plm( theDomain );
   
   #pragma omp for private(i)
   for( i=0 ; i<nz-1 ; ++i ){
      struct cell * cL = theCells+i;
      struct cell * cR = theCells+i+1;
      double r = cL->riph;
      riemann_flux( cL , cR , r );
   }
   
   for( i=0 ; i<nz-1 ; ++i ){
      struct cell * cL = theCells+i;
      struct cell * cR = theCells+i+1;
      double r = cL->riph;
      double dA = get_dA(r); 
      double dAdt = dA*dt;
      int q;
      for( q=0 ; q<NUM_Q ; ++q ){
         double flux = cL->flux_iph[q]*dAdt;
         cL->cons[q] -= flux;
         cR->cons[q] += flux;
      }
   }

}

void source( double * , double * , double , double , double );
void source_alpha( double * , double * , double * , double , double );
void calculate_mass( struct domain * );
void grav_src( struct cell * , double );

void add_source( struct domain * theDomain , double dt ){

   struct cell * theCells = theDomain->theCells;
   int nz = theDomain->nz;
   int i;
   int gravity_flag = theDomain->theParList.grav_flag;
   if( gravity_flag ) calculate_mass( theDomain );
   #pragma omp for private(i)
   for( i=0 ; i<nz ; ++i ){
      struct cell * c = theCells+i;
      double rp = c->riph;
      double rm = rp-c->dr;
      double r = get_moment_arm(rp,rm);
      double dV = get_dV(rp,rm);
      source( c->prim , c->cons , rp , rm , dV*dt );
      int inside = i>0 && i<nz-1;
      int q;
      double grad[NUM_Q];
      for( q=0 ; q<NUM_Q ; ++q ){
         if( inside ){
            struct cell * cp = theCells+i+1;
            struct cell * cm = theCells+i-1;
            double dR = .5*cp->dr + c->dr + .5*cm->dr;
            grad[q] = (cp->prim[q] - cm->prim[q])/dR;
         }else{
            grad[q] = 0.0;
         }
      }
      source_alpha( c->prim , c->cons , grad , r , dV*dt );
      if( gravity_flag ) grav_src( c , dt );
   }   


}


void longandshort( struct domain * theDomain , double * L , double * S , int * iL , int * iS ){ 

   struct cell * theCells = theDomain->theCells;
   int nz = theDomain->nz;
   double rmax = theCells[nz-1].riph;
   double rmin = theCells[0].riph;
   double R0 = theDomain->theParList.LogRadius;
   int nz0 = theDomain->theParList.Num_R;
   double dr0 = rmax/(double)nz0;
   double dx0 = log(rmax/rmin)/nz0;
   int logscale = theDomain->theParList.LogZoning;

   double Long  = 0.0; 
   double Short = 0.0; 
   int iLong  = -1;
   int iShort = -1;

   int imin = 0;
   if( logscale==1 ) imin=1;
   int imax = nz;

   int i;
   for( i=imin ; i<imax ; ++i ){
      struct cell * c = theCells+i;
      double dy = c->dr;
      double dx = dr0;
      if( logscale ) dx = c->riph*dx0;
      if( logscale==2 ) dx = (1./(double)nz0)*( c->riph-rmin + rmax/(rmax/R0-1.) )*log(rmax/R0);
      double l = dy/dx;
      double s = dx/dy;
      if( Long  < l ){ Long  = l; iLong  = i; } 
      if( Short < s ){ Short = s; iShort = i; } 
   }

   *iS = iShort;
   *iL = iLong;
   *S = Short;
   *L = Long;

}


void AMR( struct domain * theDomain ){

   double L,S;
   int iL=0;
   int iS=0;
   longandshort( theDomain , &L , &S , &iL , &iS );

   double MaxShort = theDomain->theParList.MaxShort;
   double MaxLong  = theDomain->theParList.MaxLong;

   struct cell * theCells = theDomain->theCells;
   int nz = theDomain->nz;

   int gE = theDomain->theParList.grav_e_mode;

   if( S > MaxShort ){
      printf("Short = %e #%d of %d\n",S,iS,nz);

      int iSp = iS+1;
      int iSm = iS-1;
      if( iS == nz-1 ) iSp=iS;
      if( iS == 0 ) iSm=0;
      //Possibly shift iS backwards by 1 
      double drL = theCells[iSm].dr;
      double drR = theCells[iSp].dr;
      int imin = 0;
      if( (drL<drR && iSm>imin) || iS==nz-1 ){
         --iS;
         --iSm;
         --iSp;
      }
      struct cell * c  = theCells+iS;
      struct cell * cp = theCells+iSp;

      //Remove Zone at iS+1
      c->dr   += cp->dr;
      c->riph  = cp->riph;
      c->dm   += cp->dm;
      c->miph  = cp->miph;
      int q;
      for( q=0 ; q<NUM_Q ; ++q ){
         c->cons[q]   += cp->cons[q];
         c->RKcons[q] += cp->RKcons[q];
      }
      double rp = c->riph;
      double rm = rp - c->dr;
      double dV = get_dV( rp , rm );
//      double g = 0.0;
//      if( gE == 1 ) g = get_g( c );
//      double pot = 0.0;
//      if( gE == 2 ) pot = c->pot;
      double GMr = 0.0;
      if( gE == 3 ) GMr = get_GMr( c );
      cons2prim( c->cons , c->prim , GMr , dV );
      //Shift Memory
      int blocksize = nz-iSp-1;
      memmove( theCells+iSp , theCells+iSp+1 , blocksize*sizeof(struct cell) );
      theDomain->nz -= 1;
      nz = theDomain->nz;
      theDomain->theCells = (struct cell *) realloc( theCells , nz*sizeof(struct cell) );
      theCells = theDomain->theCells;
      if( iS < iL ) iL--;

   }

   if( L > MaxLong ){

      printf("Long  = %e #%d of %d\n",L,iL,nz);
      theDomain->nz += 1;
      nz = theDomain->nz;
      theDomain->theCells = (struct cell *) realloc( theCells , nz*sizeof(struct cell) );
      theCells = theDomain->theCells;
      int blocksize = nz-iL-1;
      memmove( theCells+iL+1 , theCells+iL , blocksize*sizeof(struct cell) );

      struct cell * c  = theCells+iL;
      struct cell * cp = theCells+iL+1;

      double rp = c->riph;
      double rm = rp - c->dr;
      double r0 = pow( .5*(rp*rp*rp+rm*rm*rm) , 1./3. );
      double dm = .5*c->dm;

      c->riph  = r0;
      c->dr    = r0-rm;
//      cp->riph = rp;
      cp->dr   = rp-r0;

      c->dm    = dm;
      cp->dm   = dm;
//      cp->miph = c->miph;
      c->miph -= dm;

      int q;
      for( q=0 ; q<NUM_Q ; ++q ){
         c->cons[q]    *= .5;
         c->RKcons[q]  *= .5;
         cp->cons[q]   *= .5;
         cp->RKcons[q] *= .5;
      }

      double dV = get_dV( r0 , rm );
//      double g = 0.0;
//      double pot = 0.0;
      double GMr = 0.0;
//      if( gE == 1 ) g = get_g( c );
//      if( gE == 2 ) pot = c->pot;
      if( gE == 3 ) GMr = get_GMr( c );
      cons2prim( c->cons , c->prim , GMr , dV );
      dV = get_dV( rp , r0 );
//      if( gE == 1 ) g = get_g( cp );
//      if( gE == 2 ) pot = cp->pot;
      if( gE == 3 ) GMr = get_GMr( cp );
      cons2prim( cp->cons , cp->prim , GMr , dV );

   }

}


