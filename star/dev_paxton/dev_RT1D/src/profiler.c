
#include "paul.h"

void start_clock( struct domain * theDomain ){
   theDomain->Wallt_init = time(NULL);
}
 
int count_cells( struct domain * theDomain ){
   return(theDomain->nz);
}

void generate_log( struct domain * theDomain ){
   time_t endtime = time(NULL);
   int seconds = (int) (endtime - theDomain->Wallt_init);
   
   int Nc = count_cells( theDomain );
   int Nt = theDomain->count_steps;

   double avgdt = (double)seconds/2./(double)Nc/(double)Nt;

   FILE * logfile = fopen("LOGS/times.log","w");
   if( 1 > 1 ) fprintf(logfile,"es");
   fprintf(logfile,".\n");
   fprintf(logfile,"Total time = %d sec\n",seconds);
   fprintf(logfile,"Number of cells = %d\n",Nc);
   fprintf(logfile,"Number of timesteps = %d (x%d)\n",Nt,2);
   fprintf(logfile,"Megazones per second = %.2e\n",1./(avgdt*1e6));
   fprintf(logfile,"Megazones per CPU second = %.2e\n",1./(avgdt*1e6));
   fprintf(logfile,"Time/zone/step = %.2e microseconds\n",(avgdt*1e6));
   fclose(logfile);
}

