
int dispatch(int iop);

int main( int argc , char * argv[] ){
   int err;
   err = dispatch(0); //start
   if (err) return err;
   err = dispatch(1); //steps
   if (err) return err;
   err = dispatch(2); //finish
   return err;
}

