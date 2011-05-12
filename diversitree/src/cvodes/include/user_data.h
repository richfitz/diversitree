/* 
   This is the UserData struct that I require.  It might be nice for
   some models to have more complicated data structures, but I think
   that all models that do that also have more complicated
   everythings.
*/
typedef struct {
  int np;
  int neq;
  realtype *p; /* problem parameters */
} UserData;
