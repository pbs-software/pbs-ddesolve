
#ifndef DDE_HEADER_IN
#define DDE_HEADER_IN 1

typedef struct
{ double **buff,**gbuff,*clock,last_time,first_time;
  long offset,size,no,**lagmarker;
} histype;


#define round(a) ((a)-floor(a) <0.5 ? (int)floor(a):(int) floor(a)+1)


/***************************************************************************/
/*                       Function Prototypes                               */
/***************************************************************************/
#ifdef SUNC

void switchfunctions();
void map();
void grad();
void storehistory();
void ddeinitstate();
void output();

/*********************** Integration routines ******************************/

void rk23();
void inithisbuff();
void updatehistory();
double pastvalue();
double pastgradient();
double zeropos();
double istep();
void dde();

#else
/************************ User Supplied Routines ***************************/
/*void initcons(int *no_vars,int *no_cons);*/
void switchfunctions(double *sw,double *state,double *coeff,double time);
void map(double *state,double *coeff,double time,int swno);
void grad(double *g,double *s,double *c,double t);
void storehistory(double *his,double *ghis,double *g,double *s,double *c,double t);
void ddeinitstate(double *s,double *c,double t);
void output(double *s,double t);
void statescale(double *scale);
/*********************** Integration routines ******************************/

void rk23(double *state,double *newstate,double *g,double *newg,double *error,
	  double *coeff,int ns,double time,double dt, int clear);
void inithisbuff(int nhv,long histsize,int nlag, int clear);
void updatehistory(double *g,double *s,double *c,double t,int clear);
double pastvalue(int i,double t,int markno);
double pastgradient(int i,double t,int markno);
double zeropos(double x1,double x2,double x3,double s1,double s2,double s3);
double istep(double *sw0,double *newsws,double *s0,double *news,double *g,
	     double *newg,double *c,double *err,double t0,double t1,int nsw,
	     int ns,int *flickedswitch, int clear);
void dde(double *s,double *c,double t0,double t1,double *dt,double eps,
	 double *otimes, int no_otimes, int ns,int nsw,int nhv,
	 long hbsize,int nlag,int reset, int fixstep, int clear); /* bjc 2007-05-08*/
#endif
#endif
