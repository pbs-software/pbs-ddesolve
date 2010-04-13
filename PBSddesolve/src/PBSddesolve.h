#ifndef _DDESOLVE95_H_
#define _DDESOLVE95_H_ 1

#include <stdio.h>
#include <R.h>
#include <Rdefines.h>


typedef struct { 
	int no_var, no_otherVars;
	int nhv,nlag,nsw;
	double dt,t0,t1,tol;
	long hbsize;
	char **cname,*initialtext,*initialtitle,**cinfo;
	FILE *file;
	int quit,newrun,cont,*findex,fileno;
	double **vals, *tmp_other_vals;
	int vals_size, vals_ind;
	double current_t;
	double *otimes; /* bjc 2007-05-08*/
	int no_otimes; /* bjc 2007-05-08*/
} globaldatatype;

typedef struct {
	SEXP env;
	SEXP gradFunc;
	SEXP switchFunc;
	SEXP mapFunc;
	SEXP yinit;
	SEXP parms;
	SEXP outtimes; /* bjc 2007-05-08*/
	int useParms;
	int gradFuncListReturn;
} globalRdatatype;


globaldatatype data;
globalRdatatype r_stuff;

SEXP lang5(SEXP s, SEXP t, SEXP u, SEXP v, SEXP w);

#endif /* _DDESOLVE95_H_ */
