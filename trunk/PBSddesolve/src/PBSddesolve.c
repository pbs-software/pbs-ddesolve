#include <R.h>
#include <Rdefines.h>

#include <stdlib.h>
#include "PBSddesolve.h"
#include "ddeq.h"

#define CH_BUF_SIZE 128

int the_test_phase = 0;
int memory_freed = 1; /* when set to 0, freeglobaldata() should be called on re-entering calls to startDDE - only applicable to interupted calculations */

/*===========================================================================*/
/* The function lang5 is part of the C interface for R after version 2.11    */
/* Thanks to Martin Maechler for pointing this out to us.                    */
#ifndef lang5
SEXP lang5(SEXP s, SEXP t, SEXP u, SEXP v, SEXP w)
{
    PROTECT(s);
    s = LCONS(s, list4(t, u, v, w));
    UNPROTECT(1);
    return s;
}
#endif

/*===========================================================================*/
void output(double *s,double t)
{
	/*global_data.vals[0] is for t,
	  and [1..(no_var+1)] are reserved for s[0..no_var] vars
	  */
	int i;
	global_data.vals[0][global_data.vals_ind] = t;
	for( i = 0; i < global_data.no_var; i++ )
		global_data.vals[i+1][global_data.vals_ind] = s[i];
	
	/*ACB hack - call grad to pull out any other data*/
	/* without this call, the other values (returned by grad) won't be calculated exactly at t, but rather at t+/-delta (where delta < step size) which is used during the integration */
	if( global_data.no_otherVars > 0 ) {
		grad(NULL,s,NULL,t);	/* cause a calc exactly at t */

		/* then save the extra variables retuend in the second component of the R grad func */
		for( i = 0; i < global_data.no_otherVars; i++ )
			global_data.vals[1+global_data.no_var+i][global_data.vals_ind] = global_data.tmp_other_vals[i];
	}
	
	global_data.vals_ind++;
	
	if (global_data.vals_ind >= global_data.vals_size) {
		for(i=0;i<(1+global_data.no_var+global_data.no_otherVars);i++) {
			global_data.vals[i] = (double*)realloc(global_data.vals[i], sizeof(double)*2*global_data.vals_size);
			if (global_data.vals[i]==NULL)
				error("memory (re)allocation failed");
		}
		global_data.vals_size *= 2;
	}
}

/*===========================================================================*/
/* cont is zero for fresh run, and 1 for continuation */
void numerics(double *c,int cont, int clear)
{ 
	static double *s;
	double t0, t1, dt, *otimes; /* bjc 2007-05-08*/
	int ns, nsw, nhv, nlag, reset=1, fixstep=0, no_otimes; /* bjc 2007-05-08*/
	static int first=1;
	long hbsize;

	if(clear && first==0) { /* Bobby */
	  free(s); s = NULL; first = 1; return;
	} else if(clear) return;

	ns=global_data.no_var;
	nsw=global_data.nsw;
	nhv=global_data.nhv;
	nlag=global_data.nlag;
	t0=global_data.t0;
	t1=global_data.t1;
	dt=global_data.dt;
	hbsize=global_data.hbsize;
	otimes=global_data.otimes;
	no_otimes=global_data.no_otimes; /* bjc 2007-05-08*/

	if (cont) {
		reset=0;
	} else {
		if (!first) {
			free(s);
			/* first=0; */ /* Bobby */
  		}
		s=(double *)calloc(global_data.no_var,sizeof(double));
		first = 0; /* bobby */
		ddeinitstate(s,c,t0);
	}
	dde(s,c,t0,t1,&dt,global_data.tol,otimes,no_otimes,ns,nsw,nhv,hbsize,nlag,reset,fixstep,0); /* bjc 2007-05-08*/
	global_data.dt=dt;
}

/*===========================================================================*/
void setupglobaldata(int no_vars, int no_otherVars, int no_switch, double *settings, double *otimes, int no_otimes) /* bjc 2007-05-08*/
{ 
	int i;

	global_data.tol=settings[0];
	global_data.t0=settings[1];        /* start time */
	global_data.t1=settings[2];        /* stop time */
  
	global_data.dt=settings[3];        /* initial timestep */
	
	global_data.hbsize=settings[4];    /* how many past values to store for each history variable */
	global_data.no_var=no_vars;
  
	global_data.no_otherVars=no_otherVars;

	global_data.nsw=no_switch;          /* number of switch varaibles */  
	global_data.nhv=no_vars;         /* Number of history (lagged) variables */
	global_data.nlag=1;        /* Number of lag markers per history variable (set to 1 if unsure)*/

	/* enter out times into the data structure */
	global_data.otimes = otimes; /* bjc 2007-05-08: could be NULL*/
	global_data.no_otimes = no_otimes; /* bjc 2007-05-08: >= 0  */

	global_data.vals_size=1000; /* size will grow, this is just initial min size */
	global_data.vals_ind=0;
	global_data.vals = (double**)malloc(sizeof(double*)*(global_data.no_var+no_otherVars+1));
	if (global_data.vals==NULL)
		error("memory allocation failed");
	for(i=0;i<(global_data.no_var+no_otherVars+1);i++) {
		global_data.vals[i]=(double*)malloc(sizeof(double)*global_data.vals_size);
		if (global_data.vals[i]==NULL)
			error("memory allocation failed");
	}
	if (global_data.no_otherVars>0) {
		global_data.tmp_other_vals = (double*)malloc(sizeof(double)*global_data.no_otherVars);
		if (global_data.tmp_other_vals==NULL) {
			error("memory allocation failed");
		}
	} else {
		global_data.tmp_other_vals=NULL;
	}
}

/*===========================================================================*/
void freeglobaldata()
{
	int i;
	if (global_data.vals) {
		for(i=0;i<(global_data.no_var+global_data.no_otherVars+1);i++) {
			free(global_data.vals[i]);
		}
		free(global_data.vals);
		global_data.vals=NULL;
	}
	if (global_data.tmp_other_vals) {
		free(global_data.tmp_other_vals);
		global_data.tmp_other_vals=NULL;
	}
	/* fprintf(stdout, "freed global data\n"); */

	/* Bobby */
	/* necessary to free the memory allocated to static pointers
	   within the following functions */
	istep(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0,0,0,0,NULL,1);
	inithisbuff(0, 0, 0, 1);
	numerics(NULL, 0, 1);
	dde(NULL,NULL,0,0,NULL,0,NULL,0,0,0,0,0,0,1,0,1);
	rk23(NULL,NULL,NULL,NULL,NULL,NULL,0,0,0,1);
	updatehistory(NULL,NULL,NULL,0,1);
}

int testFunc(int no_var, double *test_vars, double t, SEXP *names, PROTECT_INDEX *names_ipx)
{
	SEXP fcall, p1, p2, result, yinit_names, y_names;
	int len, i;
	
	/* argument 1 `t' */
	PROTECT(p1=NEW_NUMERIC(1));
	memcpy(NUMERIC_POINTER(p1), &t, sizeof(double));
    
	/* argument 2 `s' */
	PROTECT(p2=NEW_NUMERIC(no_var));
	memcpy(NUMERIC_POINTER(p2), test_vars, no_var*sizeof(double));

	/* set state names */
	yinit_names = GET_NAMES(r_stuff.yinit);
	PROTECT(y_names = allocVector(STRSXP, no_var));
	if( isNull(yinit_names) == 0 ) {
		for( i = 0; i < no_var; i++ ) {
			SET_STRING_ELT(y_names, i, STRING_ELT(yinit_names, i));
		}
		setAttrib(p2, R_NamesSymbol, y_names);
	}

	/* call R user function */
	if (r_stuff.useParms)
		PROTECT(fcall = lang4(r_stuff.gradFunc, p1, p2, r_stuff.parms));
	else
    	PROTECT(fcall = lang3(r_stuff.gradFunc, p1, p2));
    PROTECT(result = eval(fcall, r_stuff.env));
    
	if (isReal(result)) {
		r_stuff.gradFuncListReturn=0;
		if (LENGTH(result)!=no_var)
			error("func return error: length of vector (%i) does not match that of initial y values (%i)\n", LENGTH(result), no_var);
		UNPROTECT(5);
    	return 0;
    } else if (isVector(result)) {
		r_stuff.gradFuncListReturn=1;
	} else {
		error("func return error: should return numeric vector or list. (got type \"%i\")\n", TYPEOF(result));
	}
	
	p1 = VECTOR_ELT(result, 0);
	
	if (LENGTH(result)!=2 && LENGTH(result)!=1)
		error("func return error: returned list should have length one or two\n", TYPEOF(p1));
	if (!isReal(p1))
		error("func return error: first element of list should be numeric. (got type \"%i\")\n", TYPEOF(p1));
	if (LENGTH(p1)!=no_var)
		error("func return error: length of first element vector (%i) does not match that of initial y values (%i)\n", LENGTH(p1), no_var);

	if (LENGTH(result)==1)
		len = 0;
	else {
		p2 = VECTOR_ELT(result, 1);
		if (!isReal(p2) && !isNull(p2))
			error("func return error: second element of list should be numeric or NULL. (got type \"%i\")\n", TYPEOF(p2));

		if (isNull(p2))
			len = 0;
		else {
			len = LENGTH(p2);
			/* return names - which stays protected even after return */
			REPROTECT((*names) = GET_NAMES(p2), *names_ipx);
		}
	}
	UNPROTECT(5);
	return(len);
}

int testSwitchFunc(int no_var, double *test_vars, double t)
{
	SEXP fcall, p1, p2, result;
	int len;
	
	if (isNull(r_stuff.switchFunc))
		return 0;
	
	/* argument 1 `t' */
	PROTECT(p1=NEW_NUMERIC(1));
	memcpy(NUMERIC_POINTER(p1), &t, sizeof(double));
    
	/* argument 2 `s' */
	PROTECT(p2=NEW_NUMERIC(no_var));
	memcpy(NUMERIC_POINTER(p2), test_vars, no_var*sizeof(double));
	
	if (r_stuff.useParms)
		PROTECT(fcall = lang4(r_stuff.switchFunc, p1, p2, r_stuff.parms));
    else
   		PROTECT(fcall = lang3(r_stuff.switchFunc, p1, p2));
    PROTECT(result = eval(fcall, r_stuff.env));
    
	if (!isReal(result)) {
		error("func return error: should return numeric vector or list. (got type \"%i\")\n", TYPEOF(result));
	}
	
	len = LENGTH(result);
	
	UNPROTECT(4);
	return(len);
}

int testMapFunc(int no_var, double *test_vars, double t, int switch_num)
{
	SEXP fcall, p1, p2, p3, result;
	int len;
	
	if (isNull(r_stuff.mapFunc))
		error("mapFunc is missing, but switchFunc was supplied. both must be defined, or both null");
	
	/* argument 1 `t' */
	PROTECT(p1=NEW_NUMERIC(1));
	memcpy(NUMERIC_POINTER(p1), &t, sizeof(double));
    
	/* argument 2 `s' */
	PROTECT(p2=NEW_NUMERIC(no_var));
	memcpy(NUMERIC_POINTER(p2), test_vars, no_var*sizeof(double));

	/* argument 3 `switchnum' */
	PROTECT(p3=NEW_NUMERIC(1));
	NUMERIC_POINTER(p3)[0] = switch_num;

	if (r_stuff.useParms)
    	PROTECT(fcall = lang5(r_stuff.mapFunc, p1, p2, p3, r_stuff.parms));
    else
		PROTECT(fcall = lang4(r_stuff.mapFunc, p1, p2, p3));	
    PROTECT(result = eval(fcall, r_stuff.env));
    
	if (!isReal(result))
		error("mapFunc return error: should return numeric vector. (got type \"%i\")\n", TYPEOF(result));
		
	len = LENGTH(result);
	if (len != no_var)
		error("mapFunc return error: should return vector of length %i \n", no_var);
	
	UNPROTECT(5);
	return(len);
}

/*===========================================================================*/
SEXP startDDE(SEXP gradFunc, SEXP switchFunc, SEXP mapFunc, SEXP env, SEXP yinit, SEXP parms, SEXP settings, SEXP outtimes)
{
	SEXP list, vect, extra_names, yinit_names, names;
	PROTECT_INDEX extra_names_ipx;
	double *p, *otimes; /* bjc 2007-05-08*/
	int i,j, no_var, no_otherVar, no_switch, no_otimes; /* bjc 2007-05-08*/
	char ch_buf[CH_BUF_SIZE];

	/* free memory on successive calls to startDDE (to prevent leaked memory when R interupts this routine) */
	if( memory_freed == 0 ) {
		memory_freed = 1;
		freeglobaldata();
	}
	
	/* save R global data for later */
	if(!isFunction(gradFunc)) error("\"gradFunc\" must be a function");
	/*TODO check switchFunc is a func, or mark as missing or null*/
	if(!isEnvironment(env)) error("\"env\" should be an environment");
	if(!isNumeric(yinit)) error("\"yinit\" should be a numeric vector");
	if(!isNumeric(settings)) error("\"settings\" should be a numeric vector");
	if(!isNumeric(outtimes) && !isNull(outtimes)) error("\"times\" should be a numeric vector or NULL"); /* bjc 2007-05-08: check times vector*/

	r_stuff.env = env;
	r_stuff.gradFunc = gradFunc;
	r_stuff.switchFunc = switchFunc;
	r_stuff.mapFunc = mapFunc;
	r_stuff.parms = parms;
	r_stuff.yinit = yinit;
	r_stuff.outtimes = outtimes; /* bjc 2007-05-08: add times to R data*/
	
	/* check if supplied function is func(y,t) or func(y,t,parms) */
	list = FORMALS(gradFunc);
	i=0;
	while (list != R_NilValue) {
		i++;
		list = CDR(list);
	}
	if (i!=2 && i!=3) error("\"gradFunc\" must be in the form func(y,t) or func(y,t,parms)");
	r_stuff.useParms = (i==3); /* only use parms if 3 arguments */
	
	no_var = LENGTH(yinit);
	
	/* test the supplied `func` and get number of other vars 
	   ***must set the_test_phase to avoid errors in pastvalue and pastgradient */
	the_test_phase=1;
	PROTECT_WITH_INDEX(extra_names=R_NilValue, &extra_names_ipx);
	no_otherVar=testFunc(no_var, NUMERIC_POINTER(yinit), NUMERIC_POINTER(settings)[1], 
	                     &extra_names, &extra_names_ipx);
	
	/* test switchfunc and get number of results returned to set nsw */
	no_switch = testSwitchFunc(no_var, NUMERIC_POINTER(yinit), NUMERIC_POINTER(settings)[1]);
	
	/*test mapfunc for each nsw and check return val length*/
	for( i = 1; i <= no_switch; i++ )
		testMapFunc(no_var, NUMERIC_POINTER(yinit), NUMERIC_POINTER(settings)[1], i);
	
	the_test_phase=0;
	
	/* print names to see if it works */
	PROTECT(names = allocVector(STRSXP, no_otherVar + no_var + 1));
	yinit_names = GET_NAMES(yinit);
	
	/* fill up names as <- c("time", names(yinit), names(extra_returned_vector)) */
	SET_STRING_ELT(names, 0, mkChar("time"));
	for(i=0;i<no_var;i++) {
		if (isNull(yinit_names)) {
			/* yuck - possible buffer overflow */
			sprintf(ch_buf, "y%i", i+1);
			SET_STRING_ELT(names, i+1, mkChar(ch_buf));
		} else {
			SET_STRING_ELT(names, i+1, STRING_ELT(yinit_names, i));
		}
	}
	for(i=0;i<no_otherVar;i++) {
		if (isNull(extra_names)) {
			/* yuck - possible buffer overflow */
			sprintf(ch_buf, "extra%i", i+1);
			SET_STRING_ELT(names, i+no_var+1, mkChar(ch_buf));
		} else {
			SET_STRING_ELT(names, i+no_var+1, STRING_ELT(extra_names, i));
		}
	}

	
	/* bjc 2007-05-08:  check that the output times are numeric and get 
	   a pointer and length */
	if (!isNumeric(outtimes)) { /* bjc 2007-05-08: if NULL*/
	    otimes = NULL;
	    no_otimes = 0;
	}
	else { /* bjc 2007-05-08: if numeric*/
	    otimes = NUMERIC_POINTER(outtimes); 
	    no_otimes = LENGTH(outtimes);
	}
	
	setupglobaldata(LENGTH(yinit), no_otherVar, no_switch, NUMERIC_POINTER(settings), otimes, no_otimes); /* bjc 2007-05-08*/
	memory_freed = 0;
	
	/* preform dde calculations */
	numerics(NUMERIC_POINTER(yinit), 0, 0);
	
	/* create list which will be base of polyset global_data.frame */
	PROTECT(list=allocVector(VECSXP, global_data.no_var+global_data.no_otherVars+1));
	
	/* room for all data (Y) AND time (T) */
	for(j=0;j<(global_data.no_var+global_data.no_otherVars+1);j++) {
		/* create numeric vector */
		PROTECT(vect=NEW_NUMERIC(global_data.vals_ind));

		/* and fill it up */
		p=NUMERIC_POINTER(vect);
		for(i=0;i<global_data.vals_ind;i++)
			p[i]=global_data.vals[j][i];

		SET_VECTOR_ELT(list, j, vect);
		UNPROTECT(1);
	}
	

	/* Set the names to the global_data.frame */
	/* if names are set as realname.y1 - it's because R is stupid and concatenates name history together:
	wn <- c( a = 1, b = 2 );
	y1 <- 5 * wn[1]; 	y2 <- wn[2] + wn[1];
	vals <- c( y1=y1, y2=y2 ); 	names( vals );
	*/
		
	setAttrib(list, R_NamesSymbol, names);

	UNPROTECT( 3 );
	freeglobaldata();
	memory_freed = 1;
	return list;
}

/*===========================================================================*/
SEXP getPastValue(SEXP t, SEXP markno)
{
	SEXP value;
	int i;

	/* return some value for `func`'s test phase
	   --this won't be used during integration */
	if (the_test_phase)
		return r_stuff.yinit;
	
	if (global_data.vals==NULL) error("pastvalue can only be called from `func` when triggered by dde solver.");
	if (!isNumeric(t)) error("\"t\" should be numeric");
	if (!isInteger(markno)) error("\"markno\" must be an integer");
	if (global_data.hbsize<=0) error("no history buffer was created. dde(...) "
	                          "should be called with hbsize>0");
	if (INTEGER_POINTER(markno)[0] >= global_data.nlag || INTEGER_POINTER(markno)[0] < 0) 
		error("markno is out of bounds and should be in 0..global_data.nlag");
	
	if (NUMERIC_POINTER(t)[0] < global_data.t0 || NUMERIC_POINTER(t)[0] >= global_data.current_t) {
		Rprintf("getvalue error: t=%.5f current integration at t=%.5f\n", NUMERIC_POINTER(t)[0], global_data.current_t);
		error("t is out of bounds and should be >= t0 and < t.\nCalling pastvalue(t) is not allowed.");
	}

	PROTECT(value=NEW_NUMERIC(global_data.no_var));
	for(i=0;i<global_data.no_var;i++) {
		NUMERIC_POINTER(value)[i] = pastvalue(i,
		                                      *NUMERIC_POINTER(t),
		                                      *INTEGER_POINTER(markno));
	}
	UNPROTECT(1);
	return(value);
}

/*===========================================================================*/
SEXP getPastGradient(SEXP t, SEXP markno)
{
	SEXP value;
	int i;

	/* return some value for `func`'s test phase
	   --this won't be used during integration */
	if (the_test_phase)
		return r_stuff.yinit;

	if (global_data.vals==NULL) error("pastgradient can only be called from `func` when triggered by dde solver.");
	if (!isNumeric(t)) error("\"t\" should be numeric");
	if (!isInteger(markno)) error("\"markno\" must be an integer");
	if (global_data.hbsize<=0) error("no history buffer was created. dde(...) "
	                          "should be called with hbsize>0");
	if (INTEGER_POINTER(markno)[0] >= global_data.nlag || INTEGER_POINTER(markno)[0] < 0) 
		error("markno is out of bounds and should be in 0..global_data.nlag");
	
	if (NUMERIC_POINTER(t)[0] < global_data.t0 || NUMERIC_POINTER(t)[0] >= global_data.current_t) 
		error("t is out of bounds and should be >= t0 and < t.\nCalling pastvalue(t) is not allowed.");

	PROTECT(value=NEW_NUMERIC(global_data.no_var));
	for(i=0;i<global_data.no_var;i++) {
		NUMERIC_POINTER(value)[i] = pastgradient(i,
		                                         *NUMERIC_POINTER(t),
		                                         *INTEGER_POINTER(markno));
	}
	UNPROTECT(1);
	return(value);
}

