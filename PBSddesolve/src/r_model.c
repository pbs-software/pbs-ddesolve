#include <math.h>
#include "ddesolve95.h"
#include "ddeq.h"


/***************************************************************************/
/*  Put global variables  here. These should never be written to from      */
/* grad() or switchfunctions(), directly or indirectly.                  	*/
/***************************************************************************/



/***************************************************************************/
/*             Problem specific routines                                   */
/***************************************************************************/

void switchfunctions(double *sw, double *s, double *c, double t)
/* This routine sets the values of the switch functions. When the switch
	functions pass through zero from positive to negative the state variables
	may be reset in function map(). The switch functions should pass smoothly
	through 0 and should never have both value and first derivative zero. The
	same switch must not pass through zero from positive to negative more than
	once in an integration timestep. An example of a switch function is:
						sw[0]=sin(pi*t/30.0)
	which passes through zero every 60 time units. Switches may include state
	variables provided the above conditions are met. Note that to get 'Solver'
	style switches define twice as many switches and let e.g. sw[1]=-sw[0] */
{
	SEXP fcall, p1, p2, result;

	if (isNull(r_stuff.switchFunc))
		return;

	/* argument 1 `t' */
	p1 = PROTECT(NEW_NUMERIC(1)); /* protect 1 */
	memcpy(NUMERIC_POINTER(p1), &t, sizeof(double));

	/* argument 2 `s' */
	p2 = PROTECT(NEW_NUMERIC(global_data.no_var)); /* protect 2 */
	memcpy(NUMERIC_POINTER(p2), s, global_data.no_var*sizeof(double));

	/* call R user function */
	if (r_stuff.useParms)
		fcall = PROTECT(lang4(r_stuff.switchFunc, p1, p2, r_stuff.parms)); /* protect 3 */
	else
		fcall = PROTECT(lang3(r_stuff.switchFunc, p1, p2)); /* protect 3 */
	result = PROTECT(eval(fcall, r_stuff.env)); /* protect 4 */

	/* copy data from R into `sw' */
	memcpy(sw, NUMERIC_POINTER(result), global_data.nsw*sizeof(double));

	UNPROTECT(4);
}

void map(double *s, double *c, double t, int swno)
/* This routine is called whenever one of the switch functions passes through
	zero. 'swno' is the number of the switch function. The state variables
	can be changed discontinuously within this routine. eg:
   if (swno==1)
	  { s[0]=coeff[1]*(s[0]);}
	time and the coefficients should not be changed.
*/
{
	SEXP fcall, p1, p2, p3, result;

	if (isNull(r_stuff.mapFunc))
		return;

	/* argument 1 `t' */
	p1 = PROTECT(NEW_NUMERIC(1)); /* protect 1 */
	memcpy(NUMERIC_POINTER(p1), &t, sizeof(double));

	/* argument 2 `s' */
	p2 = PROTECT(NEW_NUMERIC(global_data.no_var)); /* protect 2 */
	memcpy(NUMERIC_POINTER(p2), s, global_data.no_var*sizeof(double));

	/* argument 3 `switchnum' */
	p3 = PROTECT(NEW_NUMERIC(1)); /* protect 3 */
	NUMERIC_POINTER(p3)[0] = swno + 1; /* use R's index starting at 1 idiology */
	
	/* call R user function */
	if (r_stuff.useParms)
		fcall = PROTECT(lang5(r_stuff.mapFunc, p1, p2, p3, r_stuff.parms)); /* protect 4 */
	else
		fcall = PROTECT(lang4(r_stuff.mapFunc, p1, p2, p3)); /* protect 4 */
	result = PROTECT(eval(fcall, r_stuff.env)); /* protect 5 */

	/* copy returned data from R into `s' */
	memcpy(s, NUMERIC_POINTER(result), global_data.no_var*sizeof(double));

	UNPROTECT(5);
}

void grad(double *g, double *s, double *c, double t)
/* This routine must provide the gradients g for the state variables s.
	So ds[i]/dt=g[i]=fi(s,c,t) where c is the coefficient vector. lagged
   variables may be accessed here using pastvalue(i,x,j) which returns the
   ith (starting at zero) lagged variable at time x, using lag pointer k

   (lag pointers are used by pastvalue to store the history buffer location
    corresponding to a lag in order to save exectution time. For example if
    your code requires lagged varaible 0 at lags T and 2T for each gradient
    calculation then it is efficient to obtain these values using:
    pastvalue(0,t-T,0) and pastvalue(0,t-2*T,1) rather than
    pastvalue(0,t-T,0) and pastvalue(0,t-2*T,0). The latter works, it's just
    slower because more time is spent searching for lagged values)
*/
{ 
	SEXP fcall, p1, p2, result, yinit_names, names;
	int i;

	/* store current t, to prevent calls to pastvalue(t) */
	global_data.current_t = t;

	/* argument 1 `t' */
	p1 = PROTECT(NEW_NUMERIC(1)); /* protect 1 */
	memcpy(NUMERIC_POINTER(p1), &t, sizeof(double));

	/* argument 2 `s' */
	p2 = PROTECT(NEW_NUMERIC(global_data.no_var)); /* protect 2 */
	memcpy(NUMERIC_POINTER(p2), s, global_data.no_var*sizeof(double));

	/* Create the names vector. TODO: do this only once, and not in this function.
	Perhaps the testFunc section would be more appropriate */
	yinit_names = PROTECT(GET_NAMES(r_stuff.yinit)); /* protect 3 */
	names = PROTECT(allocVector(STRSXP, global_data.no_var)); /* protect 4 */
	if( isNull(yinit_names) == 0 ) {
		for( i = 0; i < global_data.no_var; i++ ) {
			SET_STRING_ELT(names, i, STRING_ELT(yinit_names, i));
		}
		setAttrib(p2, R_NamesSymbol, names);
	}

	/* call R user function */
	if (r_stuff.useParms)
		fcall = PROTECT(lang4(r_stuff.gradFunc, p1, p2, r_stuff.parms)); /* protect 5 */
	else
		fcall = PROTECT(lang3(r_stuff.gradFunc, p1, p2)); /* protect 5 */
	result = PROTECT(eval(fcall, r_stuff.env)); /* protect 6 */

	/* copy data from R into `g' */
	/* ACB: `g' can be NULL when called by output(); this is only the case when no_var > 0
	 * this is to compute the other vars at a specific time t (versus t+/-stepsize) */
	if (r_stuff.gradFuncListReturn) {
		p1 = VECTOR_ELT(result, 0);
		if( g != NULL )
			memcpy(g, NUMERIC_POINTER(p1), global_data.no_var*sizeof(double));
		if (global_data.no_otherVars>0) {
			p2 = VECTOR_ELT(result, 1);
			memcpy(global_data.tmp_other_vals, NUMERIC_POINTER(p2), global_data.no_otherVars*sizeof(double));
		}
	} else if( g != NULL ) {
		memcpy(g, NUMERIC_POINTER(result), global_data.no_var*sizeof(double));
	}
	UNPROTECT(6);
}

void storehistory(double *his, double *ghis, double *g, double *s, double *c, double t)
/* This is the routine in which the values of the history variables at time
	t are calculated and put in the array his, along with gradients in ghis,
	using state variables s, gradients of s, g, and coefficients c
	e.g. if the state variable 2 is history variable 0, you would need the line:
	his[0]=s[2];ghis[0]=g[2];
*/
{
	int i;
	for(i=0;i<global_data.no_var;i++) {
		his[i]=s[i];
		ghis[i]=g[i];
	}
}

void statescale(double *scale)
/* In this routine you can set scale factors for error control. For each
   state variable the maximum permisable error will be bounded below by the
   tolerance multiplied by scale[i]. If you don't supply values then zero will
   be used.
	Non-zero scale values are useful for variables that start at zero and
   leave zero without 3rd order continuity. */
{
	scale[0]=0.0;
}


void ddeinitstate(double *s, double *c, double t)
/* initialise state variables and any global constants here, you can use c */
{
	int i;
	for(i=0;i<global_data.no_var;i++)
		s[i]=c[i];
}

