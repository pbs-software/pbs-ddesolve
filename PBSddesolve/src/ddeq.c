/***************************************************************************/
/* This is a new dde solver, that attempts to get around the problem of    */
/* interpolating lagged variables at a lower order than the integration    */
/* scheme. The methods borrow heavily from:                                */
/* [1] C.A.H. Paul 1992 Appl.Numer.Math. 9:403-414                         */
/* [2] D.J. Higham 1993 Appl.Numer.Math. 12:315-330                        */
/* [3] Press et al. 1992 Numerical recipes in C C.U.P. chap. 16            */
/* The integration scheme used is an emdedded RK23 pair due to Fehlberg    */
/* reported in:                                                            */
/* [4]E.Hairer, S.P.Norsett & G.Wanner 1987, Solving Ordinary differential */
/*                   Equations I. springer-Verlag Berlin. p170 RKF2(3)B    */
/* How to derive better schemes is described in:                           */
/* [5] J.C. Butcher 1987 The Numerical Analysis of Ordinary Differential   */
/* Equations. John Wiley & sons, Chichester.                               */
/* Interpolation of lagged variables is by cubic hermite polynomials. This */
/* necessitates storing values and gradients of lagged variables. Some     */
/* models have to be re-cast a bit to do this, but you can always make up  */
/* a lagged quantity from the lagged state variables and the gradients of  */
/* the latter are obviously easy to find. [2] Shows why this effort is     */
/* required.                                                               */
/* Lags of less than the timestep are dealt with by extrapolating the      */
/* cubic polynomial for the last stored step, beyond that time step [1].   */
/* Switches are also dealt with. This means that lagged variables are      */
/* stored twice when a switch occurs: once before and once after. However, */
/* switch tracking is not carried out, so the step length may at times     */
/* reduce as switches in lagged variables are crossed, yielding derivative */
/* discontinuities.                                                        */
/* Stepping logic is from [3].                         		           */
/* Note that the code has no scope for specifying initial histories. This  */
/* could be changed, but initial history problems are seriously unpleasant.*/
/* NOTE: code not yet optimised for minimum no. of gradient evaluations     */
/***************************************************************************/

#define ANSI
/*#include "inverse.h"*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>


#define ERRCON 1.89e-4
#define SWITCHES

/* The following macros are for inline calculation of cubic hermite
   polynomials (HERMITE) and their gradients (GHERMITE)
   The horrible variable names are to (hopefully) ensure that the names don't
   clash with anything else - find and replace with something nicer if you
   need to check the macro. Note this can't be translated to a language
   without left to right associativity.
   Tested 4/10/95 (HERMITE), 20/11/95 (GHERMITE) 
   Modified 8/1/96 - left to right evaluation sequence IS NOT standard!! */


#define HERMITE(res,x0,x1,y0,y1,g0,g1,x) \
	HeRmItE_xx0=x-x0;HeRmItE_xx1=x-x1;HeRmItE_xx12=HeRmItE_xx1*HeRmItE_xx1;\
HeRmItE_xx02=HeRmItE_xx0*HeRmItE_xx0;\
if ((HeRmItE_h=x1-x0))\
res=(( (g0)*HeRmItE_xx0*HeRmItE_xx12 + (g1)*HeRmItE_xx1*HeRmItE_xx02 \
			+((y0)*(2.0*HeRmItE_xx0+(HeRmItE_h))*HeRmItE_xx12- \
				(y1)*(2.0*HeRmItE_xx1-HeRmItE_h)*HeRmItE_xx02)/HeRmItE_h)/ \
		(HeRmItE_h*HeRmItE_h)); else res=y0


double HeRmItE_h,HeRmItE_xx1,HeRmItE_xx12,HeRmItE_xx02,HeRmItE_xx0;


#define GHERMITE(res,x0,x1,y0,y1,g0,g1,x) \
	HeRmItE_xx0=x-x0;HeRmItE_xx1=x-x1;HeRmItE_xx12=HeRmItE_xx1*HeRmItE_xx1;\
HeRmItE_xx02=HeRmItE_xx0*HeRmItE_xx0;\
if ((HeRmItE_h=x1-x0))\
res=(( (g0)*(HeRmItE_xx12+2.0*HeRmItE_xx0*HeRmItE_xx1) \
			+ (g1)*(HeRmItE_xx02+2.0*HeRmItE_xx0*HeRmItE_xx1) \
			+ ((y0)*2.0*HeRmItE_xx1*(2.0*HeRmItE_xx0+HeRmItE_h + HeRmItE_xx1) \
				- (y1)*2.0*HeRmItE_xx0*(2.0*HeRmItE_xx1-HeRmItE_h + HeRmItE_xx0))/HeRmItE_h )\
		/ (HeRmItE_h*HeRmItE_h) ); else res=g0;

/***************** end of definitions for hermite macros *******************/


/***************************************************************************/
/*                       Global variables                                  */
/***************************************************************************/
#include "ddeq.h"

long accepted=0L,rejected=0L;
histype history;
/***************************************************************************/
/*             Routines that are not problem specific                      */
/***************************************************************************/


void rk23(double *state, double *newstate, double *g, double *newg,
	  double *errors, double *coeff, int ns, double time, double dt,
	  int clear)
	/* Takes a single integration step from time to time+dt using a 3rd order
	   embedded Runge-Kutta Fehlberg method:
	   E.Hairer, S.P.Norsett & G.Wanner 1987, Solving Ordinary differential
	   Equations I. springer-Verlag Berlin. p170 RKF2(3)B
	   The routine returns an estimated errors vector for adaptive timestepping.
	   The gradient of the state variables is to be given in function grad().
	   The routine uses the lower order scheme for updating,
	   fortunately Fehlberg optimised the coefficients for the lower order
	   scheme..... 4/10/95.
NOTE: not yet optimised for minimum gradient evaluations - see original
table of coeffs. Partially optimised 9/10/95 Only valid for ci=b4i!
Takes gradient at time in g, puts gradient at time+dt in newg - these can
be the same pointer/array */


{ static int first=1,oldns=-1;
	static double *k1,*k2,*k3,*k4,
		      /* Embedded RKF table - coded this way to save addressing time */
		      a2=0.25,  a3=27.0/40.0,
		      b21= 0.25,
		      b31=-189.0/800.0,  b32= 729.0/800.0,
		      b41= 214.0/891.0,  b42= 1.0/33.0,     b43= 650.0/891.0,
		      /* c1= 214.0/891.0,   c2= 1.0/33.0,      c3= 650.0/891.0,*/
		      cc1= 533.0/2106.0, cc3= 800.0/1053.0, cc4=-1.0/78.0;
	int i;
       
	/* Bobby */
	if(clear && first == 0) {
		free(k2);free(k3);free(k4);
		first = 1;
		return;
	} else if(clear) return;

	if ((first)||(oldns!=ns))
	{ if (!first)
		{ free(k2);free(k3);free(k4);}
		oldns=ns;first=0;
		if (ns>0)
		{ k2=(double *)calloc(ns,sizeof(double));
			k3=(double *)calloc(ns,sizeof(double));
			k4=(double *)calloc(ns,sizeof(double));
		}
	}
	k1=g;
	for (i=0;i<ns;i++) newstate[i]=state[i]+(k1[i]*b21)*dt;
	grad(k2,newstate,coeff,time+dt*a2);
	for (i=0;i<ns;i++) newstate[i]=state[i]+(k1[i]*b31+k2[i]*b32)*dt;
	grad(k3,newstate,coeff,time+dt*a3);
	for (i=0;i<ns;i++)
		newstate[i]=state[i]+(k1[i]*b41+k2[i]*b42+k3[i]*b43)*dt;

	grad(k4,newstate,coeff,time+dt);
	for (i=0;i<ns;i++)
	{ newg[i]=k4[i];
		errors[i]=state[i]+(cc1*k1[i]+cc3*k3[i]+cc4*k4[i])*dt-newstate[i];
	}
}


void inithisbuff(int nhv, long histsize, int nlag, int clear)
	/* sets up the global structure "history" and
	   sets the global long integer history.offset to zero
	   4/10/95, if it's been called before it clears up the old version first 23/11/95*/

{ static int oldnhv=0;
	int i;

	if(clear && oldnhv) { /* Bobby */
		for (i=0;i<oldnhv;i++) {
			free(history.buff[i]);
			free(history.lagmarker[i]);
			free(history.gbuff[i]);
		}
		if (oldnhv) {
			free(history.lagmarker);
			free(history.clock);
			free(history.buff);
			free(history.gbuff);
		}
		oldnhv = 0;
		return;
	}
	if(clear)
		return;

	for (i=0;i<oldnhv;i++)
	{ free(history.buff[i]);
		free(history.lagmarker[i]);
		free(history.gbuff[i]);
	}
	if (oldnhv) 
	{ free(history.lagmarker);
		free(history.clock);
		free(history.buff);
		free(history.gbuff);
	}
	oldnhv=nhv;
	if (!nhv) return;
	history.no=(long)nhv;
	history.size=histsize;
	history.lagmarker=(long **)calloc((size_t)nhv,sizeof(long *));
	for (i=0;i<nhv;i++)
		history.lagmarker[i]=(long *)calloc((size_t)nlag,sizeof(long));
	history.clock=(double *)calloc((size_t)history.size,sizeof(double));
	history.buff=(double **)calloc((size_t)nhv,sizeof(double *));
	for (i=0;i<nhv;i++)
		history.buff[i]=(double *)calloc((size_t)history.size,sizeof(double));
	history.gbuff=(double **)calloc((size_t)nhv,sizeof(double *));
	for (i=0;i<nhv;i++)
		history.gbuff[i]=(double *)calloc((size_t)history.size,sizeof(double));
	if (!history.gbuff[nhv-1])
	{ error("History buffer too long");
	}
	history.offset= -1L;
}

void updatehistory(double *g, double *s, double *c, double t, int clear)
	/* updates the history record by calling the storehistory() moving the
	   offset and updating and recording the time 4/10/95*/

{ static int first=1, oldhno=-1L;
	static double *his,*ghis;
	int i;
	
	if(clear && first==0) { free(his); free(ghis); first=1; return; }
	else if(clear) return;

	if (! history.no) return;
	if ((first)||(oldhno!=history.no))
	{ if (!first) { free(his);free(ghis);}
		first=0;his=(double *)calloc((size_t)history.no,sizeof(double));
		ghis=(double *)calloc((size_t)history.no,sizeof(double));
		oldhno=(int)history.no;
		history.first_time=t;
		history.offset=-1L;
	}
	storehistory(his,ghis,g,s,c,t);
	history.last_time=t;
	history.offset++;
	if (history.offset==history.size) {
		history.offset=0L;
	}
	history.clock[history.offset]=t;
	for (i=0;i<history.no;i++)
	{ history.buff[i][history.offset]=his[i];
		history.gbuff[i][history.offset]=ghis[i];
	}
}

double pastgradient(int i, double t, int markno)
	/* Interogates the history ringbuffers. Note that there is a fair amount of
	   assignment of one variable to another at the start: this is to save
	   on address calculation and speed up the routine. 4/10/95 (copy from
	   pastvalue 20/11/95)*/

{
	long k1,k,offset,offsetplus,size;
	double res,*y,*g,*x,x0,x1;
	y=history.buff[i];g=history.gbuff[i];
	x=history.clock; /*local pointers improve address efficiency*/
	offset=history.offset;size=history.size;
	offsetplus=offset+1L; if (offsetplus==size) offsetplus=0L;
	k=history.lagmarker[i][markno];
	k1=k+1L;if ((k1>=size)||(k1<0L)) k1=0L;
	while ((x[k1]<t)&&(k1!=offset)) { 
		k1++;if (k1==size) k1=0L;
	}
	if (k1==0L) k=size-1L; else k=k1-1L;
	while ((x[k]>t)&&(k!=offsetplus)) { 
		if (k==0L) k=size-1L; else k--;
	}
	k1=k+1L;if (k1==size) k1=0L;
	if (t<x[k]) {
		Rprintf("\nERROR: lag for variable %d too large at %g\n",i,history.last_time-t);
		error("Lag too large for history buffer - try increasing the value of 'hbsize'");
	}
	x0=x[k];x1=x[k1];
#ifdef SWITCHES  /* some code for not extrapolating through a switch */
	if ((t>x[k1])&&(x[k]==x[k1])) /* then its extrapolation through a switch */
		res=g[k1];     /* so use linear extrapolation just this once 20/11/95*/
	else
#endif
	{ GHERMITE(res,x0,x1,y[k],y[k1],g[k],g[k1],t);}    /* 20/11/95*/
	history.lagmarker[i][markno]=k;
	return(res);
}


double pastvalue(int i, double t, int markno)
	/* Interogates the history ringbuffers. Note that there is a fair amount of
	   assignment of one variable to another at the start: this is to save
	   on address calculation and speed up the routine. 4/10/95*/

{ long k1,k,offset,offsetplus,size;
	double res,*y,*g,*x,x0,x1;
	y=history.buff[i];g=history.gbuff[i];
	x=history.clock; /*local pointers improve address efficiency*/
	offset=history.offset;size=history.size;
	if (x[offset]==t) return(y[offset]);
	offsetplus=offset+1L; if (offsetplus==size) offsetplus=0L;
	k=history.lagmarker[i][markno];
	k1=k+1L;if ((k1>=size)||(k1<0L)) k1=0L;
	while ((x[k1]<t)&&(k1!=offset))
	{ k1++;if (k1==size) k1=0L;}
	if (k1==0L) k=size-1L; else k=k1-1L;
	while ((x[k]>t)&&(k!=offsetplus))
	{ if (k==0L) k=size-1L; else k--;}
	k1=k+1L;if (k1==size) k1=0L;
	if (t<x[k])
	{
		Rprintf("\nERROR: lag for variable %d too large at %g\nx[k]=%g   k=%ld   t=%g\n",i,
				history.last_time-t,x[k],k,t);
		error("Lag too large for history buffer - try increasing the value of 'hbsize'");
	}
	x0=x[k];x1=x[k1];
#ifdef SWITCHES  /* some code for not extrapolating through a switch */
	if ((t>x[k1])&&(x[k]==x[k1])) /* then its extrapolation through a switch */
		res=y[k1]+(t-x[k1])*g[k1]; /* so use linear extrapolation just this once */
	else
#endif
	{ HERMITE(res,x0,x1,y[k],y[k1],g[k],g[k1],t);}
	/*  if (x[k1]==x[k]) return(y[k1]);*/
	/*printf("\nx0=%g   x1=%g  yk=%g  yk1=%g\n  gk=%g  gk1=%g  res=%g",
	  x0,x1,y[k],y[k1],g[k],g[k1],res); */
	history.lagmarker[i][markno]=k;
	return(res);
}


double zeropos(double x1, double x2, double x3, double s1, double s2, double s3)
	/* finds the root in [x1,x3] of a quadratic passing through the (xi,si)s
	   it is assumed that s3<s1*/


{ double z,y,zpy,a,b,c,d,a1,b1,c1,p;
	int ok=1;
	static int first=1;
	static double udge;
	if (first)
	{ first=0;
		udge=1.00000001;
	}
	z=x3-x2;y=x2-x1;zpy=z+y;
	if (z==0.0||y==0.0) error("Error in switching: zero switch interval");
	a1=a=s2;c1=c=(z*s1+y*s3-zpy*s2)/(zpy*z*y);b1=b=(s2-s1)/y+c*y;
	d=b*b-4.0*a*c;c*=2.0;
	p= -a/b; /* linear only approximation - in case c numerically zero */
	if (c==0.0) a=p;
	else
	{ if (d>=0.0)
		{ d=sqrt(d);a=(-b+d)/c;b=(-b-d)/c;
			if ((b>=-y)&&(b<=z)) a=b; else
				if ((a<-y)||(a>z)) ok=0;
		}
		if ((d<0.0)||(!ok))
		{ if (-s3<s1) a=z; else a=-y;}
		z=a1+a*b1+a*a*c1;
		d=a1+p*b1+p*p*c1;
		if (fabs(z)>fabs(d)) a=p; /* check that linear interpolation is not better */
	}
	a+=x2;
	if (a>x3) a=x3;
	if (a<=x1)
	{ if (a==0.0) a=udge-1.0; else if (a<0.0) a/=udge; else a*=udge;}
	return(a);
}


double istep(double *sw0, double *newsws, double *s0, double *news, double *g,
	     double *newg, double *c, double *err, double t0, double t1,
	     int nsw, int ns, int *flickedswitch, int clear)
	/* executes RK23 step to next switch or target depending on which comes first
	   If step is to the first switch then the number of that switch is returned
	   in flickedswitch, but map() is not called.
	   Returns how far it got. 5/10/95. g is assumed to contain the gradient at
	   t0, it will contain the gradient at time t1 on exit. This improves
	   efficiency by making use of info. calculated in the previous step of rk23.
	   */

{ 
	static int first=1,nsold,nswold,*flicked;
	static double *err1,*s1,*s2,*sw1,*sw2;
	int k,i,switches=0;
	double zp,dt,sp2,sp1,minp,udge,ds;
	if(clear && first == 0) { /* Bobby */
		free(sw1);free(s1);free(sw2);free(s2); free(flicked); free(err1);
		first = 1;
		return 0;
	} else if (clear) return 0;
	if ((first)||(ns!=nsold)||(nsw!=nswold)) {
		if (!first) {
			free(sw1);free(s1);free(sw2);free(s2);
			/* Bobby */
			free(flicked); free(err1);
		}
		first=0;
		sw1=(double *)calloc(nsw,sizeof(double));
		sw2=(double *)calloc(nsw,sizeof(double));
		s1=(double *)calloc(ns,sizeof(double));
		s2=(double *)calloc(ns,sizeof(double));
		err1=(double *)calloc(ns,sizeof(double));
		flicked=(int *)calloc(nsw,sizeof(int));
		nsold=ns;nswold=nsw;
	}
	dt=t1-t0;
	rk23(s0,news,g,newg,err,c,ns,t0,dt,0);
	if (nsw) 
		switchfunctions(newsws,news,c,t1);
	for (i=0;i<nsw;i++) {                 /* are there any switches */
		if ((sw0[i]>0.0)&&(newsws[i]<=0.0)) {
			flicked[switches]=i;
			switches++;
		}
	}

	if (!switches) {  /* No switches so its an ordinary step */
		*flickedswitch=-1;
		return(t1);
	}

	/* Logic for stepping to first switch */
	sp1=t0+dt*0.5;
	for (k=0;k<100;k++) /* if k gets to 100 routine fails */
	{ 
	  rk23(s0,s1,g,newg,err,c,ns,t0,sp1-t0,0); /* step to approx. 1st switch position */
		switchfunctions(sw1,s1,c,sp1);

		switches=0;
		for (i=0;i<nsw;i++)     /* are there any switches ? MACRO after debug*/
			if ((sw0[i]>0.0)&&(sw1[i]<=0.0))
			{ flicked[switches]=i;switches++;}

		if ((k)&&(switches==1))  /* MACRO after debug */
		{ *flickedswitch=flicked[0];
			for (i=0;i<ns;i++) news[i]=s1[i];
			for (i=0;i<nsw;i++) newsws[i]=sw1[i];
			return(sp1);
		}

		rk23(s1,s2,newg,newg,err1,c,ns,sp1,t1-sp1,0);/* step to end of interval */
		switchfunctions(sw2,s2,c,t1);

		for (i=0;i<nsw;i++)     /* are there any switches ? MACRO after debug*/
			if ((sw1[i]>0.0)&&(sw2[i]<=0.0))
			{ 
				flicked[switches]=i;
				switches++;
			}

		if (!switches)  /* MACRO after debug */
		{ 
			*flickedswitch=-1;
			for (i=0;i<ns;i++)
			{ news[i]=s2[i]; err[i]=sqrt(err[i]*err[i]+err1[i]*err1[i]);}
			for (i=0;i<nsw;i++) newsws[i]=sw2[i];
			return(t1);
		}

		/* having got this far switch positions must be estimated */

		/* locate the first switch */
		sp2=t1;minp=t1;
		for (i=0;i<switches;i++)
		{ 
			if ((t0==t1)||(sp1==t1)||(t0==sp1)) 
				zp=t1; 
			else
				zp=zeropos(t0,sp1,t1,sw0[flicked[i]],sw1[flicked[i]],sw2[flicked[i]]);
			if (zp<minp) { 
				sp2=minp;minp=zp;
			}
		}
		sp1=minp;udge=1e-9;ds=sp2-sp1;
		if (ds>0.0)
			do { 
				sp1+=udge*ds;udge*=10.0;
			} while (sp1==minp);
	}
	error("Problem with switch logic");return(t0);

}

/* void dde(...)
     s          State variables
     c          Coefficients
     t0, t1     Start and stop times
     dt         Pointer to initial timestep (returns final step - which
                is step that would have been used if t1 not reached!)
                max step set to 100 times this, minimum set to 1e-9 of this
     eps        Fractional tolerance for adaptive stepping
     dout       Interval for output via user routine output(). Every
                time a step passes 1 or more times of the form t0+i*dout
                output() is called once. Hence output is only roughly
                regular. dout=0.0 for no output.
     otimes     An ordered array of times >= t0 at which
                the state variables should be recorded for output
     no_otimes  Length of otimes
     hbsize     The number of elements to store in the history buffer
     nsw        Numbwer of switches
     ns         Number of state variables
     nhv        Number of lagged variables
     reset      Set to 0 not to reset, to 1 to reset
     nlag       Number of place markers per history variable
     fixstep    Set to 0 for adaptive stepping, or 1 for fixed timestep
*/
void dde(double *s, double *c, double t0, double t1, double *dt, double eps,
	 double *otimes, int no_otimes, int ns, int nsw, int nhv, long hbsize,
	 int nlag, int reset, int fixstep, int clear)
{
  double D,Da,errmax,rerr,target,t=0;
  static double *err,*newsws,*sws,*news,*newg,*dum,*sp,*nswp,*swp,*nsp,*e0,*scale;
	static double *g,mindt,maxdt;
	static double oldt, *sout, *olds, *oldg; /* bjc 2007-05-08*/
	static int first=1;
	long i,iout;
	int swi=0;

	if(clear && first == 0) { /* Bobby */
		free(g); first = 1; 
	} 

	/* clear if allocated when interrupted */
	/* Bobby: not sure how it is getting interrupted but it is */
	if(clear) { /* Bobby */
		if(swp) { free(swp); swp = NULL; }
		if(nswp) { free(nswp); nswp = NULL; }
		if(nsp) { free(nsp); nsp = NULL; }
		if(err) { free(err); err = NULL; }
		if(e0) { free(e0); e0 = NULL; }
		if(newg) { free(newg); newg = NULL; }
		if(scale) { free(scale); scale = NULL; }
		if(olds) { free(olds); olds = NULL; }
		if(oldg) { free(oldg); oldg = NULL; }
		if(sout) { free(sout); sout = NULL; }
		return;
	}

	swp = nswp = nsp = err = e0 = newg = scale = NULL; /* Bobby */
	sout = olds = oldg = NULL; /* Bobby */

	nswp=newsws=(double *)calloc(nsw,sizeof(double));
	swp=sws=(double *)calloc(nsw,sizeof(double));
	nsp=news=(double *)calloc(ns,sizeof(double));
	newg=(double *)calloc(ns,sizeof(double));
	err=(double *)calloc(ns,sizeof(double));
	e0=(double *)calloc(ns,sizeof(double));
	scale=(double *)calloc((size_t)ns,sizeof(double));
	olds=(double *)calloc(ns, sizeof(double)); /* bjc 2007-05-08*/
	oldg=(double *)calloc(ns, sizeof(double)); /* bjc 2007-05-08*/
	statescale(scale);
	if (nsw) switchfunctions(sws,s,c,t0);
	if (reset) { 
		if (first) 
			first=0;
		else 
			free(g);
		mindt=(*dt)*1e-9;
		maxdt=(*dt)*100.0;
		g=(double *)calloc(ns,sizeof(double));
		sout=(double *)calloc(ns, sizeof(double)); /* bjc 2007-05-08*/
		inithisbuff(nhv,hbsize,nlag,0);
		grad(g,s,c,t0);
		updatehistory(g,s,c,t0,0);
	}

	t=t0;
	sp=s;
	D = (*dt);

	iout = 0L;
	if (iout < no_otimes && otimes[iout] == t0) {
		output(s, otimes[iout++]); 
	} /* bjc 2007-05-08*/

	while (t0<t1) {
		if (t0+D>t1) { 
			target=t1;
		} else { 
			target=t+D;
		}
		/* bjc 2007-05-08: copy current to old */
		for( i=0; i<ns; i++ ) {
			oldg[i]=g[i];
			olds[i]=s[i];
		}
		oldt=t; 
		t=istep(sws,newsws,s,news,g,newg,c,err,t0,target,nsw,ns,&swi,0);
		errmax=0.0;
		if ((!fixstep)&&(D>mindt)) { /* error control: see Press et al. 1992 p718*/
			for (i=0;i<ns;i++) /*e0[i]=eps*(fabs(s[i])+fabs((t-t0)*g[i])+1e-20);*/
				e0[i]=eps*(fabs(s[i])+fabs((t-t0)*(g[i]+newg[i])*0.5)+scale[i]);
			for (i=0;i<ns;i++) {
				if (err[i]<1e-150&&err[i]>-1e-150) rerr=0.0;else
					rerr=err[i]/e0[i];
				rerr=fabs(rerr);
				if (rerr>errmax) 
					errmax=rerr;
			}
		}
		if (errmax<1.0) { 
			accepted++;
			/*printf("\nt = %g    dt = %g",t,*dt);*/
			dum=s;s=news;news=dum;dum=sws;sws=newsws;newsws=dum;
			dum=g;g=newg;newg=dum;
			updatehistory(g,s,c,t,0);
			/* bjc 2007-05-08: new fancy outputting */
			/* Interpolate using HERMITE between current and last accepted points */
			while( iout < no_otimes && otimes[iout] <= t ) { /* outputting results */
				for (i=0;i<ns;i++) {
					HERMITE(sout[i], oldt, t, olds[i], s[i], oldg[i], g[i], otimes[iout]);
				}
				output(sout,otimes[iout++]);
			}

			t0=t;
			if ( swi > -1 ) { 
				output(s,t);
				map(s,c,t,swi);
				output(s,t);
				/* printf("\nStep at t=%g\n",t); */
				grad(g,s,c,t);updatehistory(g,s,c,t,0);
			} else {
				/* increase stepsize */
				if ((!fixstep)&&(t<t1))
					D = (errmax > ERRCON ? 0.95*D*pow(errmax,-0.2) : 5.0*D);
				if (D>maxdt) 
					D=maxdt;
			}
		} else { 
			Da=t-t0; /* Step actually achieved */
			rejected++;
			/* Shrink D from Da */
#ifndef STEPOFF
			D=0.95*Da*pow(errmax,-0.25);
			D=( D < 0.1*Da ? 0.1*Da : D );
			if (D < 1e-14*(*dt)) {/* printf("\nStepsize failure\n");*/}
#endif
			t=t0;
		}
	}
	(*dt)=D;
	/* bjc 2007-05-08: the final output is already covered above */
	/* if (dout && (otimes[no_otimes-1]==t1)) output(s,t1); */
	for (i=0;i<ns;i++) sp[i]=s[i]; /* copying results to correct address */
	free(swp);free(nswp);free(nsp);free(err);free(e0);free(newg);free(scale);
	swp = nswp = nsp = err = e0 = newg = scale = NULL; /* Bobby */
	free(sout); free(olds); free(oldg); /* bjc 2007-05-08*/
	sout = olds = oldg = NULL; /* Bobby */
}


