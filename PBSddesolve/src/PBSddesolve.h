/*=============================================================================
  Copyright (C) 2007-2019 Fisheries and Oceans Canada

  This file is part of PBSddesolve.

  PBSddesolve is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  PBSddesolve is distributed WITHOUT ANY WARRANTY;
  without even the implied warranty of MERCHANTABILITY or 
  FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  Free Software Foundation, Inc.
  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
=============================================================================*/

#ifndef _PBSDDESOLVE_H_
#define _PBSDDESOLVE_H_

#include <R.h>
#include <Rdefines.h>

/*-----------------------------------------------------------------------------
  startDDE:
    Simon Wood’s (1999) numerical routines produce the core functionality of
    PBSddesolve.

  Arguments:
    gradfunc   = gradient function;
    switchfunc = optional function that determines conditions when the DDE
                 system experiences switches;
    mapfunc    = optional function associated with switchfunc that describes
                 how DDE system changes (poss.discontinuously) at switch times;
    env        = environment (probably new);
    y          = vector of initial values for the states (also determines n);
    parms      = optional vector of parameters to pass to gradfunc;
    settings   = vector of:
                 tol  = scalar that sets the max. error tolerated in solution;
                 from = first times value;
                 to   = maximum times value;
                 dt   = maximum initial time step used in constructing the
                        numerical solution;
                 hbsize = history buffer size required for retaining lagged
                          state variable values;
    times      = numeric vector of explicit times at which the solution
                 should be obtained.
  ---------------------------------------------------------------------------*/
SEXP startDDE(SEXP gradFunc, SEXP switchFunc, SEXP mapFunc, SEXP env,
              SEXP yinit, SEXP parms, SEXP settings, SEXP outtimes);


/*-----------------------------------------------------------------------------
  getPastValue:
    These routines provides access to variable history at lagged times.

  Arguments:
    t      = specific time in history;
    markno = used for optimization when more than one lag time is used.
  ---------------------------------------------------------------------------*/
SEXP getPastValue(SEXP t, SEXP markno);


/*-----------------------------------------------------------------------------
  getPastGradient:
    These routines provides access to variable history at lagged times.

  Arguments:
    t      = specific time in history;
    markno = used for optimization when more than one lag time is used.
  ---------------------------------------------------------------------------*/
SEXP getPastGradient(SEXP t, SEXP markno);


#endif /* _PBSDDESOLVE_H_ */
