#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "PBSddesolve.h"
#include "ddeq.h"

static const R_CMethodDef CEntries[]  = {
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"startDDE",         (DL_FUNC) &startDDE, 8},         /* PBSddesolve.c: 309 */
    {"getPastValue",     (DL_FUNC) &getPastValue, 2},     /* PBSddesolve.c: 442 */
    {"getPastGradient",  (DL_FUNC) &getPastGradient, 2},  /* PBSddesolve.c: 474 */
    {NULL, NULL, 0}
};

void 
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
R_init_PBSddesolve(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

