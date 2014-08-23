/*-----------------------------------------------------------------------------------------
 * 2010-01-03 hxsim_ppoe1.c: my MEX-file for MATLAB, to replace
 * the simulation loops part in phxsim_sim.m. the random number
 * generator in this program("rans.h & rans.c" is from http://www.cs.wm.edu/~va/software/park/
 * ------------------------------------------------------------------------------------------*/

#include "mex.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "rngs.h"

#define MODULUS    2147483647 /* DON'T CHANGE THIS VALUE                  */
#define MULTIPLIER 48271      /* DON'T CHANGE THIS VALUE                  */
#define CHECK      399268537  /* DON'T CHANGE THIS VALUE                  */
#define STREAMS    256        /* # of streams, DON'T CHANGE THIS VALUE    */
#define A256       22925      /* jump multiplier, DON'T CHANGE THIS VALUE */
#define DEFAULT    123456789  /* initial seed, use 0 < DEFAULT < MODULUS  */

static long seed[STREAMS] = {DEFAULT};  /* current state of each stream   */
static int  stream        = 0;          /* stream index, 0 is the default */
static int  initialized   = 0;          /* test for stream initialization */


double Random(void)
/* ----------------------------------------------------------------
 * Random returns a pseudo-random real number uniformly distributed
 * between 0.0 and 1.0.
 * ----------------------------------------------------------------
 */
{
    const long Q = MODULUS / MULTIPLIER;
    const long R = MODULUS % MULTIPLIER;
    long t;
    
    t = MULTIPLIER * (seed[stream] % Q) - R * (seed[stream] / Q);
    if (t > 0)
        seed[stream] = t;
    else
        seed[stream] = t + MODULUS;
    return ((double) seed[stream] / MODULUS);
}


void PlantSeeds(long x)
/* ---------------------------------------------------------------------
 * Use this function to set the state of all the random number generator
 * streams by "planting" a sequence of states (seeds), one per stream,
 * with all states dictated by the state of the default stream.
 * The sequence of planted states is separated one from the next by
 * 8,367,782 calls to Random().
 * ---------------------------------------------------------------------
 */
{
    const long Q = MODULUS / A256;
    const long R = MODULUS % A256;
    int  j;
    int  s;
    
    initialized = 1;
    s = stream;                            /* remember the current stream */
    SelectStream(0);                       /* change to stream 0          */
    PutSeed(x);                            /* set seed[0]                 */
    stream = s;                            /* reset the current stream    */
    for (j = 1; j < STREAMS; j++) {
        x = A256 * (seed[j - 1] % Q) - R * (seed[j - 1] / Q);
        if (x > 0)
            seed[j] = x;
        else
            seed[j] = x + MODULUS;
    }
}


void PutSeed(long x)
/* ---------------------------------------------------------------
 * Use this function to set the state of the current random number
 * generator stream according to the following conventions:
 *    if x > 0 then x is the state (unless too large)
 *    if x < 0 then the state is obtained from the system clock
 *    if x = 0 then the state is to be supplied interactively
 * ---------------------------------------------------------------
 */
{
    char ok = 0;
    
    if (x > 0)
        x = x % MODULUS;                       /* correct if x is too large  */
    if (x < 0)
        x = ((unsigned long) time((time_t *) NULL)) % MODULUS;
    if (x == 0)
        while (!ok) {
        printf("\nEnter a positive integer seed (9 digits or less) >> ");
        scanf("%ld", &x);
        ok = (0 < x) && (x < MODULUS);
        if (!ok)
            printf("\nInput out of range ... try again\n");
        }
    seed[stream] = x;
}


void GetSeed(long *x)
/* ---------------------------------------------------------------
 * Use this function to get the state of the current random number
 * generator stream.
 * ---------------------------------------------------------------
 */
{
    *x = seed[stream];
}


void SelectStream(int index)
/* ------------------------------------------------------------------
 * Use this function to set the current random number generator
 * stream -- that stream from which the next random number will come.
 * ------------------------------------------------------------------
 */
{
    stream = ((unsigned int) index) % STREAMS;
    if ((initialized == 0) && (stream != 0))   /* protect against        */
        PlantSeeds(DEFAULT);                     /* un-initialized streams */
}


void TestRandom(void)
/* ------------------------------------------------------------------
 * Use this (optional) function to test for a correct implementation.
 * ------------------------------------------------------------------
 */
{
    long   i;
    long   x;
    double u;
    char   ok = 0;
    
    SelectStream(0);                  /* select the default stream */
    PutSeed(1);                       /* and set the state to 1    */
    for(i = 0; i < 10000; i++)
        u = Random();
    GetSeed(&x);                      /* get the new state value   */
    ok = (x == CHECK);                /* and check for correctness */
    
    SelectStream(1);                  /* select stream 1                 */
    PlantSeeds(1);                    /* set the state of all streams    */
    GetSeed(&x);                      /* get the state of stream 1       */
    ok = ok && (x == A256);           /* x should be the jump multiplier */
    if (ok)
        printf("\n The implementation of rngs.c is correct.\n\n");
    else
        printf("\n\a ERROR -- the implementation of rngs.c is not correct.\n\n");
}


/*********************************************************************************************************/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
        const mxArray *prhs[]) {
    double deltaT, START, END, Dfraction, flagSim;
    double Hfraction;
    double *pf_U, *pf_I, *pf_Ix, *pf_N, *kcHD, *kcDH, *k_fm, *FoldStates, *HDmatrix;
    double *pf, ran;
    int simSteps;
    int step, M, N, i, j;
    
    /* Check for proper number of arguments. */
    if (nrhs != 15)
        mexErrMsgTxt("15 input arguments required.");
    else if (nlhs > 0)
        mexErrMsgTxt("Too many output arguments.");
    
    simSteps = mxGetScalar(prhs[0]);
    deltaT = mxGetScalar(prhs[1]); 
    START = mxGetScalar(prhs[2]);
    END = mxGetScalar(prhs[3]);
    Dfraction = mxGetScalar(prhs[4]);
    Hfraction = 1 - Dfraction;
    pf_U = mxGetPr(prhs[5]);
    pf_I = mxGetPr(prhs[6]);
    pf_Ix = mxGetPr(prhs[7]);
    pf_N = mxGetPr(prhs[8]);
    kcHD = mxGetPr(prhs[9]);
    kcDH = mxGetPr(prhs[10]);
    k_fm = mxGetPr(prhs[11]);  /* [kUI, kIU, kIIx, kIxI, kIN, kNI] */
    FoldStates = mxGetPr(prhs[12]);
    HDmatrix = mxGetPr(prhs[13]);
    flagSim = mxGetScalar(prhs[14]);
    
    M = mxGetM(prhs[13]);
    N = mxGetN(prhs[13]);
    
    mexPrintf("current simSteps = %d\n running...\n", simSteps);
    for (step=1; step<=simSteps; step++) {
        
        /*SelectStream(step%256); */     
        
        for (i=0; i<M; i++) {
            switch ((int)FoldStates[i]) {
                case 0:
                    pf=pf_U;
                case 1:
                    pf=pf_I;
                case 2:
                    pf=pf_Ix;
                case 3:
                    pf=pf_N;
            }
            
            for (j=(int)START; j<=(int)END; j++) {
                ran=Random();
                switch ((int)*(HDmatrix+i+M*(j-(int)START))) {    /*All MATLAB data is stored columnwise*/
                    case 0: /* H */
                        if (ran <= (1-exp(-kcHD[j-1]*deltaT/pf[j-1]))*Dfraction)
                            *(HDmatrix+i+M*(j-(int)START))=1; /* H->D */
                    case 1: /* D */
                        if (ran <= (1-exp(-kcDH[j-1]*deltaT/pf[j-1]))*Hfraction)
                            *(HDmatrix+i+M*(j-(int)START))=0; /* D->H */
                }
            }
            
            if (flagSim!=0) {
                ran=Random();
                switch ((int)FoldStates[i]) {                       
                    case 0: /* U */
                        if (ran <= 1-exp(-k_fm[0]*deltaT))
                            FoldStates[i]=1; /* U->I */
                    case 1: /* I */ {
                        if (ran <= 1-exp(-k_fm[1]*deltaT))
                            FoldStates[i]=0; /* I->U */
                        if (ran > 1-exp(-k_fm[1]*deltaT) && ran <= (1-exp(-k_fm[1]*deltaT))+(1-exp(-k_fm[4]*deltaT)))
                            FoldStates[i]=3; /* I->N */
                        if (ran>(1-exp(-k_fm[1]*deltaT))+(1-exp(-k_fm[4]*deltaT)) && ran<=(1-exp(-k_fm[1]*deltaT))+(1-exp(-k_fm[4]*deltaT))+(1-exp(-k_fm[2]*deltaT)))
                            FoldStates[i]=2; /* I->Ix */}
                    case 2: /* Ix */
                        if (ran <= 1-exp(-k_fm[3]*deltaT))
                            FoldStates[i]=1; /* Ix->I */
                    case 3: /* N */
                        if (ran <= 1-exp(-k_fm[5]*deltaT))
                            FoldStates[i]=1; /* N->I */
                }
            }
        }
    }
}
