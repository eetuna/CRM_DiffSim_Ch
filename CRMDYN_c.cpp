#include <iostream>
#include <cmath>
#include <chrono> 
#include "CRMDYN.hpp"
using namespace std::chrono;

#define EPS 1e-5
#define NUM_DYN_STATE  (NUM_ACT_SET*6 + NUM_FLEX_SEG*3 + NUM_ACT_SET*6 + NUM_ACT_SET* (3+9) )    // (v, w, u, n, m, p, R)

/*   Copyright 2005-2015 The MathWorks, Inc. */
/*   Written by Peter Lindskog. */

/* Include libraries. */
#include "/usr/local/MATLAB/R2022b/extern/include/mex.h"

/* Specify the number of outputs here. */
#define NY 3


template <typename T>
void printMatrix(T *p, int D1, int D2, const char *text) {
    std::cout << text << " --" << std::endl;
    for (int i=0; i<D1; i++) {
        for (int j=0; j<D2; j++) {
            std::cout << *(p+i*D2+j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "----" << std::endl;
}

/* State equations. */
void compute_dx(double *dx, double t, double *x, double *u, double **p,
                const mxArray *auxvar)
{

    /** Retrieve x. **/
    double v_L_pre[NUM_ACT_SET][3], w_L_pre[NUM_ACT_SET][3],  u0_initialguess[NUM_FLEX_SEG][3], nL_initialguess[NUM_ACT_SET][3], mL_initialguess[NUM_ACT_SET][3], pL_pre[NUM_ACT_SET][3], RL_pre[NUM_ACT_SET][9];
    // Define initial guesses to be used when solving boundary value problem
    for (int i = 0; i < NUM_ACT_SET; ++i) {
        for (int j = 0; j < 3; ++j) {
            v_L_pre[i][j] = x[j+i*3];
        }
    }

    for (int i = 0; i < NUM_ACT_SET; ++i) {
        for (int j = 0; j < 3; ++j) {
            w_L_pre[i][j] = x[NUM_ACT_SET * 3 + j+i*3];
        }
    }

    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
        for (int j = 0; j < 3; ++j) {
            u0_initialguess[i][j] = x[NUM_ACT_SET*6 + j + 3*i];

        }
    }


    for (int i = 0; i < NUM_ACT_SET; ++i) {
        for (int j = 0; j < 3; ++j) {
            mL_initialguess[i][j] = x[NUM_ACT_SET*6 + NUM_FLEX_SEG*3 + j+ 3*i];
        }
    }

    for (int i = 0; i < NUM_ACT_SET; ++i) {
        for (int j = 0; j < 3; ++j) {
            nL_initialguess[i][j] = x[NUM_ACT_SET*9 + NUM_FLEX_SEG*3 + j+ 3*i];
        }
    }
    for (int i = 0; i < NUM_ACT_SET; ++i) {
        for (int j = 0; j < 3; ++j) {
            pL_pre[i][j] = x[NUM_ACT_SET*12 + NUM_FLEX_SEG*3 + j+i*3];
        }
    }

//    int ind_r = NUM_ACT_SET*15 + 3;
    for (int i = 0; i < NUM_ACT_SET; ++i) {
        for (int j = 0; j < 9; ++j) {
            RL_pre[i][j] = x[NUM_ACT_SET*15 + NUM_FLEX_SEG*3 + j+i*9];
        }
    }

//    std::cout <<" v_L_pre : " << v_L_pre[0] << " " << v_L_pre[1] << " " << v_L_pre[2] << std::endl;
//    std::cout <<" w_L_pre: " << w_L_pre[0] << " " << w_L_pre[1] << " " << w_L_pre[2] << std::endl;
//
////    std::cout << "u0_initialguess: " << u0_initialguess[0] << " " << u0_initialguess[1] << " " << u0_initialguess[2] <<  std::endl;
//    std::cout << "pL_pre: " << pL_pre[0] << " " << pL_pre[1] << " " << pL_pre[2] <<  std::endl;
//
//
    /** Retrieve u. **/
    double ActuationCurrents[NUM_ACT_SET][3];
    for (int i = 0; i < NUM_ACT_SET; i++)	for (int j = 0; j < 3; j++)	ActuationCurrents[i][j] = u[i * 3 + j];
//    std::cout << "ActuationCurrents: " << ActuationCurrents[0][0] << " " << ActuationCurrents[0][1] << " " << ActuationCurrents[0][2] <<  std::endl;

    /** Retrieve model parameters. **/
    double damping_[NUM_ACT_SET][6];
//    for (int i = 0; i < 6; ++i) {
//        damping_[i] = p[0][i];
//    }
    for (int i = 0; i < NUM_ACT_SET; ++i) {
        damping_[i][0] = damping_[i][1] = p[0][0 +i*4];
        damping_[i][2] = p[0][1+i*4];
        damping_[i][3] = damping_[i][4] = p[0][2+i*4];
        damping_[i][5] = p[0][3+i*4];
    }

//    std::cout <<"damping_ v : " << damping_[0] << " " << damping_[1] << " " << damping_[2] << std::endl;
//    std::cout <<"damping_ w: " << damping_[3] << " " << damping_[4] << " " << damping_[5] << std::endl;
    double Delta_T = p[1][0];

    // Declare the output variables for BVP
    // calculated curvature at the catheter base
    double u0_calc[NUM_FLEX_SEG][3];
    double nL_calc[NUM_ACT_SET][3];
    double mL_calc[NUM_ACT_SET][3];
    // calculated contstraint force at the catheter tip (this will be used when ContactMode == ContactModeType::FIXED_TIP)
    double ftip_calc[3];
    // numerical nonlinear equation solver diagnostic outputs
    int localmin;

//    std::cout <<" damping_: " << damping_[0] << " " << damping_[1] << " " << damping_[2] << " " << damping_[3] << " " << damping_[4] << " " << damping_[5] << std::endl;
//    std::cout <<" Delta_T: " << Delta_T << std::endl;

    double oRlist[NUM_FLEX_SEG] = { p[2][0], p[2][0]}; //{ 1.5875, 1.5875 };
    // Inner radii of each of the flexible segments - unit: mm
    double iRlist[NUM_FLEX_SEG] = { p[2][1], p[2][1]}; // { 0.9906, 0.9906 };
    // Young's moduli of each of the flexible segments - unit: ??
    double YoungModlist[NUM_FLEX_SEG] = {p[3][0], p[3][0]}; // { 4.9105, 4.9105 };
    // Shear moduli of each of the flexbile segments - unit: ??
    double ShearModlist[NUM_FLEX_SEG] = {p[3][1],p[3][1] };// { 1.6545, 1.6545 };
    // Alignment angles for the side coils, for each of the actuator sets - unit: radians
    double CoilAlignmentAngles[NUM_ACT_SET][2] ={ {p[4][0], p[4][1]}  };// { {-0.0871, -0.3934}  };
    // Coil turn area matrices for each of the actuator sets  - unit: mm2
    double CoilTurnAreaMat[NUM_ACT_SET][9] =  { { p[5][0], 0.0, 0.0, 0.0, p[5][1], 0.0, 0.0, 0.0, p[5][2] } };// { { 1.44, 0.0, 0.0, 0.0, 1.3851, 0.0, 0.0, 0.0, 1.60 } };
    // Lengths of each of the catheter segments - unit: mm
    double SegmentLengths[NUM_SEGMENTS] = { 19.85,  18.3, 59.40 };

    // Mass of each of the actuator sets - unit: ??
    double ActMass[NUM_ACT_SET] ={p[6][0]};// { 7.7736e-5};
//    std::cout <<" ActMass: " << ActMass[0] << std::endl;

    // Array of lambda values for marker locations  (distance from the tip to each of the markers) - unit: mm
    double MarkerLoc[NUM_LOCALIZATION_MARKERS] = { 2.18, 18.79, 50.06, 101.52, 104.02 };
    // Length density for each of the catheter segments
    double rho[NUM_SEGMENTS] = { 2.8814e-7, 2.8814e-7, 2.8814e-7 };
    // Local curvature in unloaded configuration for each of the flexible segments; NUM_FLEX_SEG*3 long array corresponding to NUM_FLEX_SEG many 3x1 vectors
    double ustarlist[NUM_FLEX_SEG][3] = { 0.0148626102272827, 0.0094448815853795, 0, 0.000799849479553773, -0.0001892201697658481, 0};//0.00567, -0.00822, 0.0};


    // *** Catheter Configuration in spatial coordinates
    // B0 field vector of the MRI scanner (in spatial coordinates) - unit: Tesla
    double B0[3] = { 0.0, 3.0, 0.0 };
    // Gravity vector - unit: ??
    double gravity[3] = { 0.0, 0.0, 9.81 };
    // Catheter entry point coordinates (in spatial coordinates) - unit: mm
    double p0[3] = { 0.0, 0.0, 0.0 };
    // Catheter entry point orientation (3x3 rotation matrix describing catheter entry point frame orientation relative to the spatial frame stored in row major order) - unit: unitless
    double R0[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

    // *** Other External variables
    // specify if catheter is in free space or if the catheter tip is constrained to a contact point
    ContactModeType ContactMode = ContactModeType::FREE_TIP;
    // External point force (in spatial coordinates) applied at the tip of the catheter (\lambda = 0)  - unit: ??
    //   (this will be used when ContactMode == ContactModeType::FREE_TIP)
    double TipForce[3] = { 0.0, 0.0, 0.0 };
    // The spatial coordinates of the point where the catheter tip is constrained to be
    //   (this will be used when ContactMode == ContactModeType::FIXED_TIP)
    double TipConstraintPoint[3] = { 0.0, 0.0, 0.0 };


    // *** Numerical Computation Parameters
    // Stepsize used in numerical integration along the length of the catheter during IVP - unit: mm
    double IntegrationStepSize = 0.2;

    /**
     * These are hard coded, need to revise these later
     */
    // we are adding the tubing mass of the coil section to the total mass of actuator:
//        double tubing_mass = rho[0] * SegmentLengths[1];
//    ActMass[0] += tubing_mass;
    double ActInertia[NUM_ACT_SET][9];
    for (int i = 0; i < NUM_ACT_SET; ++i)
    {
        double I_zz = 0.5 * (ActMass[i]) * (oRlist[0] * oRlist[0] + iRlist[0] * iRlist[0]);
        double I_xx = 0.25 * (ActMass[i]) * (oRlist[0] * oRlist[0] + iRlist[0] * iRlist[0]) + 1.0 / 12 * (ActMass[i]) * SegmentLengths[2*i+1] * SegmentLengths[2*i+1];
        ActInertia[i][0] = I_xx; ActInertia[i][1] = 0.0; ActInertia[i][2] = 0.0;
        ActInertia[i][3] = 0.0; ActInertia[i][4] = I_xx; ActInertia[i][5] = 0.0;
        ActInertia[i][6] = 0.0; ActInertia[i][7] = 0.0; ActInertia[i][8] = I_zz;
    }


    // Define variables of convenience that will be use to package catheter physical parameters and spatial configuration parameters
    CRMCatheterModelParams CathParams;		// Catheter physical parameters will be packaged into this structure
    CatheterConfiguration CathConfig;		// Catheter spatial configuration parameters will be packages into this structure


    CRMShootingMethodBVP_Prep(B0,gravity, p0, R0,SegmentLengths, MarkerLoc,
                              iRlist, oRlist, YoungModlist, ShearModlist, ustarlist,
                              CoilAlignmentAngles, CoilTurnAreaMat, rho, ActMass, ActInertia,
                              CathParams, CathConfig);


    double ftip_initialguess[3] = { 0.0, 0.0, 0.0 };

    /**
      * initial curvature is assumed constant curvature, we are initializing the full length u_history,
      * Prepare for interpolation of dynamic insertion length
      */
    double SegEnds[NUM_SEGMENTS+1];
    SegEnds[0] = 0.0;
    for (int i = 1; i < NUM_SEGMENTS+1; ++i) {
        SegEnds[i] = SegEnds[i-1] + SegmentLengths[NUM_SEGMENTS - i]; // SegmentLengths is ordered from the distal
    }
//    std::cout << "segends: " << SegEnds[0] << " " << SegEnds[1] << " " << SegEnds[2] << " " << SegEnds[3] << std::endl;

    int SegSteps[NUM_FLEX_SEG];
    double h0[NUM_FLEX_SEG];
    for (int i=0; i<NUM_FLEX_SEG; i++) { 		// flexible catheter segment
        // Calculate the number of integration steps based on the given IntegrationStepSize
        SegSteps[i]=int(ceil( (SegEnds[2*i+1]-SegEnds[2*i]) / IntegrationStepSize ) );
        h0[i] = dVal((SegEnds[2*i+1]-SegEnds[2*i]))/(SegSteps[i]*1.0);
    }


    /** Compute next state **/
    CRMShootingMethodParams<double> BVPParams{};

    // Declare output variables
    // catheter shape state at the entry point of the catheter
    //   states are packed u[0..2], R[0..8], p[0..2](R: 3x3 matrix stored in row major order R11 R12 R13 R21 R22 R23 R31 R32 R33)
    double xf[NUM_STATES];
    // moment residual at the catheter tip - this should converge to {0,0,0} if the catheter is at its equilibrium configuration
//    std::cout << "ActuationCurrents: " << ActuationCurrents[0][0] << " " << ActuationCurrents[0][1] << " " << ActuationCurrents[0][2] <<  std::endl;
//
//
//    std::cout << "pL_pre: " << pL_pre[0] << " " << pL_pre[1] << " " << pL_pre[2] <<  std::endl;
//    std::cout << "RL: " << std::endl;
//    std::cout <<  RL_pre[0] << " " << RL_pre[1] << " " << RL_pre[2] <<  std::endl;
//    std::cout <<  RL_pre[3] << " " << RL_pre[4] << " " << RL_pre[5] <<  std::endl;
//    std::cout <<  RL_pre[6] << " " << RL_pre[7] << " " << RL_pre[8] <<  std::endl;

    double InsertedLength =  0.0;//98.5; // u[NUM_CONTROL -1 ];
    for (int i = 0; i < NUM_SEGMENTS; ++i) {
        InsertedLength += SegmentLengths[i]; // we are not controlling this now
    }


//        ActuationCurrents[0][2] = 0.2;
    CRMConstructShootingMethodParamSet<double>(CathParams, CathConfig, InsertedLength, ActuationCurrents, ContactMode,
                                               TipConstraintPoint, TipForce, IntegrationStepSize,
                                               v_L_pre, w_L_pre, pL_pre,  RL_pre, damping_, Delta_T,
                                               BVPParams);


    DynamicsBVP(BVPParams, u0_initialguess, mL_initialguess, nL_initialguess, ftip_initialguess,
                u0_calc, mL_calc, nL_calc, ftip_calc, localmin);

    double out_ReportedMarkerPos[NUM_LOCALIZATION_MARKERS][3], x_coil[NUM_COIL_STATES];

    DYNSolverIVP(BVPParams, u0_calc, mL_calc, nL_calc, ftip_calc,
                 true, xf, x_coil,out_ReportedMarkerPos);

    /** Output report **/
//    for (int i = 0; i < NUM_DYN_STATE; ++i) {
//        if(i < 6){   dx[i] = x_coil[i]; } //v, w
//        else if (i < 9) dx[i] = u0_calc[i - 6];
//        else if (i < 12) dx[i] = mL_calc[0][i - 9];
//        else if (i < 15) dx[i] = nL_calc[0][i - 12];
//        else if (i < 18) dx[i] = x_coil[i - 9]; //p
//        else dx[i] = x_coil[i - 9];//R
//    }

    for (int i = 0; i < 6; ++i) {
        dx[i] = x_coil[i];
    }
    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
        for (int j = 0; j < 3; ++j) {
            dx[6 + j + 3*i] = u0_calc[i][j];
        }
    }
    for (int i = 0; i < 3; ++i) {
        dx[6 + NUM_FLEX_SEG*3 +i ]= mL_calc[0][i];
    }
    for (int i = 0; i < 3; ++i) {
        dx[6 + NUM_FLEX_SEG*3 + 3 + i ]= nL_calc[0][i];
    }
    for (int i = 0; i < 3; ++i) {
        dx[6 + NUM_FLEX_SEG*3 + 6 + i ]= x_coil[i+6];
    }
    for (int i = 0; i < 9; ++i) {
        dx[6 + NUM_FLEX_SEG*3 + 9 + i ]= x_coil[i+9];
    }



//    std::cout << "-------------------------------" << std::endl;
//        std::cout << "ActuationCurrents: " << ActuationCurrents[0][0] << " " << ActuationCurrents[0][1] << " " << ActuationCurrents[0][2] <<  std::endl;
//
//    std::cout <<" v : " << x_coil[0] << " " << x_coil[1] << " " << x_coil[2] << std::endl;
//    std::cout <<" w: " << x_coil[3] << " " << x_coil[4] << " " << x_coil[5] << std::endl;
//
//        std::cout << "pL_pre: " << x_coil[6] << " " << x_coil[7] << " " << x_coil[8] <<  std::endl;
//        std::cout << "RL: " << std::endl;
//        std::cout <<  x_coil[9] << " " <<x_coil[10] << " " << x_coil[11] <<  std::endl;
//        std::cout <<  x_coil[12] << " " << x_coil[13] << " " << x_coil[14] <<  std::endl;
//        std::cout <<  x_coil[15] << " " << x_coil[16] << " " << x_coil[17] <<  std::endl;


}


/** Output equation. **/
void compute_y(double *y, double t, double *x, double *u, double **p,
               const mxArray *auxvar)
{

    for (int i = 0; i < 3; ++i) {
        y[i] = x[i + 15]; //p
    }

//    for (int i = 0; i < 3; ++i) {
//        y[i+3] = x[i + 24]; //Rz
//    }

}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* Declaration of input and output arguments. */
    double *x, *u, **p, *dx, *y, *t;
    int     i, np;
    size_t  nu, nx;
    const mxArray *auxvar = NULL; /* Cell array of additional data. */

    if (nrhs < 3) {
        mexErrMsgIdAndTxt("IDNLGREY:ODE_FILE:InvalidSyntax",
                          "At least 3 inputs expected (t, u, x).");
    }

    /* Determine if auxiliary variables were passed as last input.  */
    if ((nrhs > 3) && (mxIsCell(prhs[nrhs-1]))) {
        /* Auxiliary variables were passed as input. */
        auxvar = prhs[nrhs-1];
        np = nrhs - 4; /* Number of parameters (could be 0). */
    } else {
        /* Auxiliary variables were not passed. */
        np = nrhs - 3; /* Number of parameters. */
    }

    /* Determine number of inputs and states. */
    nx = mxGetNumberOfElements(prhs[1]); /* Number of states. */
    nu = mxGetNumberOfElements(prhs[2]); /* Number of inputs. */

    /* Obtain double data pointers from mxArrays. */
    t = mxGetPr(prhs[0]);  /* Current time value (scalar). */
    x = mxGetPr(prhs[1]);  /* States at time t. */
    u = mxGetPr(prhs[2]);  /* Inputs at time t. */

    p = static_cast<double **>(mxCalloc(np, sizeof(double*)));
    for (i = 0; i < np; i++) {
        p[i] = mxGetPr(prhs[3+i]); /* Parameter arrays. */
    }

    /* Create matrix for the return arguments. */
    plhs[0] = mxCreateDoubleMatrix(nx, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(NY, 1, mxREAL);
    dx      = mxGetPr(plhs[0]); /* State derivative values. */
    y       = mxGetPr(plhs[1]); /* Output values. */

    /*
      Call the state and output update functions.

      Note: You may also pass other inputs that you might need,
      such as number of states (nx) and number of parameters (np).
      You may also omit unused inputs (such as auxvar).

      For example, you may want to use orders nx and nu, but not time (t)
      or auxiliary data (auxvar). You may write these functions as:
          compute_dx(dx, nx, nu, x, u, p);
          compute_y(y, nx, nu, x, u, p);
    */

    /* Call function for state derivative update. */
    compute_dx(dx, t[0], x, u, p, auxvar);

    /* Call function for output update. */
    compute_y(y, t[0], x, u, p, auxvar);

    /* Clean up. */
    mxFree(p);
}
