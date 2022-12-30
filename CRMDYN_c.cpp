#include <iostream>
#include <chrono>
#include "CRMDYN.hpp"
using namespace std::chrono;

#define EPS 1e-5

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
    double v_L_pre[3], w_L_pre[3], pL_pre[3], RL_pre[9];
    // Define initial guesses to be used when solving boundary value problem
    for (int i = 0; i < 3; ++i) {
        v_L_pre[i] = x[i];
        w_L_pre[i] = x[i+3];
    }

    double initial_guesses[NUM_RESIDUAL];
    for (int i = 0; i < 3; ++i) {
        initial_guesses[i] = x[i+6];
        initial_guesses[i+3] = x[i+9];
    }
    if (NUM_RESIDUAL > 6){
        for (int i = 0; i < 3; ++i) {
            initial_guesses[i+6] = x[i+12];
            initial_guesses[i+9] = x[i+15];
        }
    }

    for (int i = 0; i < 3; ++i) {
        pL_pre[i] = x[i+18];
    }
    for (int i = 0; i < 9; ++i) {
        RL_pre[i] = x[i+21];
    }

//    double h0[NUM_FLEX_SEG] = {x[30], x[31]};


//std:: cout << "SegSteps_test " << SegSteps_test[0] << " " << SegSteps_test[1] << std::endl;
//        std:: cout << "h0_test " << h0_test[0] << " " << h0_test[1] << std::endl;

    int SegSteps[NUM_FLEX_SEG] = {(int)p[0][0], (int)p[0][1]};

    int length_0 = ( SegSteps[0] + 1 ) * 3 ;

    paramHistory u_history[NUM_FLEX_SEG], v_history[NUM_FLEX_SEG], w_history[NUM_FLEX_SEG];
    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
        int length = (SegSteps[i] + 1) * 3;
        u_history[i].length = length;
        for (int j = 0; j < (SegSteps[i] + 1); ++j) {
            u_history[i].data[j * 3] = x[30+ i * length_0 + j * 3];
            u_history[i].data[j * 3 + 1] = x[30+ i * length_0 + j * 3+1];
            u_history[i].data[j * 3 + 2] = x[30+ i * length_0 + j * 3+2];
        }
    }
    int length_out = ( SegSteps[0] + 1 ) * 3 + ( SegSteps[1] + 1 ) * 3;
    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
        int length = (SegSteps[i] + 1) * 3;
        v_history[i].length = length;
        for (int j = 0; j < (SegSteps[i] + 1); ++j) {
            v_history[i].data[j * 3] = x[30+length_out + i * length_0 + j * 3];
            v_history[i].data[j * 3 + 1] = x[30+length_out + i * length_0 + j * 3 + 1];
            v_history[i].data[j * 3 + 2] = x[30+length_out + i * length_0 + j * 3 + 2];
        }
    }
    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
        int length = (SegSteps[i] + 1) * 3;
        w_history[i].length = length;
        for (int j = 0; j < (SegSteps[i] + 1); ++j) {
            w_history[i].data[j * 3] = x[30+ 2* length_out + i * length_0 + j * 3];
            w_history[i].data[j * 3 + 1] = x[30+ 2* length_out + i * length_0 + j * 3 + 1];
            w_history[i].data[j * 3 + 2] = x[30+ 2* length_out + i * length_0 + j * 3 + 2];
        }
    }

    /** Retrieve u. **/
    double ActuationCurrents[NUM_ACT_SET][3];
    for (int i = 0; i < NUM_ACT_SET; i++)	for (int j = 0; j < 3; j++)	ActuationCurrents[i][j] = u[i * 3 + j];

    /** Retrieve model parameters. **/
    double damping_tubing[3], damping_coil[6];
    for (int i = 0; i < 3; ++i) {
        damping_tubing[i] = p[1][i];
    }

    damping_coil[0] = damping_coil[1] = p[2][0]; damping_coil[2] =p[2][1];
    damping_coil[3] = damping_coil[4] = p[2][2]; damping_coil[5] =p[2][3];


    double Delta_T = p[3][0];

    double oRlist[NUM_FLEX_SEG] ={p[4][0], p[4][0]}; // { 1.5875, 1.5875 };
    // Inner radii of each of the flexible segments - unit: mm
    double iRlist[NUM_FLEX_SEG] ={p[4][1], p[4][1]};// { 0.9906, 0.9906 };
    // Young's moduli of each of the flexible segments - unit: ??
    double YoungModlist[NUM_FLEX_SEG] = {p[5][0], p[5][0]};// { 5.3948, 5.3948 };// { 4.9105, 4.9105 };
    // Shear moduli of each of the flexbile segments - unit: ??
    double ShearModlist[NUM_FLEX_SEG] = {p[5][1], p[5][1]};//{ 2.3881, 2.3881 };//{ 1.6545, 1.6545 };
    // Alignment angles for the side coils, for each of the actuator sets - unit: radians
    double CoilAlignmentAngles[NUM_ACT_SET][2] = {p[6][0], p[6][1]}; //{ {0.0, 0.0}  };// { {-0.0871, -0.3934}  };
    // Coil turn area matrices for each of the actuator sets  - unit: mm2
    double CoilTurnAreaMat[NUM_ACT_SET][9] ={ { p[7][0], 0.0, 0.0, 0.0, p[7][1], 0.0, 0.0, 0.0, p[7][2] } };;// { { 1.44, 0.0, 0.0, 0.0, 1.3851, 0.0, 0.0, 0.0, 1.60 } };
    // Lengths of each of the catheter segments - unit: mm
    double SegmentLengths[NUM_SEGMENTS] = {p[8][0], p[8][1], p[8][2]};// { 19.85,  18.3, 59.40 };

    // Mass of each of the actuator sets - unit: ??
    double ActMass[NUM_ACT_SET] = { p[9][0]};// { 7.7736e-5};

    // *** Numerical Computation Parameters
    // Stepsize used in numerical integration along the length of the catheter during IVP - unit: mm
    double IntegrationStepSize = p[10][0];// 0.2;

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


    // we are adding the tubing mass of the coil section to the total mass of actuator:
    double ActInertia[NUM_ACT_SET][9];
    for (int i = 0; i < NUM_ACT_SET; ++i)
    {
        double I_zz = 0.5 * (ActMass[i]) * (oRlist[0] * oRlist[0] + iRlist[0] * iRlist[0]);
        double I_xx = 0.25 * (ActMass[i]) * (oRlist[0] * oRlist[0] + iRlist[0] * iRlist[0]) + 1.0 / 12 * (ActMass[i]) * SegmentLengths[2*i+1] * SegmentLengths[2*i+1];
        ActInertia[i][0] = I_xx; ActInertia[i][1] = 0.0; ActInertia[i][2] = 0.0;
        ActInertia[i][3] = 0.0; ActInertia[i][4] = I_xx; ActInertia[i][5] = 0.0;
        ActInertia[i][6] = 0.0; ActInertia[i][7] = 0.0; ActInertia[i][8] = I_zz;
    }

    /**
      * initial curvature is assumed constant curvature, we are initializing the full length u_history,
      * Prepare for interpolation of dynamic insertion length
      */
    double v0[3] = {0.0,0.0,0.0};              //initial linear velocity in local frame of entry point
    double w0[3] = {0.0,0.0,0.0};

    // Area Moment of Inertia for a hollow cylindrical section
    // https://en.wikipedia.org/wiki/Second_moment_of_area
    // https://www.engineeringtoolbox.com/area-moment-inertia-d_1328.html
    double tubingInertia[NUM_FLEX_SEG][9];
    for (int i = 0; i < NUM_FLEX_SEG; ++i)
    {
        double I_xx =  M_PI * (oRlist[i] * oRlist[i] * oRlist[i] * oRlist[i] - iRlist[i] * iRlist[i] * iRlist[i] * iRlist[i] )  * 0.25;
        double I_zz = 2 * I_xx;
        tubingInertia[i][0] = I_xx; tubingInertia[i][1] = 0.0; tubingInertia[i][2] = 0.0;
        tubingInertia[i][3] = 0.0; tubingInertia[i][4] = I_xx; tubingInertia[i][5] = 0.0;
        tubingInertia[i][6] = 0.0; tubingInertia[i][7] = 0.0; tubingInertia[i][8] = I_zz;
    }

//    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
//        for (int j = 0; j < ( SegSteps[i] + 1 ) ; ++j) {
//            std::cout << " w_history " << w_history[i].data[j*3] << " " << w_history[i].data[j*3+1] << " " << w_history[i].data[j*3+2] << std::endl;
//            std::cout <<"end of first round: " << j*3+2 << std::endl;
//        }
//        std::cout << "length_of w: " <<  w_history[i].length << std::endl;
//
//    }

    // Define variables of convenience that will be use to package catheter physical parameters and spatial configuration parameters
    CRMCatheterModelParams CathParams;		// Catheter physical parameters will be packaged into this structure
    CatheterConfiguration CathConfig;		// Catheter spatial configuration parameters will be packages into this structure

    // Package catheter physical parameters and spatial configuration parameters
    //   this step would typically needs to be executed only once
    CRMShootingMethodBVP_Prep(B0,gravity, p0, R0, v0, w0,
                              SegmentLengths, MarkerLoc, iRlist, oRlist, YoungModlist, ShearModlist, ustarlist,
                              CoilAlignmentAngles, CoilTurnAreaMat, rho, ActMass, ActInertia, tubingInertia,
                              CathParams, CathConfig);

    double InsertedLength  = 0.0;
    for (int i = 0; i < NUM_SEGMENTS; ++i) {
        InsertedLength += SegmentLengths[i];
    }
    /** Compute next state **/
    CRMShootingMethodParams<double> BVPParams{};

    double ftip_initialguess[3] = { 0.0, 0.0, 0.0 };
    // Declare the output variables for BVP
    double out_calc[NUM_RESIDUAL];
    // calculated contstraint force at the catheter tip (this will be used when ContactMode == ContactModeType::FIXED_TIP)
    double ftip_calc[3];
    // numerical nonlinear equation solver diagnostic outputs
    int localmin;
    // Declare output variables
    // catheter shape state at the entry point of the catheter
    //   states are packed u[0..2], R[0..8], p[0..2](R: 3x3 matrix stored in row major order R11 R12 R13 R21 R22 R23 R31 R32 R33)
    double xf[NUM_STATES];
    // moment residual at the catheter tip - this should converge to {0,0,0} if the catheter is at its equilibrium configuration
    double residual[NUM_RESIDUAL];
    // spatial coordinates of the localization markers
    double ReportedMarkerPos[NUM_LOCALIZATION_MARKERS][3];

    CRMConstructShootingMethodParamSet<double>(CathParams, CathConfig, InsertedLength, ActuationCurrents, ContactMode, TipConstraintPoint, TipForce,
                                               IntegrationStepSize, v_L_pre, w_L_pre, u_history, v_history, w_history, Delta_T, damping_tubing, damping_coil, BVPParams);

    // Cosserat Rod Model - Solve the Boundary Value Problem to calculate the equilibrium configuration of the catheter
    CRMShootingMethodBVP(BVPParams, initial_guesses, ftip_initialguess,
                         out_calc, ftip_calc, localmin);

    paramHistory u_history_update[NUM_FLEX_SEG], v_history_update[NUM_FLEX_SEG], w_history_update[NUM_FLEX_SEG];
    double v_L_update[3], w_L_update[3], pL_update[3], RL_update[9], h0_new[NUM_FLEX_SEG];

    // Cosserat Rod Model - Solve the Initial Value Problem to calculate the shape of the catheter
    CRMSolverIVP(BVPParams, out_calc, ftip_calc, true, xf, residual, ReportedMarkerPos, u_history_update, v_history_update, w_history_update,
                 v_L_update, w_L_update, pL_update, RL_update, h0_new);

//
//        std::cout << "ActuationCurrents: " << ActuationCurrents[0][0] << " " << ActuationCurrents[0][1] << " " << ActuationCurrents[0][2] <<  std::endl;
//        std::cout << "pL_pre: " << pL_pre[0] << " " << pL_pre[1] << " " << pL_pre[2] <<  std::endl;
//        std::cout << "RL: " << std::endl;
//        std::cout <<  RL_pre[0] << " " << RL_pre[1] << " " << RL_pre[2] <<  std::endl;
//        std::cout <<  RL_pre[3] << " " << RL_pre[4] << " " << RL_pre[5] <<  std::endl;
//        std::cout <<  RL_pre[6] << " " << RL_pre[7] << " " << RL_pre[8] <<  std::endl;

    /** Output report **/
    for (int i = 0; i < 3; ++i) {
        if(i < 3){   dx[i] = v_L_update[i]; } //v, w
        else if (i < 6) dx[i+3] = w_L_update[i - 3];
        else if (i < 18) dx[i+6] = out_calc[i - 6]; //curvatures and internal forces
        else if (i < 21) dx[i+18] = pL_update[i - 18];
        else if (i < 30) dx[i+21] = RL_update[i - 21]; //p
    }
//    dx[30] = h0_new[0]; dx[31] = h0_new[1];

    double param_u[length_out], param_v[length_out], param_w[length_out];
    int ind_seg = 0;
    for (int i = 0; i < SegSteps[ind_seg] + 1; ++i) {
        param_u[ 3* i ] = u_history_update[ind_seg].data[i*3];
        param_u[ 3* i + 1] = u_history_update[ind_seg].data[i*3+1];
        param_u[ 3* i + 2] = u_history_update[ind_seg].data[i*3+2];

        param_v[ 3* i ] = v_history_update[ind_seg].data[i*3];
        param_v[ 3* i + 1] = v_history_update[ind_seg].data[i*3+1];
        param_v[ 3* i + 2] = v_history_update[ind_seg].data[i*3+2];

        param_w[ 3* i ] = w_history_update[ind_seg].data[i*3];
        param_w[ 3* i + 1] = w_history_update[ind_seg].data[i*3+1];
        param_w[ 3* i + 2] = w_history_update[ind_seg].data[i*3+2];

    }
    ind_seg=1;
    for (int i = 0; i < SegSteps[ind_seg] + 1; ++i) {
        param_u[ length_0 + 3* i ] = u_history_update[ind_seg].data[i*3];
        param_u[ length_0 + 3* i + 1] = u_history_update[ind_seg].data[i*3+1];
        param_u[ length_0 + 3* i + 2] = u_history_update[ind_seg].data[i*3+2];

        param_v[ length_0+ 3* i ] = v_history_update[ind_seg].data[i*3];
        param_v[ length_0+ 3* i + 1] = v_history_update[ind_seg].data[i*3+1];
        param_v[ length_0+ 3* i + 2] = v_history_update[ind_seg].data[i*3+2];

        param_w[ length_0 + 3* i ] = w_history_update[ind_seg].data[i*3];
        param_w[ length_0 + 3* i + 1] = w_history_update[ind_seg].data[i*3+1];
        param_w[ length_0 + 3* i + 2] = w_history_update[ind_seg].data[i*3+2];
    }

    for (int i = 0; i <length_out; ++i) {dx[i+ 30] = param_u[i];}
    for (int i = 0; i <length_out ; ++i) {dx[i+ 30 + length_out] = param_v[i];}
    for (int i = 0; i <length_out ; ++i) {dx[i+ 30+ length_out *2] = param_w[i];}

//        for (int i = 0; i < NUM_FLEX_SEG; ++i) {
//            for (int j = 0; j < ( SegSteps[i] + 1 ) ; ++j) {
//                std::cout << " w_history_update " << w_history_update[i].data[j*3] << " " << w_history_update[i].data[j*3+1] << " " << w_history_update[i].data[j*3+2] << std::endl;
//                std::cout <<"end of first round: " << j*3+2 << std::endl;
//            }
//            std::cout << "length_of w_history_update : " <<  w_history_update[i].length << std::endl;
//
//        }

}


/** Output equation. **/
void compute_y(double *y, double t, double *x, double *u, double **p,
               const mxArray *auxvar)
{

    for (int i = 0; i < 3; ++i) {
        y[i] = x[i + 18]; //p
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

////test before run on Matlab
//int main(int argc, char** argv){
//
//    double v_L_pre[3] = { 0.0, 0.0, 0.0 };
//    double w_L_pre[3] = { 0.0, 0.0, 0.0 };
//    // Define initial guesses to be used when solving boundary value problem
//    double u0_initialguess[3] = { 0.0, 0.0, 0.0};
//    // initial guess for the contstraint force at the catheter tip (this will be used when ContactMode == ContactModeType::FIXED_TIP)
//    double nL_initialguess[3] = { 0.0, 0.0, 0.0 };
//    double mL_initialguess[3] = { 0.0, 0.0, 0.0 };
//
//
//    double p_0[3] = {-0.407795, -0.301721, 69.4976};
//    double R_0[9] = {0.999946, -3.98697e-05, -0.0103811,
//                     -3.98697e-05, 0.999971, -0.00768085,
//                    0.0103811, 0.00768085, 0.999917};
//
//
//    double x[NUM_DYN_STATE],  dx[NUM_DYN_STATE], u[NUM_CONTROL];
//    double **p;
//
//    p = (double **) malloc(1 * sizeof(double)); //just damping gmat for now
//    p[0] = (double *) malloc(6*sizeof(double));
//
//    for (int i = 0; i < NUM_DYN_STATE; ++i) {
//        if(i < 3) x[i] = v_L_pre[i];
//        else if (i < 6) x[i] = w_L_pre[i-3];
//        else if (i < 9) x[i] = u0_initialguess[i-6];
//        else if (i < 12) x[i] = nL_initialguess[i-9];
//        else if (i < 15) x[i] = mL_initialguess[i-12];
//        else if (i < 18) x[i] = p_0[i-15];
//        else x[i] = R_0[i-18];
//    }
//
//    for (int i = 0; i < 3; ++i) {
//        u[i] = 0.0;
//    }
//    u[NUM_CONTROL-1] = 98.5;
//
//    double t = 0;
//    for (int i = 0; i < 3; ++i) {
//        p[0][i] = 1.0; //v
//    }
//    for (int i = 0; i < 3; ++i) {
//        p[0][i+3] = 0.05; //w
//    }
//    const mxArray *auxvar = NULL;
//
//    compute_dx(dx, t, x, u, static_cast<double **>(p), auxvar);
//
//    std::cout << "pL: " << dx[15] << " " << dx[16] << " " << dx[17] <<  std::endl;
//
//    std::cout << "RL: " << std::endl;
//    std::cout <<  dx[18] << " " << dx[19] << " " << dx[20] <<  std::endl;
//    std::cout <<  dx[21] << " " << dx[22] << " " << dx[23] <<  std::endl;
//    std::cout <<  dx[24] << " " << dx[25] << " " << dx[26] <<  std::endl;
//
//
//
//
//}