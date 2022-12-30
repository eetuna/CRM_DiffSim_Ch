#pragma once
#include "/usr/local/MATLAB/R2022b/extern/include/mexAdapter.hpp"
#include "/usr/local/MATLAB/R2022b/extern/include/mex.hpp"

#include <cmath>
#include "CRMDYN.hpp"

#define M_PI 3.14159265358979323846
#define _USE_MATH_DEFINES


using matlab::mex::ArgumentList; // added
using namespace matlab::data;	 // added


//#define MIN(a, b) (((a) < (b)) ? (a) : (b))
//#define POW4(a) ((a) * (a) * (a) * (a))

class MexFunction : public matlab::mex::Function {
    ArrayFactory factory;

public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
    {
        // *** Physical Description of the Catheter
        // IMPORTANT NOTE: For now, most proximal segment is assumed to be always flexible
        //    and the flexible and rigid segments are assumed to be alternating
        //    most distal segment can be flexible or rigid
        // Items are listed in distal-to-proximal order (starting from the tip of the catheter towards the base)
        // Outer radii of each of the flexible segments - unit: mm
        /** Retrieve everything from inputs **/
        /** Retrieve x. **/
        double v_L_pre[3], w_L_pre[3];
        // Define initial guesses to be used when solving boundary value problem
        for (int i = 0; i < 3; ++i) {
            v_L_pre[i] = inputs[0][i];
            w_L_pre[i] = inputs[1][i];
        }

        double initial_guesses[NUM_RESIDUAL];
        for (int i = 0; i < 3; ++i) {
            initial_guesses[i] = inputs[2][i];
            initial_guesses[i+3] = inputs[3][i];
        }
        if (NUM_RESIDUAL > 6){
            for (int i = 0; i < 3; ++i) {
                initial_guesses[i+6] = inputs[4][i];
                initial_guesses[i+9] = inputs[5][i];
            }
        }

        /** Retrieve u. **/
        double ActuationCurrents[NUM_ACT_SET][3];
        for (int i = 0; i < NUM_ACT_SET; i++)	for (int j = 0; j < 3; j++)	ActuationCurrents[i][j] = inputs[6][i * 3 + j];

        /** Retrieve model parameters. **/
        double damping_tubing[3], damping_coil[6];
        for (int i = 0; i < 3; ++i) {
            damping_tubing[i] = inputs[7][i];
        }

        damping_coil[0] = damping_coil[1] = inputs[8][0]; damping_coil[2] =inputs[8][1];
        damping_coil[3] = damping_coil[4] = inputs[8][2]; damping_coil[5] =inputs[8][3];

        double Delta_T = inputs[9][0];

        double oRlist[NUM_FLEX_SEG] ={inputs[10][0], inputs[10][0]}; // { 1.5875, 1.5875 };
        // Inner radii of each of the flexible segments - unit: mm
        double iRlist[NUM_FLEX_SEG] ={inputs[10][1], inputs[10][1]};// { 0.9906, 0.9906 };
        // Young's moduli of each of the flexible segments - unit: ??
        double YoungModlist[NUM_FLEX_SEG] = {inputs[11][0], inputs[11][0]};// { 5.3948, 5.3948 };// { 4.9105, 4.9105 };
        // Shear moduli of each of the flexbile segments - unit: ??
        double ShearModlist[NUM_FLEX_SEG] = {inputs[11][1], inputs[11][1]};//{ 2.3881, 2.3881 };//{ 1.6545, 1.6545 };
        // Alignment angles for the side coils, for each of the actuator sets - unit: radians
        double CoilAlignmentAngles[NUM_ACT_SET][2] = {inputs[12][0], inputs[12][1]}; //{ {0.0, 0.0}  };// { {-0.0871, -0.3934}  };
        // Coil turn area matrices for each of the actuator sets  - unit: mm2
        double CoilTurnAreaMat[NUM_ACT_SET][9] ={ { inputs[13][0], 0.0, 0.0, 0.0, inputs[13][1], 0.0, 0.0, 0.0, inputs[13][2] } };;// { { 1.44, 0.0, 0.0, 0.0, 1.3851, 0.0, 0.0, 0.0, 1.60 } };
        // Lengths of each of the catheter segments - unit: mm
        double SegmentLengths[NUM_SEGMENTS] = {inputs[14][0], inputs[14][1], inputs[14][2]};// { 19.85,  18.3, 59.40 };

        // Mass of each of the actuator sets - unit: ??
        double ActMass[NUM_ACT_SET] = { inputs[15][0]};// { 7.7736e-5};

        // *** Numerical Computation Parameters
        // Stepsize used in numerical integration along the length of the catheter during IVP - unit: mm
        double IntegrationStepSize = inputs[16][0];// 0.2;

        // Array of lambda values for marker locations  (distance from the tip to each of the markers) - unit: mm
        double MarkerLoc[NUM_LOCALIZATION_MARKERS] = { 2.18, 18.79, 50.06, 101.52, 104.02 };
        // Length density for each of the catheter segments
        double rho[NUM_SEGMENTS] = { 2.8814e-7, 2.8814e-7, 2.8814e-7 };
        // Local curvature in unloaded configuration for each of the flexible segments; NUM_FLEX_SEG*3 long array corresponding to NUM_FLEX_SEG many 3x1 vectors
        double ustarlist[NUM_FLEX_SEG][3] = { 0.000148626102272827, 0.00094448815853795, 0, 0.000799849479553773, -0.0001892201697658481, 0};//0.00567, -0.00822, 0.0};

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


        int SegSteps[NUM_FLEX_SEG] = {inputs[17][0], inputs[17][1]};
//        double h0[NUM_FLEX_SEG] = {inputs[18][0], inputs[18][1]};

        int length_0 = ( SegSteps[0] + 1 ) * 3 ;

//std:: cout << "SegSteps_test " << SegSteps_test[0] << " " << SegSteps_test[1] << std::endl;
//        std:: cout << "h0_test " << h0_test[0] << " " << h0_test[1] << std::endl;

        paramHistory u_history[NUM_FLEX_SEG], v_history[NUM_FLEX_SEG], w_history[NUM_FLEX_SEG];
        for (int i = 0; i < NUM_FLEX_SEG; ++i) {
            int length = ( SegSteps[i] + 1 ) * 3 ;
            u_history[i].length = length;
            for (int j = 0; j < ( SegSteps[i] + 1 ) ; ++j) {
                u_history[i].data[j*3] = inputs[18][i*length_0 + j*3];
                u_history[i].data[j*3+1] = inputs[18][i*length_0 + j*3+1];
                u_history[i].data[j*3+2] = inputs[18][i*length_0 + j*3+2];
            }

            v_history[i].length = length;
            for (int j = 0; j < ( SegSteps[i] + 1 ) ; ++j) {
                v_history[i].data[j*3] = inputs[19][i*length_0 + j*3];
                v_history[i].data[j*3+1] = inputs[19][i*length_0 + j*3+1];
                v_history[i].data[j*3+2] = inputs[19][i*length_0 + j*3+2];
            }

            w_history[i].length = length;
            for (int j = 0; j < ( SegSteps[i] + 1 ) ; ++j) {
                w_history[i].data[j*3] = inputs[20][i*length_0 + j*3];
                w_history[i].data[j*3+1] = inputs[20][i*length_0 + j*3+1];
                w_history[i].data[j*3+2] = inputs[20][i*length_0 + j*3+2];
            }
        }

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
        outputs[0] = factory.createArray<double>({1, 3}, {v_L_update[0], v_L_update[1], v_L_update[2]});
        outputs[1] = factory.createArray<double>({1, 3}, {w_L_update[0], w_L_update[1], w_L_update[2]});
        outputs[2] = factory.createArray<double>({1, 3}, {out_calc[0], out_calc[1], out_calc[2]}); //u0
        outputs[3] = factory.createArray<double>({1, 3}, {out_calc[3], out_calc[4], out_calc[5]}); //n0
        outputs[4] = factory.createArray<double>({1, 3}, {out_calc[6], out_calc[7], out_calc[8]}); // uL
        outputs[5] = factory.createArray<double>({1, 3}, {out_calc[9], out_calc[10], out_calc[11]}); //nL
        outputs[6] = factory.createArray<double>({1, 3}, {pL_update[0], pL_update[1], pL_update[2]}); //number of flexible segments
        outputs[7] = factory.createArray<double>({1, 9}, {RL_update[0], RL_update[1],RL_update[2],RL_update[3],RL_update[4],RL_update[5],RL_update[6],RL_update[7],RL_update[8]}); //number of flexible segments
//        outputs[8] = factory.createArray<double>({1, 2}, {h0_new[0], h0_new[1]}); //number of flexible segments


        size_t length_out = (SegSteps[0] + 1)*3 + (SegSteps[1] + 1)*3;

        TypedArray<double> param_u = factory.createArray<double>({1, length_out }); //
        TypedArray<double> param_v = factory.createArray<double>({1, length_out }); //
        TypedArray<double> param_w = factory.createArray<double>({1, length_out }); //

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

        outputs[8] = param_u;
        outputs[9] = param_v;
        outputs[10] = param_w;

    }

};

