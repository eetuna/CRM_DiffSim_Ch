#pragma once
#include "mexAdapter.hpp"
#include "mex.hpp"

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
        double v_L_pre[NUM_ACT_SET][3], w_L_pre[NUM_ACT_SET][3],  u0_initialguess[NUM_FLEX_SEG][3], nL_initialguess[NUM_ACT_SET][3], mL_initialguess[NUM_ACT_SET][3], pL_pre[NUM_ACT_SET][3], RL_pre[NUM_ACT_SET][9];
        // Define initial guesses to be used when solving boundary value problem

        for (int i = 0; i < NUM_ACT_SET; ++i) {
            for (int j = 0; j < 3; ++j) {
                v_L_pre[i][j] = inputs[0][j+i*3];
                w_L_pre[i][j] = inputs[1][j+i*3];
            }
        }
        for (int i = 0; i < NUM_FLEX_SEG; ++i) {
            for (int j = 0; j < 3; ++j) {
                u0_initialguess[i][j] = inputs[2][j+3*i];
            }
        }

        for (int i = 0; i < NUM_ACT_SET; ++i) {
            for (int j = 0; j < 3; ++j) {
                mL_initialguess[i][j] = inputs[3][j+3*i];
                nL_initialguess[i][j] = inputs[4][j+3*i];
            }
        }

        for (int i = 0; i < NUM_ACT_SET; ++i) {
            for (int j = 0; j < 3; ++j) {pL_pre[i][j] = inputs[5][j+i*3];}
        }

        for (int i = 0; i < NUM_ACT_SET; ++i) {
            for (int j = 0; j < 9; ++j) {RL_pre[i][j] = inputs[6][j+i*9];}
        }

        // Declare the output variables for BVP
        // calculated curvature at the catheter base
        double u0_calc[NUM_FLEX_SEG][3];
        double nL_calc[NUM_ACT_SET][3];
        double mL_calc[NUM_ACT_SET][3];
        // calculated contstraint force at the catheter tip (this will be used when ContactMode == ContactModeType::FIXED_TIP)
        double ftip_calc[3];
        // numerical nonlinear equation solver diagnostic outputs
        int localmin;

        /** Retrieve u. **/
        double ActuationCurrents[NUM_ACT_SET][3];

        for (int i = 0; i < NUM_ACT_SET; i++)	for (int j = 0; j < 3; j++)	ActuationCurrents[i][j] = inputs[7][i * 3 + j];

        /** Retrieve model parameters. **/
        double damping_[NUM_ACT_SET][6];

        for (int i = 0; i < NUM_ACT_SET; ++i) {
            damping_[i][0] = damping_[i][1] = inputs[8][0+i*4]; damping_[i][2] =inputs[8][1+i*4];
            damping_[i][3] = damping_[i][4] = inputs[8][2+i*4]; damping_[i][5] =inputs[8][3+i*4];
        }

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
        double SegmentLengths[NUM_SEGMENTS] = { 19.85,  18.3, 59.40 };

        // Mass of each of the actuator sets - unit: ??
        double ActMass[NUM_ACT_SET] = { inputs[14][0]};// { 7.7736e-5};

        // Array of lambda values for marker locations  (distance from the tip to each of the markers) - unit: mm
        double MarkerLoc[NUM_LOCALIZATION_MARKERS] = { 2.18, 18.79, 50.06, 101.52, 104.02 };
        // Length density for each of the catheter segments
        double rho[NUM_SEGMENTS] ={ 2.8814e-7, 2.8814e-7, 2.8814e-7 }; // {1.1701e-6, 1.1701e-6, 1.1701e-6}; 
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
        double residual[NUM_RESIDUAL];

        std::cout << "ActuationCurrents: " << ActuationCurrents[0][0] << " " << ActuationCurrents[0][1] << " " << ActuationCurrents[0][2] <<  std::endl;

        double InsertedLength =  0.0;//98.5; // u[NUM_CONTROL -1 ];
        for (int i = 0; i < NUM_SEGMENTS; ++i) {
            InsertedLength += SegmentLengths[i]; // we are not controlling this now
        }


        CRMConstructShootingMethodParamSet<double>(CathParams, CathConfig, InsertedLength, ActuationCurrents, ContactMode,
                                                   TipConstraintPoint, TipForce, IntegrationStepSize,
                                                   v_L_pre, w_L_pre, pL_pre,  RL_pre, damping_, Delta_T,
                                                   BVPParams);


        DynamicsBVP(BVPParams, u0_initialguess, mL_initialguess, nL_initialguess, ftip_initialguess,
                    u0_calc, mL_calc, nL_calc, ftip_calc, localmin);

//        std::cout << "out_u0: " << u0_calc[0] << " " << u0_calc[1] << " " << u0_calc[2] <<  std::endl;
//        std::cout << "out_mL: " << mL_calc[0] << " " << mL_calc[1] << " " << mL_calc[2] <<  std::endl;
//        std::cout << "out_nL: " << nL_calc[0] << " " << nL_calc[1] << " " << nL_calc[2] <<  std::endl;

        double out_ReportedMarkerPos[NUM_LOCALIZATION_MARKERS][3], x_coil[NUM_COIL_STATES];

        DYNSolverIVP(BVPParams, u0_calc, mL_calc, nL_calc, ftip_calc,
                     true, xf, x_coil,out_ReportedMarkerPos);

        /** Output report **/
        outputs[0] = factory.createArray<double>({1, 3}, {x_coil[0], x_coil[1], x_coil[2]});
        outputs[1] = factory.createArray<double>({1, 3}, {x_coil[3], x_coil[4], x_coil[5]});
        outputs[2] = factory.createArray<double>({1, NUM_FLEX_SEG*3}, {u0_calc[0][0], u0_calc[0][1], u0_calc[0][2], u0_calc[1][0], u0_calc[1][1], u0_calc[1][2]});
        outputs[3] = factory.createArray<double>({1, 3}, {mL_calc[0][0], mL_calc[0][1], mL_calc[0][2]});
        outputs[4] = factory.createArray<double>({1, 3}, {nL_calc[0][0], nL_calc[0][1], nL_calc[0][2]});
        outputs[5] = factory.createArray<double>({1, 3}, {x_coil[6], x_coil[7], x_coil[8]});
        outputs[6] = factory.createArray<double>({1, 9}, {x_coil[9], x_coil[10], x_coil[11],x_coil[12], x_coil[13], x_coil[14],x_coil[15], x_coil[16], x_coil[17]});
        outputs[7] = factory.createArray<double>({1, 3}, {xf[12], xf[13], xf[14]});

    }

};

