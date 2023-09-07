#include <iostream>
#include <cmath>
#include <chrono> 
#include<string.h>

#include "CRMTest.h"
#include "CRMDYN.hpp"
using namespace std::chrono;

#define EPS 1e-5


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

template <int D1, int D2>
bool MatrixEqual(double A[D1][D2], double B[D1][D2], double eps) {
	for (int i=0; i<D1;i++)
		for (int j=0; j<D2; j++)
			if (fabs(A[i][j]-B[i][j])<eps) {} else return false;
	return true;
}

void wait_for_enter(const std::string &msg) {
    std::cout << msg << std::endl;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

int RunTests(void);
int RunExample(void);
int RunDynamicsExample(std::vector<std::vector<double> > &output_tip_pts);


int main(int argc, char** argv) {

	RunExample();

}

int RunExample(void) {

	std::cout << "### CRM Forward Kinematics Examples... " << std::endl;
	std::cout << std::endl << "Free Space Deflection Example: " << std::endl << std::endl;

//    double youngs =  8.2224;
//    double shears = 1.7553;
//
//    // *** Physical Description of the Catheter
//    // IMPORTANT NOTE: For now, most proximal segment is assumed to be always flexible
//    //    and the flexible and rigid segments are assumed to be alternating
//    //    most distal segment can be flexible or rigid
//    // Items are listed in distal-to-proximal order (starting from the tip of the catheter towards the base)
//    // Outer radii of each of the flexible segments - unit: mm
//    double oRlist[NUM_FLEX_SEG] = { 1.5875, 1.5875 };
//    // Inner radii of each of the flexible segments - unit: mm
//    double iRlist[NUM_FLEX_SEG] = { 0.9906, 0.9906 };
//    // Young's moduli of each of the flexible segments - unit: Mpa
//    double YoungModlist[NUM_FLEX_SEG] =  { youngs, youngs };//
//    // Shear moduli of each of the flexbile segments - unit: Mpa
//    double ShearModlist[NUM_FLEX_SEG] = { shears, shears};//
//    // Alignment angles for the side coils, for each of the actuator sets - unit: radians
//    double CoilAlignmentAngles[NUM_ACT_SET][2] = { { -0.1190, -3.0760}   };// { {-0.0871, -0.3934}  };
//    // Coil turn area matrices for each of the actuator sets  - unit: mm2
//    double CoilTurnAreaMat[NUM_ACT_SET][9] = { { 1.8367,  0.0, 0.0, 0.0, 2.8032, 0.0, 0.0, 0.0,  1.8448}};
//    // Lengths of each of the catheter segments - unit: mm
//    double SegmentLengths[NUM_SEGMENTS] = {19.85,  18.3, 159.40};// { 18.72, 16.14, 109.40};
//
//    // Mass of each of the actuator sets - unit: kg
//    double ActMass[NUM_ACT_SET] = { 8.2859e-06};
//    // Array of lambda values for marker locations  (distance from the tip to each of the markers) - unit: mm
//    double MarkerLoc[NUM_LOCALIZATION_MARKERS] = { 2.18, 18.79, 50.06, 101.52, 104.02 };
//    // Length density for each of the catheter segments
//    double rho[NUM_SEGMENTS] = { 2.8814e-7, 2.8814e-7, 2.8814e-7 };
//    // unit (1/ mm) Local curvature in unloaded configuration for each of the flexible segments; NUM_FLEX_SEG*3 long array corresponding to NUM_FLEX_SEG many 3x1 vectors
//    double ustarlist[NUM_FLEX_SEG][3] = { 0.0000148626102272827, 0.000094448815853795, 0, 0.0007399849479553773, -0.0002292201697658481, 0 };//0.00567, -0.00822, 0.0};
////    double ustarlist[NUM_FLEX_SEG][3] = { 0.0, 0.0, 0, 0.0, -0.0, 0};//0.00567, -0.00822, 0.0};
//
//
//    double InsertedLength  = 70.03562;
//    // Actuation currents for each of the coils for each of the coil sets - unit: A
//    double ActuationCurrents[NUM_ACT_SET][3] = {{ 0.0182649, -0.28875, 0.528}};
//
//
//    double damping[NUM_ACT_SET][6] = {{12.1761626666366,12.1761626666366,
//                                       284.429938756989,
//                                       0.0304776127617393,0.0304776127617393,
//                                       0.00502712804532508}};
//    double DELTA_T = 0.05;
//
//double xf[NUM_STATES] = {      0.0148634,0.00944569,-1.33222e-07,
//                                         0.984631,
//                                         0.0255816,
//                                         0.172762,
//                                         0.0178539,
//                                         0.969287,
//                                         -0.245282,
//                                         -0.173731,
//                                         0.244597,
//                                         0.953934,
//                                         1.08571,
//    0.231991,
//    97.1879
//    };
//    double pL[NUM_ACT_SET][3] = {-0.432543, 1.80066, 68.467 };
//    double RL[NUM_ACT_SET][9] = {0.999937, 0.00149663, -0.011162,
//                                 -0.000976585, 0.99892, 0.0464513,
//                                 0.0112194 ,-0.0464374, 0.998858};
//    double v_L_pre[NUM_ACT_SET][3] = {{ 0.0, 0.0, 0.0 }};
//    double w_L_pre[NUM_ACT_SET][3] = {{ 0.0, 0.0, 0.0 }};
//    // Define initial guesses to be used when solving boundary value problem
//    double nL_initialguess[NUM_ACT_SET][3] = {{ 0.0, 0.0, 0.0 }};
//    double mL_initialguess[NUM_ACT_SET][3] = {{ 0.0, 0.0, 0.0 }};


    /**
     * 2 coil sets
     */
    double youngs =  8.2224;
    double shears = 1.7553;

    // *** Physical Description of the Catheter
    // IMPORTANT NOTE: For now, most proximal segment is assumed to be always flexible
    //    and the flexible and rigid segments are assumed to be alternating
    //    most distal segment can be flexible or rigid
    // Items are listed in distal-to-proximal order (starting from the tip of the catheter towards the base)
    // Outer radii of each of the flexible segments - unit: mm
    double oRlist[NUM_FLEX_SEG] = { 1.5875, 1.5875, 1.5875 };
    // Inner radii of each of the flexible segments - unit: mm
    double iRlist[NUM_FLEX_SEG] = { 0.9906, 0.9906, 0.9906 };
    // Young's moduli of each of the flexible segments - unit: Mpa
    double YoungModlist[NUM_FLEX_SEG] =  { youngs, youngs, youngs };//
    // Shear moduli of each of the flexbile segments - unit: Mpa
    double ShearModlist[NUM_FLEX_SEG] = { shears, shears, shears};//
    // Alignment angles for the side coils, for each of the actuator sets - unit: radians
    double CoilAlignmentAngles[NUM_ACT_SET][2] = { { -0.1190, -3.0760},  { -0.1190, -3.0760}    };// { {-0.0871, -0.3934}  };
    // Coil turn area matrices for each of the actuator sets  - unit: mm2
    double CoilTurnAreaMat[NUM_ACT_SET][9] = { { 1.8367,  0.0, 0.0, 0.0, 2.8032, 0.0, 0.0, 0.0,  1.8448}, { 1.8367,  0.0, 0.0, 0.0, 2.8032, 0.0, 0.0, 0.0,  1.8448}};
    // Lengths of each of the catheter segments - unit: mm
    double SegmentLengths[NUM_SEGMENTS] = { 25.72, 16.14, 15.395 , 15.61, 109.40};

    // Mass of each of the actuator sets - unit: kg
    double ActMass[NUM_ACT_SET] = { 8.2859e-06, 8.2859e-06};
    // Array of lambda values for marker locations  (distance from the tip to each of the markers) - unit: mm
    double MarkerLoc[NUM_LOCALIZATION_MARKERS] = { 2.18, 18.79, 50.06, 101.52, 104.02 };
    // Length density for each of the catheter segments
    double rho[NUM_SEGMENTS] = { 2.8814e-7, 2.8814e-7, 2.8814e-7, 2.8814e-7, 2.8814e-7 };
    // unit (1/ mm) Local curvature in unloaded configuration for each of the flexible segments; NUM_FLEX_SEG*3 long array corresponding to NUM_FLEX_SEG many 3x1 vectors
    double ustarlist[NUM_FLEX_SEG][3] = { 0.0000148626102272827, 0.000094448815853795, 0, 0.0007399849479553773, -0.0002292201697658481, 0, 0.0007399849479553773, -0.0002292201697658481, 0 };


    double InsertedLength  =     86;
    // Actuation currents for each of the coils for each of the coil sets - unit: A
//    double ActuationCurrents[NUM_ACT_SET][3] = {{0.000891613642404,-0.015136231914507,0.119640152093819},{0.000949165551101,0.004289529301964,0.197909928730404}};

    double ActuationCurrents[NUM_ACT_SET][3] = {{0.0,-0.0,0.119640152093819},{0.0,0.0,0.19}};

    double damping[NUM_ACT_SET][6] = {{120.1761626666366,120.1761626666366,
                                              284.429938756989,
                                              0.304776127617393,0.304776127617393,
                                              0.00502712804532508}, {120.1761626666366,120.1761626666366,
                                              284.429938756989,
                                              0.304776127617393,0.304776127617393,
                                              0.00502712804532508}};
    double DELTA_T = 0.05;

    double xf_pre[NUM_STATES] , pL[NUM_ACT_SET][3] ,RL[NUM_ACT_SET][9] ;
    double v_L_pre[NUM_ACT_SET][3] = {{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 }};
    double w_L_pre[NUM_ACT_SET][3] = {{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 }};
    // Define initial guesses to be used when solving boundary value problem
    double nL_initialguess[NUM_ACT_SET][3] = {{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 }};
    double mL_initialguess[NUM_ACT_SET][3] = {{ 0.0, 0.0, 0.0 },{ 0.0, 0.0, 0.0 }};


//    double nL_initialguess[NUM_ACT_SET][3] = {{ 2.10841e-05, -0.0131396 ,-0.127776 }, {8.64775e-05, -0.0591842 ,-0.153273}};
//    double mL_initialguess[NUM_ACT_SET][3] = {{ -0.442178, -0.000670959 ,-0.000497768 },{ -0.885441 ,-0.00146746, 0.00025844}};

    /**
     * configuration
     */
	// *** Catheter Configuration in spatial coordinates
	// B0 field vector of the MRI scanner (in spatial coordinates) - unit: Tesla
	double B0[3] = { 0.0, 3.0, 0.0 };
	// Gravity vector - unit: m/s^2
	double gravity[3] = { 0.0, 0.0, 9.81 };
	// Catheter entry point coordinates (in spatial coordinates) - unit: mm
	double p0[3] = { 0.0, 0.0, 0.0 };
	// Catheter entry point orientation (3x3 rotation matrix describing catheter entry point frame orientation relative to the spatial frame stored in row major order) - unit: unitless
	double R0[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

	// *** Other External variables
	// specify if catheter is in free space or if the catheter tip is constrained to a contact point
	ContactModeType ContactMode = ContactModeType::FREE_TIP;
	// External point force (in spatial coordinates) applied at the tip of the catheter (\lambda = 0)  - unit: N
	//   (this will be used when ContactMode == ContactModeType::FREE_TIP)
	double TipForce[3] = { 0.0, 0.0, 0.0 };
	// The spatial coordinates of the point where the catheter tip is constrained to be
	//   (this will be used when ContactMode == ContactModeType::FIXED_TIP)
	double TipConstraintPoint[3] = { 0.0, 0.0, 0.0 };

//	// *** Control Inputs
//	// Inserted Length of the catheter (length of the catheter that is inside the heart chamber) - unit: mm
//    double InsertedLength  = 0.0;
//    for (int i = 0; i < NUM_SEGMENTS; ++i) {
//        InsertedLength += SegmentLengths[i];
//    }

    // *** Numerical Computation Parameters
    // Stepsize used in numerical integration along the length of the catehter during IVP - unit: mm
    double IntegrationStepSize = 0.2;

    /**
     * These are hard coded, need to revise these later
     */
    // we are adding the tubing mass of the coil section to the total mass of actuator unit(kg * mm^2)
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
	CRMCatheterModelParams CathParams{};		// Catheter physical parameters will be packaged into this structure
	CatheterConfiguration CathConfig{};		// Catheter spatial configuration parameters will be packages into this structure


	// Package catheter physical parameters and spatial configuration parameters
	//   this step would typically needs to be executed only once
//	CRMShootingMethodBVP_Prep(B0, gravity, p0, R0,
//		SegmentLengths, MarkerLoc, iRlist, oRlist, YoungModlist, ShearModlist, ustarlist,
//		CoilAlignmentAngles, CoilTurnAreaMat, rho, ActMass,
//		CathParams, CathConfig);

    CRMShootingMethodBVP_Prep(B0,gravity, p0, R0,SegmentLengths, MarkerLoc,
                              iRlist, oRlist, YoungModlist, ShearModlist, ustarlist,
                              CoilAlignmentAngles, CoilTurnAreaMat, rho, ActMass, ActInertia,
                              CathParams, CathConfig);

	// Define initial guesses to be used when solving boundary value problem
    double u0_initialguess[3] = { 0.000709849479553773, -0.0002292201697658481, 0};
	// initial guess for the contstraint force at the catheter tip (this will be used when ContactMode == ContactModeType::FIXED_TIP)
	double ftip_initialguess[3] = { 0.0, 0.0, 0.0 };

	// Declare the output variables for BVP
	// calculated curvature at the catheter base
	double u0_calc[3];
	// calculated contstraint force at the catheter tip (this will be used when ContactMode == ContactModeType::FIXED_TIP)
	double ftip_calc[3];
	// numerical nonlinear equation solver diagnostic outputs
	int localmin;


    /**
     * Initial BVP
     */


    // Define the variable convenience used to package the arguments passed to BVP Solver
	CRMShootingMethodParams<double> BVPParams{};

    // Declare output variables
    // catheter shape state at the entry point of the catheter
    //   states are packed u[0..2], R[0..8], p[0..2](R: 3x3 matrix stored in row major order R11 R12 R13 R21 R22 R23 R31 R32 R33)
    // moment residual at the catheter tip - this should converge to {0,0,0} if the catheter is at its equilibrium configuration
    double residual[NUM_RESIDUAL];
    // spatial coordinates of the localization markers
    double ReportedMarkerPos[NUM_LOCALIZATION_MARKERS][3];

    CRMConstructShootingMethodParamSet<double>(CathParams, CathConfig, InsertedLength, ActuationCurrents, ContactMode,
                                               TipConstraintPoint, TipForce, IntegrationStepSize,
                                               v_L_pre, w_L_pre, pL,  RL, damping, DELTA_T,
                                               BVPParams);

    CRM_FKState(BVPParams, u0_initialguess, ftip_initialguess, xf_pre, pL, RL, ReportedMarkerPos, localmin);
    /**
     * TEST DYN
     */
    for (int i = 0; i < NUM_ACT_SET; ++i) {
        std::cout << "pL: " << pL[i][0] << " " << pL[i][1] << " " << pL[i][2] <<  std::endl;
        std::cout << "RL: " << std::endl;
        std::cout <<  RL[i][0] << " " << RL[i][1] << " " << RL[i][2] <<  std::endl;
        std::cout <<  RL[i][3] << " " << RL[i][4] << " " << RL[i][5] <<  std::endl;
        std::cout <<  RL[i][6] << " " << RL[i][7] << " " << RL[i][8] <<  std::endl;
    }

    std::cout << "xf_pre: " << std::endl;

    for (int i = 0; i < 15; ++i) {
        std::cout << xf_pre[i] << std::endl;
    }
//
//    std::cout << "ActuationCurrents: " << ActuationCurrents[0][0] << " " << ActuationCurrents[0][1] << " " << ActuationCurrents[0][2] <<  std::endl;

    double xf[NUM_STATES];
    double out_u0[3], out_nL[NUM_ACT_SET][3], out_mL[NUM_ACT_SET][3], out_tau[NUM_ACT_SET][3];
    double x_coil[NUM_ACT_SET][NUM_COIL_STATES];

    for (int k = 0; k < 4; ++k) {
        // Populate the variable of convenience for conveniently passing lots of arguments to the BVP Solver
        //   this step would typicall need to be executed every time BVP needs to be calculated as actuation variables would change at every time step
        CRMConstructShootingMethodParamSet<double>(CathParams, CathConfig, InsertedLength, ActuationCurrents, ContactMode,
                                                   TipConstraintPoint, TipForce, IntegrationStepSize,
                                                   v_L_pre, w_L_pre, pL,  RL, damping, DELTA_T,
                                                   BVPParams);


        DynamicsBVP(BVPParams, xf_pre, mL_initialguess, nL_initialguess, ftip_initialguess,
                    out_u0, out_mL, out_nL, out_tau, ftip_calc, localmin);
//
//        std::cout << "out_u0: " << out_u0[0] << " " << out_u0[1] << " " << out_u0[2] <<  std::endl;
//        for (int i = 0; i < NUM_ACT_SET; ++i) {
//            std::cout << "out_mL: " << out_mL[i][0] << " " << out_mL[i][1] << " " << out_mL[i][2] <<  std::endl;
//            std::cout << "out_nL: " << out_nL[i][0] << " " << out_nL[i][1] << " " << out_nL[i][2] <<  std::endl;
//        }

        DYNSolverIVP(BVPParams, out_u0, out_mL, out_nL, out_tau, ftip_calc,
                     true, xf, x_coil,ReportedMarkerPos);

//    double v_L_pre_[NUM_ACT_SET][3], w_L_pre_[NUM_ACT_SET][3], pL_[NUM_ACT_SET][3], RL_[NUM_ACT_SET][9], xf_pre_[NUM_STATES];
        for (int j = 0; j < NUM_ACT_SET; ++j) {
            for (int i = 0; i < 3; ++i) {
                v_L_pre[j][i] = x_coil[j][i];
                w_L_pre[j][i] = x_coil[j][i+3];
                pL[j][i] = x_coil[j][i+6];
            }

            for (int i = 0; i < 9; ++i) {
                RL[j][i] = x_coil[j][i+9];
            }
//            printMatrix(RL[j], 1, 9, "R at Coil");
//            std::cout << " mew coil" << std::endl;

//            for (int i = 0; i < 3; ++i) {
//                mL_initialguess[j][i] = out_mL[j][i];
//                nL_initialguess[j][i] = out_nL[j][i];
//            }

        }

        for (int i = 0; i < NUM_STATES; ++i) {
            xf_pre[i]= xf[i];
        }
        // Print outputs

        std::cout << " *********TEST DYN Configuration********* " << std::endl;

        printMatrix(out_u0, 1, 3, "Calculated curvature at base");
        printMatrix(out_nL[0], 1, 3, "Calculated internal force at L");
        printMatrix(out_mL[0], 1, 3, "Calculated internal moment force at L");

        printMatrix(ftip_calc, 1, 3, "Calculated Tip Force");
        printMatrix(&(xf[12]), 1, 3, "Catheter Tip Position");
        printMatrix(&(xf[3]), 1, 9, "Catheter Tip Orientation");

        for (int i = 0; i < NUM_ACT_SET; ++i) {
            printMatrix(v_L_pre[i], 1, 3, "v at Coil");
            printMatrix(w_L_pre[i], 1, 3, "w at Coil");
            printMatrix(pL[i], 1, 3, "p at Coil");
            printMatrix(RL[i], 1, 9, "R at Coil");
        }

        std::cout << "--localmin--" << localmin << std::endl;

    }

	return (localmin);
}
