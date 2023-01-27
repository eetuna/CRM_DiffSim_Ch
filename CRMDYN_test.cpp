#include <iostream>
#include <cmath>
#include <chrono> 

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

//    const char *path_to_input = "/home/ranhao/Documents/CRM_Dynamics/input_file/circle_01_new.play";
//
//    FILE *file = fopen(path_to_input, "r");
//    if (file == NULL) {
//        printf("Impossible to load input currents! Are you in the right path ?\n");
//    }
//    std::vector< float > input_vec, time_vec;
//    float x;
//    std::ifstream File;
//    File.open(path_to_input);
//    while (File >> x) {
//        input_vec.push_back(x);
//        if (File.peek() == ',')
//            File.ignore();
//    }
//
//    int SIZE_CUR_SET = 10;// input_vec.size() / 9;
//    double input_currents[SIZE_CUR_SET][3], time_step_[SIZE_CUR_SET];
//
//    for (int i = 0; i < SIZE_CUR_SET; ++i) {
//        for (int j = 0; j < 3; ++j) {
//            input_currents[i][j] = 0.001 * input_vec[3 * (3 * i + j) + 1];
//        }
//        time_step_[i] = 0.001 * input_vec[9*(i+1)-1];
//    }

//    for (int i = 0; i < SIZE_CUR_SET; ++i) {
//        std::cout << "input_currents: " << input_currents[i][0] << " " << input_currents[i][1]  << " " << input_currents[i][2] << std::endl;
//        std::cout << time_step_[i] << std::endl;
//    }
//    std::cout << "size_input: " << SIZE_CUR_SET << std::endl;

//    std::vector< std::vector< double > > tip_pts;
//    tip_pts.resize(SIZE_CUR_SET);
//    RunDynamicsExample(SIZE_CUR_SET, input_currents, time_step_, tip_pts);


	RunExample();
//    std::vector< std::vector< double > > tip_pts;
//    RunDynamicsExample(tip_pts);


}

int RunExample(void) {

	std::cout << "### CRM Forward Kinematics Examples... " << std::endl;
	std::cout << std::endl << "Free Space Deflection Example: " << std::endl << std::endl;

    // *** Physical Description of the Catheter
    // IMPORTANT NOTE: For now, most proximal segment is assumed to be always flexible
    //    and the flexible and rigid segments are assumed to be alternating
    //    most distal segment can be flexible or rigid
    // Items are listed in distal-to-proximal order (starting from the tip of the catheter towards the base)
    // Outer radii of each of the flexible segments - unit: mm
    double oRlist[NUM_FLEX_SEG] = { 1.5875, 1.5875 };
    // Inner radii of each of the flexible segments - unit: mm
    double iRlist[NUM_FLEX_SEG] = { 0.9906, 0.9906 };
    // Young's moduli of each of the flexible segments - unit: ??
    double YoungModlist[NUM_FLEX_SEG] =  { 5.3948, 5.3948 };// { 4.9105, 4.9105 };
    // Shear moduli of each of the flexbile segments - unit: ??
    double ShearModlist[NUM_FLEX_SEG] = { 2.3881, 2.3881 };//{ 1.6545, 1.6545 };
    // Alignment angles for the side coils, for each of the actuator sets - unit: radians
    double CoilAlignmentAngles[NUM_ACT_SET][2] = { {0.0, 0.0}  };// { {-0.0871, -0.3934}  };
    // Coil turn area matrices for each of the actuator sets  - unit: mm2
    double CoilTurnAreaMat[NUM_ACT_SET][9] = { { 1.3851,  0.0, 0.0, 0.0, 1.44, 0.0, 0.0, 0.0, 1.60 } };
    // Lengths of each of the catheter segments - unit: mm
    double SegmentLengths[NUM_SEGMENTS] = { 19.85,  18.3, 59.80 };

    // Mass of each of the actuator sets - unit: ??
    double ActMass[NUM_ACT_SET] = { 5.7736e-5};
    // Array of lambda values for marker locations  (distance from the tip to each of the markers) - unit: mm
    double MarkerLoc[NUM_LOCALIZATION_MARKERS] = { 2.18, 18.79, 50.06, 101.52, 104.02 };
    // Length density for each of the catheter segments
    double rho[NUM_SEGMENTS] = { 2.8814e-7, 2.8814e-7, 2.8814e-7 };
    // Local curvature in unloaded configuration for each of the flexible segments; NUM_FLEX_SEG*3 long array corresponding to NUM_FLEX_SEG many 3x1 vectors
    double ustarlist[NUM_FLEX_SEG][3] = { 0.0000148626102272827, 0.000094448815853795, 0, 0.0007399849479553773, -0.0002292201697658481, 0};//0.00567, -0.00822, 0.0};


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

	// *** Control Inputs
	// Inserted Length of the catheter (length of the catheter that is inside the heart chamber) - unit: mm
    double InsertedLength  = 0.0;
    for (int i = 0; i < NUM_SEGMENTS; ++i) {
        InsertedLength += SegmentLengths[i];
    }

    // *** Numerical Computation Parameters
    // Stepsize used in numerical integration along the length of the catehter during IVP - unit: mm
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
    double v_L_pre[NUM_ACT_SET][3] = { 0.0, 0.0, 0.0 };
    double w_L_pre[NUM_ACT_SET][3] = { 0.0, 0.0, 0.0 };

	// Define variables of convenience that will be use to package catheter physical parameters and spatial configuration parameters
	CRMCatheterModelParams CathParams{};		// Catheter physical parameters will be packaged into this structure
	CatheterConfiguration CathConfig{};		// Catheter spatial configuration parameters will be packages into this structure


	// Package catheter physical parameters and spatial configuration parameters
	//   this step would typically needs to be executed only once

    CRMShootingMethodBVP_Prep(B0,gravity, p0, R0,SegmentLengths, MarkerLoc,
                              iRlist, oRlist, YoungModlist, ShearModlist, ustarlist,
                              CoilAlignmentAngles, CoilTurnAreaMat, rho, ActMass, ActInertia,
                              CathParams, CathConfig);

	// Define initial guesses to be used when solving boundary value problem
    double u0_initialguess[NUM_FLEX_SEG][3] = { 0.000709849479553773, -0.0002292201697658481, 0, 0.0000148626102272827, 0.000094448815853795, 0};
	// initial guess for the contstraint force at the catheter tip (this will be used when ContactMode == ContactModeType::FIXED_TIP)
	double ftip_initialguess[3] = { 0.0, 0.0, 0.0 };
    // Define initial guesses to be used when solving boundary value problem
    double nL_initialguess[NUM_ACT_SET][3] = { 0.0, 0.0, 0.0 };

    double mL_initialguess[NUM_ACT_SET][3] = { 0.0, 0.0, 0.0 };
	// Declare the output variables for BVP
	// calculated curvature at the catheter base
	double u0_calc[NUM_FLEX_SEG][3];
	// calculated contstraint force at the catheter tip (this will be used when ContactMode == ContactModeType::FIXED_TIP)
	double ftip_calc[3];
	// numerical nonlinear equation solver diagnostic outputs
	int localmin;


    /**
     * Initial BVP
     */
    // Actuation currents for each of the coils for each of the coil sets - unit: A
    double ActuationCurrents[NUM_ACT_SET][3] = {0.0, 0.0, 0.0};// { {0.100, 0.100, 0.100}};
//    InsertedLength = 98.5;

    double damping[NUM_ACT_SET][6] = {0, 0, 0, 0, 0, 0};
    double DELTA_T = 0.01;

    // Define the variable convenience used to package the arguments passed to BVP Solver
	CRMShootingMethodParams<double> BVPParams{};

    // Declare output variables
    // catheter shape state at the entry point of the catheter
    //   states are packed u[0..2], R[0..8], p[0..2](R: 3x3 matrix stored in row major order R11 R12 R13 R21 R22 R23 R31 R32 R33)
    // moment residual at the catheter tip - this should converge to {0,0,0} if the catheter is at its equilibrium configuration
    double residual[NUM_RESIDUAL];
    // spatial coordinates of the localization markers
    double ReportedMarkerPos[NUM_LOCALIZATION_MARKERS][3];
    double xf[NUM_STATES];

    // Catheter entry point coordinates (in spatial coordinates) - unit: mm
    double pL[NUM_ACT_SET][3] = { 0.0, 0.0, 0.0 };
    // Catheter entry point orientation (3x3 rotation matrix describing catheter entry point frame orientation relative to the spatial frame stored in row major order) - unit: unitless
    double RL[NUM_ACT_SET][9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
//     Populate the variable of convenience for conveniently passing lots of arguments to the BVP Solver
//	   this step would typicall need to be executed every time BVP needs to be calculated as actuation variables would change at every time step
	CRMConstructShootingMethodParamSet<double>(CathParams, CathConfig, InsertedLength, ActuationCurrents, ContactMode,
                                               TipConstraintPoint, TipForce, IntegrationStepSize,
                                               v_L_pre, w_L_pre, pL,  RL, damping, DELTA_T,
                                               BVPParams);


    // Cosserat Rod Model - Solve the Boundary Value Problem to calculate the equilibrium configuration of the catheter
	CRMShootingMethodBVP(BVPParams, u0_initialguess, u0_calc, ftip_calc, localmin);


//
//    // output for time advance, not used in BVP, just placeholders

//
//    // Cosserat Rod Model - Solve the Initial Value Problem to calculate the shape of the catheter
    CRMSolverIVP(BVPParams, u0_calc, ftip_calc, false, xf, residual, pL, RL, ReportedMarkerPos);

//    std::cout << " ----Initial Configuration under NO actuation------ " << std::endl;
//    printMatrix(u0_calc, 1, 3, "Calculated curvature at base");
//
//    printMatrix(ftip_calc, 1, 3, "Calculated Tip Force");
//    printMatrix(&(xf[12]), 1, 3, "Catheter Tip Position");
//    printMatrix(&(xf[3]), 1, 9, "Catheter Tip Orientation");

//    double pL[3] = {-0.711731, -1.6718, 69.4687};
//    double RL[9] = {0.999836, -0.00038574, -0.0181181,
//                    -0.00038574 ,0.999094 ,-0.042558,
//    0.0181181, 0.042558, 0.99893};

//    for (int i = 0; i < 3; ++i) {
//        pL[i] = pL_[i];
//    }
//    for (int i = 0; i < 9; ++i) {
//        RL[i] = RL_[i];
//    }


    std::cout << "pL: " << pL[0][0] << " " << pL[0][1] << " " << pL[0][2] <<  std::endl;
    std::cout << "RL: " << std::endl;
    std::cout <<  RL[0][0] << " " << RL[0][1] << " " << RL[0][2] <<  std::endl;
    std::cout <<  RL[0][3] << " " << RL[0][4] << " " << RL[0][5] <<  std::endl;
    std::cout <<  RL[0][6] << " " << RL[0][7] << " " << RL[0][8] <<  std::endl;


    /**
     * TEST DYN
     */
    // Actuation currents for each of the coils for each of the coil sets - unit: A

    ActuationCurrents[0][2] = 0.02;

    std::cout << "ActuationCurrents: " << ActuationCurrents[0][0] << " " << ActuationCurrents[0][1] << " " << ActuationCurrents[0][2] <<  std::endl;

    // Populate the variable of convenience for conveniently passing lots of arguments to the BVP Solver
    //   this step would typicall need to be executed every time BVP needs to be calculated as actuation variables would change at every time step
    CRMConstructShootingMethodParamSet<double>(CathParams, CathConfig, InsertedLength, ActuationCurrents, ContactMode,
                                               TipConstraintPoint, TipForce, IntegrationStepSize,
                                               v_L_pre, w_L_pre, pL,  RL, damping, DELTA_T,
                                               BVPParams);

    double out_u0[NUM_FLEX_SEG][3], out_nL[NUM_ACT_SET][3], out_mL[NUM_ACT_SET][3];
    DynamicsBVP(BVPParams, u0_initialguess, mL_initialguess, nL_initialguess, ftip_initialguess,
                   out_u0, out_mL, out_nL, ftip_calc, localmin);


    std::cout << "out_u0 0: " << out_u0[0][0] << " " << out_u0[0][1] << " " << out_u0[0][2] <<  std::endl;
    std::cout << "out_u0 1: " << out_u0[1][0] << " " << out_u0[1][1] << " " << out_u0[1][2] <<  std::endl;

//    std::cout << "out_mL: " << out_mL[0] << " " << out_mL[1] << " " << out_mL[2] <<  std::endl;
//    std::cout << "out_nL: " << out_nL[0] << " " << out_nL[1] << " " << out_nL[2] <<  std::endl;



    double x_coil[NUM_COIL_STATES];
    DYNSolverIVP(BVPParams, out_u0, out_mL, out_nL, ftip_calc,
                    true, xf, x_coil,ReportedMarkerPos);

    double v_L_update[3], w_L_update[3], pL_update[3], RL_update[9];
    for (int i = 0; i < 3; ++i) {
        v_L_update[i] = x_coil[i];
        w_L_update[i] = x_coil[i+3];
        pL_update[i] = x_coil[i+6];
    }

    for (int i = 0; i < 9; ++i) {
        RL_update[i] = x_coil[i+9];
    }
	// Print outputs
    std::cout << " *********TEST DYN Configuration********* " << std::endl;
//    printMatrix(&(out_u0[0]), 1, 3, "Calculated curvature at base");
//    printMatrix(&(out_nL[0]), 1, 3, "Calculated internal force at L");
//    printMatrix(&(out_mL[0]), 1, 3, "Calculated internal moment force at L");

    printMatrix(ftip_calc, 1, 3, "Calculated Tip Force");
	printMatrix(&(xf[12]), 1, 3, "Catheter Tip Position");
	printMatrix(&(xf[3]), 1, 9, "Catheter Tip Orientation");

//    printMatrix(&(xf[15]), 1, 3, "v at tip");
//    printMatrix(&(xf[18]), 1, 3, "w at tip");
//
    printMatrix(v_L_update, 1, 3, "v at Coil");
    printMatrix(w_L_update, 1, 3, "w at Coil");
    printMatrix(pL_update, 1, 3, "p at Coil");
    printMatrix(RL_update, 1, 9, "R at Coil");


//    printMatrix(residual, 1, 3, "residual at coil");

    for (int i = 0; i < NUM_LOCALIZATION_MARKERS; ++i) {
        std::cout << "ReportedMarkerPos: " << ReportedMarkerPos[i][0] << " " << ReportedMarkerPos[i][1] << " " << ReportedMarkerPos[i][2] << std::endl;
    }

    std::cout << "----" << std::endl;

	return (localmin);
}
