#include <iostream>
#include <cmath>
#include <chrono> 

#include "CRMTest.h"
#include "CRMDYN.hpp"
using namespace std::chrono;


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
int RunFKExample(void);

int main(int argc, char** argv) {

	RunExample();
//
//	RunFKExample();
//
//	std::cout << std::endl << std::endl << "### Running Tests ..." << std::endl << std::endl;
//	int FailFlag = RunTests();
//	if (FailFlag == 0) {
//		std::cout << std::endl << std::endl << "###" << std::endl;
//		std::cout << "ALL TESTS PASSED!..." << std::endl;
//		std::cout << "###" << std::endl;
//	}
//	else {
//		std::cout << std::endl << std::endl << "###" << std::endl;
//		std::cout << FailFlag << " TEST(S) FAILED!..." << std::endl;
//		std::cout << "###" << std::endl;
//	}
//	return (FailFlag);



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
	double YoungModlist[NUM_FLEX_SEG] = { 5.3948, 5.3948 };
	// Shear moduli of each of the flexbile segments - unit: ??
	double ShearModlist[NUM_FLEX_SEG] = { 2.3881, 2.3881 };
	// Alignment angles for the side coils, for each of the actuator sets - unit: radians
	double CoilAlignmentAngles[NUM_ACT_SET][2] = { {0.0, 0.0}  };
	// Coil turn area matrices for each of the actuator sets  - unit: mm2
	double CoilTurnAreaMat[NUM_ACT_SET][9] = { { 1.44, 0.0, 0.0, 0.0, 1.3851, 0.0, 0.0, 0.0, 1.60 } };
	// Lengths of each of the catheter segments - unit: mm
	double SegmentLengths[NUM_SEGMENTS] = { 10.72, 16.14, 15.395 };

    double SegEnds[NUM_SEGMENTS] = { 10.72, 26.86, 42.255 };

    // Mass of each of the actuator sets - unit: ??
	double ActMass[NUM_ACT_SET] = { 7.7736e-5};
	// Array of lambda values for marker locations  (distance from the tip to each of the markers) - unit: mm
	double MarkerLoc[NUM_LOCALIZATION_MARKERS] = { 2.18, 18.79, 50.06, 101.52, 104.02 };
	// Length density for each of the catheter segments
	double rho[NUM_SEGMENTS] = { 2.8814e-7, 2.8814e-7, 2.8814e-7 };
	// Local curvature in unloaded configuration for each of the flexible segments; NUM_FLEX_SEG*3 long array corresponding to NUM_FLEX_SEG many 3x1 vectors
	double ustarlist[NUM_FLEX_SEG][3] = { 0.000128, -0.000173, 0.0, 0.00567, -0.00822, 0.0};

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
	double InsertedLength = 104.0;
	// Actuation currents for each of the coils for each of the coil sets - unit: A
	double ActuationCurrents[NUM_ACT_SET][3] = { {0.100, 0.100, 0.100}};

	// *** Numerical Computation Parameters
	// Stepsize used in numerical integration along the length of the catehter during IVP - unit: mm
	double IntegrationStepSize = 0.2;

    double ActInertia[NUM_ACT_SET][9];
    for (int i = 0; i < NUM_ACT_SET; ++i)
    {
        double I_zz = 0.5 * ActMass[i] * (oRlist[0] * oRlist[0] + iRlist[0] * iRlist[0]);
        double I_xx = 0.25 * ActMass[i] * (oRlist[0] * oRlist[0] + iRlist[0] * iRlist[0]) + 1 / 12 * ActMass[i] * SegEnds[2*i+1];
        ActInertia[i][0] = I_xx; ActInertia[i][1] = 0.0; ActInertia[i][2] = 0.0;
        ActInertia[i][3] = 0.0; ActInertia[i][4] = I_xx; ActInertia[i][5] = 0.0;
        ActInertia[i][6] = 0.0; ActInertia[i][7] = 0.0; ActInertia[i][8] = I_zz;
    }
    double v_L_pre[3] = { 0.0, 0.0, 0.0 };
    double w_L_pre[3] = { 0.0, 0.0, 0.0 };

    double v0[3] = {0.0,0.0,0.0};              //initial linear velocity in local frame of entry point
    double w0[3] = {0.0,0.0,0.0};

	// Define variables of convenience that will be use to package catheter physical parameters and spatial configuration parameters
	CRMCatheterModelParams CathParams;		// Catheter physical parameters will be packaged into this structure
	CatheterConfiguration CathConfig;		// Catheter spatial configuration parameters will be packages into this structure


	// Package catheter physical parameters and spatial configuration parameters
	//   this step would typically needs to be executed only once
//	CRMShootingMethodBVP_Prep(B0, gravity, p0, R0,
//		SegmentLengths, MarkerLoc, iRlist, oRlist, YoungModlist, ShearModlist, ustarlist,
//		CoilAlignmentAngles, CoilTurnAreaMat, rho, ActMass,
//		CathParams, CathConfig);

    CRMShootingMethodBVP_Prep(B0,gravity, p0, R0, v0, w0,
                              SegmentLengths, MarkerLoc, iRlist, oRlist, YoungModlist, ShearModlist, ustarlist,
                              CoilAlignmentAngles, CoilTurnAreaMat, rho, ActMass, ActInertia,
                              CathParams, CathConfig);

	// Define initial guesses to be used when solving boundary value problem
	// initial guess for the curvature at the catheter base
	double u0_initialguess[3] = { 0.0, 0.0, 0.0 };
	// initial guess for the contstraint force at the catheter tip (this will be used when ContactMode == ContactModeType::FIXED_TIP)
	double ftip_initialguess[3] = { 0.0, 0.0, 0.0 };

	// Declare the output variables for BVP
	// calculated curvature at the catheter base
	double u0_calc[3];
	// calculated contstraint force at the catheter tip (this will be used when ContactMode == ContactModeType::FIXED_TIP)
	double ftip_calc[3];
	// numerical nonlinear equation solver diagnostic outputs
	int localmin;

	// Define the variable convenience used to package the arguments passed to BVP Solver
	CRMShootingMethodParams<double> BVPParams;
	// Populate the variable of convenience for conveniently passing lots of arguments to the BVP Solver
	//   this step would typicall need to be executed every time BVP needs to be calculated as actuation variables would change at every time step
	CRMConstructShootingMethodParamSet<double>(CathParams, CathConfig, InsertedLength, ActuationCurrents, ContactMode, TipConstraintPoint, TipForce, IntegrationStepSize, v_L_pre, w_L_pre, BVPParams);

    // Cosserat Rod Model - Solve the Boundary Value Problem to calculate the equilibrium configuration of the catheter
	CRMShootingMethodBVP(BVPParams, u0_initialguess, ftip_initialguess, u0_calc, ftip_calc, localmin);

	// Declare input variables for IVP
	// catheter shape state at the entry point of the catheter
	//   states are packed u[0..2], R[0..8], p[0..2](R: 3x3 matrix stored in row major order R11 R12 R13 R21 R22 R23 R31 R32 R33)
	double x0[NUM_STATES];
	for (int i = 0; i < 3; i++) x0[i] = u0_calc[i];
	for (int i = 0; i < 9; i++) x0[i + 3] = R0[i];
	for (int i = 0; i < 3; i++) x0[i + 3 + 9] = p0[i];

	// Declare output variables
	// catheter shape state at the entry point of the catheter
	//   states are packed u[0..2], R[0..8], p[0..2](R: 3x3 matrix stored in row major order R11 R12 R13 R21 R22 R23 R31 R32 R33)
	double xf[NUM_STATES];
	// moment residual at the catheter tip - this should converge to {0,0,0} if the catheter is at its equilibrium configuration
	double residual[3];
	// spatial coordinates of the localization markers
	double ReportedMarkerPos[NUM_LOCALIZATION_MARKERS][3];

	// Cosserat Rod Model - Solve the Initial Value Problem to calculate the shape of the catheter 
	CRMSolverIVP(BVPParams, x0, ftip_calc, false, xf, residual, ReportedMarkerPos);

	// Print outputs
	std::cout << "localmin  " << localmin << std::endl;
	printMatrix(u0_calc, 1, EQNDIMENSION, "Calculated curvature at base");

	printMatrix(ftip_calc, 1, 3, "Calculated Tip Force");
	printMatrix(&(xf[12]), 1, 3, "Catheter Tip Position");
	printMatrix(&(xf[3]), 1, 9, "Catheter Tip Orientation");

    printMatrix(&(xf[15]), 1, 3, "v");
    printMatrix(&(xf[18]), 1, 3, "w");

    printMatrix(residual, 1, 3, "Force residual at coil");
    printMatrix(&(residual[3]), 1, 3, "Moment residual at coil");

    std::cout << "----" << std::endl;

	return (localmin);
}

//
//int RunFKExample(void) {
//
//	std::cout << "### CRM Forward Kinematics Examples... " << std::endl;
//	std::cout << std::endl << "Free Space Deflection Example: " << std::endl << std::endl;
//
//	// *** Physical Description of the Catheter
//	// IMPORTANT NOTE: For now, most proximal segment is assumed to be always flexible
//	//    and the flexible and rigid segments are assumed to be alternating
//	//    most distal segment can be flexible or rigid
//	// Items are listed in distal-to-proximal order (starting from the tip of the catheter towards the base)
//	// Outer radii of each of the flexible segments - unit: mm
//	double oRlist[NUM_FLEX_SEG] = { 1.5875, 1.5875, 1.5875 };
//	// Inner radii of each of the flexible segments - unit: mm
//	double iRlist[NUM_FLEX_SEG] = { 0.9906, 0.9906, 0.9906 };
//	// Young's moduli of each of the flexible segments - unit: ??
//	double YoungModlist[NUM_FLEX_SEG] = { 5.3948, 5.3948, 5.3948 };
//	// Shear moduli of each of the flexbile segments - unit: ??
//	double ShearModlist[NUM_FLEX_SEG] = { 2.3881, 2.3881, 2.3881 };
//	// Alignment angles for the side coils, for each of the actuator sets - unit: radians
//	double CoilAlignmentAngles[NUM_ACT_SET][2] = { {0.0, 0.0} , {0.0, 0.0} };
//	// Coil turn area matrices for each of the actuator sets  - unit: mm2
//	double CoilTurnAreaMat[NUM_ACT_SET][9] = { { 1.44, 0.0, 0.0, 0.0, 1.3851, 0.0, 0.0, 0.0, 1.60 } , { 1.65, 0.0, 0.0, 0.0, 1.326, 0.0, 0.0, 0.0, 1.56 } };
//	// Lengths of each of the catheter segments - unit: mm
//	double SegmentLengths[NUM_SEGMENTS] = { 10.72, 16.14, 15.395, 15.61, 60.00 };
//	// Mass of each of the actuator sets - unit: ??
//	double ActMass[NUM_ACT_SET] = { 7.7736e-5, 8.0157e-5 };
//	// Array of lambda values for marker locations  (distance from the tip to each of the markers) - unit: mm
//	double MarkerLoc[NUM_LOCALIZATION_MARKERS] = { 2.18, 18.79, 50.06, 101.52, 104.02 };
//	// Length density for each of the catheter segments
//	double rho[NUM_SEGMENTS] = { 2.8814e-7, 2.8814e-7, 2.8814e-7, 2.8814e-7, 2.8814e-7 };
//	// Local curvature in unloaded configuration for each of the flexible segments; NUM_FLEX_SEG*3 long array corresponding to NUM_FLEX_SEG many 3x1 vectors
//	double ustarlist[NUM_FLEX_SEG][3] = { 0.000128, -0.000173, 0.0, 0.00567, -0.00822, 0.0, -0.00250, 0.00188, 0.0 };
//
//	// *** Catheter Configuration in spatial coordinates
//	// B0 field vector of the MRI scanner (in spatial coordinates) - unit: Tesla
//	double B0[3] = { 0.0, 3.0, 0.0 };
//	// Gravity vector - unit: ??
//	double gravity[3] = { 0.0, 0.0, 9.81 };
//	// Catheter entry point coordinates (in spatial coordinates) - unit: mm
//	double p0[3] = { 0.0, 0.0, 0.0 };
//	// Catheter entry point orientation (3x3 rotation matrix describing catheter entry point frame orientation relative to the spatial frame stored in row major order) - unit: unitless
//	double R0[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
//
//	// *** Other External variables
//	// specify if catheter is in free space or if the catheter tip is constrained to a contact point
//	ContactModeType ContactMode = ContactModeType::FREE_TIP;
//	// External point force (in spatial coordinates) applied at the tip of the catheter (\lambda = 0)  - unit: ??
//	//   (this will be used when ContactMode == ContactModeType::FREE_TIP)
//	double TipForce[3] = { 0.0, 0.0, 0.0 };
//	// The spatial coordinates of the point where the catheter tip is constrained to be
//	//   (this will be used when ContactMode == ContactModeType::FIXED_TIP)
//	double TipConstraintPoint[3] = { 0.0, 0.0, 0.0 };
//
//	// *** Control Inputs
//	// Inserted Length of the catheter (length of the catheter that is inside the heart chamber) - unit: mm
//	double InsertedLength = 104.0;
//	// Actuation currents for each of the coils for each of the coil sets - unit: A
//	double ActuationCurrents[NUM_ACT_SET][3] = { {0.100, 0.100, 0.100}, {0.100, 0.100, 0.100} };
//
//	// *** Numerical Computation Parameters
//	// Stepsize used in numerical integration along the length of the catehter during IVP - unit: mm
//	double IntegrationStepSize = 0.2;
//
//	// Define variables of convenience that will be use to package catheter physical parameters and spatial configuration parameters
//	CRMCatheterModelParams CathParams;		// Catheter physical parameters will be packaged into this structure
//	CatheterConfiguration CathConfig;		// Catheter spatial configuration parameters will be packages into this structure
//
//	// Package catheter physical parameters and spatial configuration parameters
//	//   this step would typically needs to be executed only once
//	CRMShootingMethodBVP_Prep(B0, gravity, p0, R0,
//		SegmentLengths, MarkerLoc, iRlist, oRlist, YoungModlist, ShearModlist, ustarlist,
//		CoilAlignmentAngles, CoilTurnAreaMat, rho, ActMass,
//		CathParams, CathConfig);
//
//	// Define initial guesses to be used when solving boundary value problem
//	// initial guess for the curvature at the catheter base
//	double u0_initialguess[3] = { 0.0, 0.0, 0.0 };
//	// initial guess for the contstraint force at the catheter tip (this will be used when ContactMode == ContactModeType::FIXED_TIP)
//	double ftip_initialguess[3] = { 0.0, 0.0, 0.0 };
//
//	// Declare the output variables for BVP
//	// calculated curvature at the catheter base
//	double u0_calc[3];
//	// calculated contstraint force at the catheter tip (this will be used when ContactMode == ContactModeType::FIXED_TIP)
//	double ftip_calc[3];
//	// numerical nonlinear equation solver diagnostic outputs
//	int localmin;
//
//	//
//	// New Forward Kinematics Function
//	//
//
//	double ReportedMarkerPos[NUM_LOCALIZATION_MARKERS][3];  // temporary storage for storing localization marker positions
//		// The Forward Kinematics Function uses a different parameter structure than BVP
//	CRMForwardKinematicsData FKParams;
//	FKParams.CathConfig = &CathConfig;
//	FKParams.CathParams = &CathParams;
//	FKParams.ContactMode = ContactMode;
//	FKParams.FinalValueOnly = false;  // we want the localization coil locations, too
//	FKParams.ReportedMarkerPos = &ReportedMarkerPos;
//	int outputdim;
//	if (ContactMode == ContactModeType::FREE_TIP) outputdim = 15; else outputdim = 18;
//	for (int i = 0; i < 3; i++) FKParams.TipConstraintPoint[i] = TipConstraintPoint[i];
//	for (int i = 0; i < 3; i++) FKParams.TipForce[i] = TipForce[i];
//	for (int i = 0; i < 3; i++) FKParams.u0_initialguess[i] = u0_initialguess[i];
//	for (int i = 0; i < 3; i++) FKParams.ftip_initialguess[i] = ftip_initialguess[i];
//	FKParams.IntegrationStepSize = IntegrationStepSize;
//	// and inputs
//	double control_inputs[NUM_ACT_SET * 3 + 1];  // actuator currents (distal to proximal) followed by inserted length
//	for (int i = 0; i < NUM_ACT_SET; i++)  for (int j = 0; j < 3; j++) control_inputs[i*3+j] = ActuationCurrents[i][j];
//	control_inputs[NUM_ACT_SET * 3] = InsertedLength;
//	double * output_values = new double [outputdim]; // 18 for FIXED_TIP // p[0..2],R[0..8],u0[0..2](,ftip[0..2])  R: in row major order
//
//	// Multiple repetitions to more reliably measure time
//	int REPS = 100;
//	// Get starting timepoint
//	auto start = high_resolution_clock::now();
//
//	// Cosserat Rod Model - Solve the Forward Kinematics
//	for (int cnt = 0; cnt < REPS; cnt++)
//		CRM_ForwardKinematics<double>(control_inputs, output_values, FKParams);
//
//	// Get ending timepoint
//	auto stop = high_resolution_clock::now();
//	// Get duration. Substart timepoints to
//	// get duration. To cast it to proper unit
//	// use duration cast method
//	auto duration = duration_cast<microseconds>(stop - start);
//	std::cout << std::endl << "Average time taken by BVP Solution in " << REPS << " repetitions: " << duration.count() / REPS << " microseconds" << std::endl;
//
//	printMatrix(output_values, 1, outputdim, "FK Output -- p,R,u0");
//
//	delete[] output_values;
//	return 0;
//}
//
//
//int RunTests (void) {
//
//	int FailFlag = 0;
//
//	std::cout << "######## CRM Tools Test Results " << std::endl;
//
//	//
//	//  Test code for CRMSolverIVP  (Solver for Initial Value Problem)
//	//
//	std::cout << "### CRMSolverIVP Test: " << std::endl;
//
//	// set up the test problem inputs
//	double IntegrationStepSize = 0.2;
//	double InsertedLength = 104.02;
//	double B0[3] = { 0.0, 3.0, 0.0 };			// B0 field vector of the MRI scanner (in spatial coordinates)
//	double x0[NUM_STATES] = { 0.0,0.0,0.0,   1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,   0.0, 0.0, 0.0 };
//	bool FinalStepOnly = false;
//	const int Na = NUM_ACT_SET; 				// Number of actuator sets
//	const int Nc = NUM_FLEX_SEG;				// Number of flexible catheter segments
//	const int Ns = NUM_SEGMENTS;  				// NumSegments: Ns = 2 * Na + 1   // derived quantity
//	const int Nm = NUM_LOCALIZATION_MARKERS;	// Number of Markers
//	//	For all parameters below, segments and actuator units are numbered/ordered from the tip of the catheter towards the base
//	// Array of lambda values for segment endpoints; Ns long array
//	double SegEnds[Ns] = { 10.72, 26.86, 42.255, 57.865, 104.02 };
//	// Array of lambda values for marker locations
//	double MarkerLoc[Nm] = { 2.18, 18.79, 50.06, 101.52, 104.02 };
//	// Catheter Rigidity Matrix; (Na+1)*9 long array, (Na+1) 3x3 matrices stored in row major order
//	double Klist[Nc][9] = { {22.8304161859862, 0.0, 0.0, 0.0, 22.8304161859862, 0.0, 0.0, 0.0, 20.2125442625320}, {22.8304161859862, 0.0, 0.0, 0.0, 22.8304161859862, 0.0, 0.0, 0.0, 20.2125442625320}, {22.8304161859862, 0.0, 0.0, 0.0, 22.8304161859862, 0.0, 0.0, 0.0, 20.2125442625320} };
//	// Inverses of K matrices; (Na+1)*9 long array, (Na+1) 3x3 matrices stored in row major order  // derived quantity
//	double Kinvlist[Nc][9] = { {0.0438012164059375, 0.0, 0.0, 0.0, 0.0438012164059375, 0.0, 0.0, 0.0, 0.0494742268470230}, {0.0438012164059375, 0.0, 0.0, 0.0, 0.0438012164059375, 0.0, 0.0, 0.0, 0.0494742268470230}, {0.0438012164059375, 0.0, 0.0, 0.0, 0.0438012164059375, 0.0, 0.0, 0.0, 0.0494742268470230} };
//	// Local curvature in unloaded configuration for each of the flexible segments; (Na+1)*3 long array, (Na+1) 3x1 vectors
//	double ustarlist[Nc][3] = { 0.000127537459505030, -0.000173238845806410, 0.0, 0.00566574659169662, -0.00822266255995812, 0.0, -0.00250030480901699, 0.00188434148765547, 0.0 }; //{-0.0025, 0.0019, 0.0, 0.0057, -0.0082, 0.0, 0.00013, -0.00017, 0.0};
//	// Actuator magnetization moments; Na*3 long array, Na 3x1 vectors; MagMoment = CoilAlignMat * CoilTurnAreaMat * ActuationCurrentVector
//	double MagMomentlist[Na][3] = { 0.0, 0.0, 0.1600, 0.0, 0.0, 0.0 };
//	// Reciprocal of \Delta s (= \Delta \lambda) used in discretizing fcum  ( deltasinv = 1 / (L/N) = N/L)
//	double dlambdainv = 1.0;
//	double FullLength = SegEnds[Ns - 1];
//	// 3*(N+1) by 1 array (grouped by 3 doubles) storing cumulative external force integrated from \lambda = index * \Delta\lambda to the catheter tip (\lambda=0)
//	double fcumIVP[NUM_FCUM_LAMBDA + 1][3] = {};
//	double ftip[3] = { 0.0, 0.0, 0.0 };
//
//	// declare the output variables
//	double xf[NUM_STATES];
//	double residual[3];
//	double ReportedMarkerPos[NUM_LOCALIZATION_MARKERS][3];
//
//	// Cosserat Rod Model - Integrator for Solving the Initial Value Problem
//	CRMSolverIVP<double>(x0, IntegrationStepSize, InsertedLength, dlambdainv, SegEnds, MarkerLoc, Klist, Kinvlist, ustarlist, MagMomentlist, fcumIVP, ftip, B0, FinalStepOnly, xf, residual, ReportedMarkerPos);
//
//	// compare the results to the reference values
//	double ExpectedStates[NUM_STATES] = { 0.0225, -0.0024, -0.0074, 0.982564773234176, 0.029114762388711, -0.183626784558004, -0.094107711560958, 0.929675871807347, -0.356155178150316, 0.160344017817671, 0.367226228689942, 0.916206686783898, -5.497054669396932, -5.570783408417309, 103.0917629459139 }; //  ???f, ???f, ???f,  ???f, ???f, ???f, ???f, ???f, ???f, ???f, ???f, ???f,   ???f, ???f, ???f };
//	double ExpectedMarkerPos[NUM_LOCALIZATION_MARKERS][3] = { {-5.1038, -4.8439, 101.0747}, {-2.4420, -1.9709, 84.9666}, {0.0000, 0.0000, 53.9600}, {0.0000, 0.0000, 2.5000}, {0, 0, 0} };
//	double errorv2[NUM_STATES], ssqerrorv2;
//	double errorv3[NUM_LOCALIZATION_MARKERS][3], ssqerrorv3;
//
//	ssqerrorv2 = 0.0;
//	for (int i = 0; i < NUM_STATES; i++) ssqerrorv2 += pow((errorv2[i] = ExpectedStates[i] - xf[i]), 2);
//	ssqerrorv3 = 0.0;
//	for (int i = 0; i < NUM_LOCALIZATION_MARKERS; i++)
//		for (int j = 0; j < 3; j++)
//			ssqerrorv3 += pow((errorv3[i][j] = ExpectedMarkerPos[i][j] - ReportedMarkerPos[i][j]), 2);
//
//
//	printMatrix(ExpectedStates, 1, NUM_STATES, "ExpectedStates");
//	printMatrix(xf, 1, NUM_STATES, "xf");
//	printMatrix(errorv2, 1, NUM_STATES, "errorv2");
//	std::cout << "Error Norm:" << sqrt(ssqerrorv2) << std::endl;
//	std::cout << "----" << std::endl;
//	printMatrix(&(ExpectedMarkerPos[0][0]), NUM_LOCALIZATION_MARKERS, 3, "ExpectedMarkerPos");
//	printMatrix(&(ReportedMarkerPos[0][0]), NUM_LOCALIZATION_MARKERS, 3, "ReportedMarkerPos");
//	printMatrix(&(errorv3[0][0]), NUM_LOCALIZATION_MARKERS, 3, "errorv3");
//	std::cout << "Error Norm:" << sqrt(ssqerrorv3) << std::endl;
//	std::cout << "----" << std::endl;
//
//	if ((sqrt(ssqerrorv2) < 0.2) && (sqrt(ssqerrorv3) < 0.2)) {
//		std::cout << "Test passed for Inserted Length =" << InsertedLength << std::endl;
//		std::cout << "###" << std::endl;
//	}
//	else {
//		std::cout << "TEST FAIL: " << "Results do not match!" << std::endl;
//		std::cout << "###" << std::endl;
//		FailFlag++;
//	}
//
//
//	//
//	//  Test code for CRMShootingMethodBVP (Shooting method for solving the boundary value problem)
//	//
//
//	std::cout << "### CRMShootingMethodBVP Test: " << std::endl;
//	std::cout << std::endl << "Free Space Deflection Test: " << std::endl;
//
//	// set up the test problem inputs
//	double InLength = 104.0;
//	double oRlist[NUM_FLEX_SEG] = { 1.5875, 1.5875, 1.5875 };
//	double iRlist[NUM_FLEX_SEG] = { 0.9906, 0.9906, 0.9906 };
//	double YoungModlist[NUM_FLEX_SEG] = { 5.3948, 5.3948, 5.3948 };
//	double ShearModlist[NUM_FLEX_SEG] = { 2.3881, 2.3881, 2.3881 };
//	double CoilAlignmentAngles[2][2] = { {0.0, 0.0} , {0.0, 0.0} };
//	double CoilTurnAreaMat[2][9] = { { 1.44, 0.0, 0.0, 0.0, 1.3851, 0.0, 0.0, 0.0, 1.60 } , { 1.65, 0.0, 0.0, 0.0, 1.326, 0.0, 0.0, 0.0, 1.56 } };
//	double ActuationCurrents[NUM_ACT_SET][3] = { {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1} };
//	double SegmentLengths[NUM_SEGMENTS] = { 10.72, 16.14, 15.395, 15.61, 46.135 };
//	double ActMass[NUM_ACT_SET] = {7.7736e-5, 8.0157e-5};
//	double gravity[3] = { 0.0, 0.0, 9.81 };
//	double rho[NUM_SEGMENTS] = { 2.8814e-7, 2.8814e-7, 2.8814e-7, 2.8814e-7, 2.8814e-7 };
//
//	ContactModeType ContactMode = ContactModeType::FREE_TIP;
//	double TipConstraintPoint[3] = { 0.0, 0.0, 0.0 };
//	double TipForce[3] = { 0.0, 0.0, 0.0 };
//
//	// setup the input parameters
//	CRMCatheterModelParams CathParams;
//	CatheterConfiguration CathConfig;
//	CRMShootingMethodParams<double> BVPParams;
//
//	CRMShootingMethodBVP_Prep(B0,gravity,&(x0[12])/*p0*/,&(x0[3])/*R0*/,
//			SegmentLengths, MarkerLoc, iRlist, oRlist, YoungModlist, ShearModlist, ustarlist,
//			CoilAlignmentAngles, CoilTurnAreaMat, rho, ActMass,
//			CathParams, CathConfig);
//
//	CRMConstructShootingMethodParamSet<double>(CathParams, CathConfig, InLength, ActuationCurrents, ContactMode, TipConstraintPoint, TipForce, IntegrationStepSize, BVPParams);
//
//	double u0_initialguess[3]={0.0, 0.0, 0.0};
//	double ftip_initialguess[3] = { 0.0, 0.0, 0.0 };
//
//	// declare the output variables
//    double u0_calc[3];
//	double ftip_calc[3];
//	int localmin;
//	double residualReal[EQNDIMENSION], residualref[EQNDIMENSION], residualinitial[EQNDIMENSION];
//	double RefMarkerPos[NUM_LOCALIZATION_MARKERS][3], InitialMarkerPos[NUM_LOCALIZATION_MARKERS][3];
//
//	int REPS = 100;
//
//	// Get starting timepoint
//	auto start = high_resolution_clock::now();
//
//    // Cosserat Rod Model - Solve the Boundary Value Problem to calculate the equilibrium configuration of the catheter
//	for (int cnt=0; cnt<REPS; cnt++)
//		CRMShootingMethodBVP (	BVPParams, u0_initialguess, ftip_initialguess, u0_calc, ftip_calc, localmin	);
//
//	// Get ending timepoint
//	auto stop = high_resolution_clock::now();
//
//	// Get duration. Substart timepoints to
//	// get duration. To cast it to proper unit
//	// use duration cast method
//	auto duration = duration_cast<microseconds>(stop - start);
//
//	std::cout << std::endl << "Average time taken by BVP Solution in " << REPS << " repetitions: " << duration.count()/REPS << " microseconds" << std::endl;
//
//    double error3[EQNDIMENSION], expectedu0[EQNDIMENSION]={0.0, 0.0, 0.0};
//	double xfReal[NUM_STATES];
//	double x0Real[NUM_STATES] = { 0.0,0.0,0.0,   1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,   0.0, 0.0, 0.0 };
//
//	for(int i=0; i<3; i++) x0Real[i]=u0_calc[i];
//	CRMSolverIVP(BVPParams, x0Real, ftip_calc, false, xfReal, residualReal, ReportedMarkerPos );
//
//	for(int i=0; i<3; i++) x0Real[i]=u0_initialguess[i];
//	CRMSolverIVP(BVPParams, x0Real, ftip_calc, false, xfReal, residualinitial, InitialMarkerPos );
//
//	for(int i=0; i<3; i++) x0Real[i]=0.0;
//	CRMSolverIVP(BVPParams, x0Real, ftip_calc, false, xfReal, residualref, RefMarkerPos );
//
//
//	std::cout << "localmin  " << localmin << std::endl;
//	mSub_AB<3,1>(u0_calc, expectedu0, error3 );
//
//	printMatrix(expectedu0,1,EQNDIMENSION,"Expected u0");
//	printMatrix(u0_calc,1,EQNDIMENSION,"Calculated u0");
//	printMatrix(error3,1,EQNDIMENSION,"error");
//	printMatrix(residualReal,1,EQNDIMENSION,"residual");
//	printMatrix(residualinitial,1,EQNDIMENSION,"residualinitial");
//	printMatrix(residualref,1,EQNDIMENSION,"residualref");
//	std::cout << "Error Norm:" << sqrt(vNormSq<EQNDIMENSION>(error3)) << std::endl;
//	std::cout << "----" << std::endl;
//	printMatrix(&(ReportedMarkerPos[0][0]), NUM_LOCALIZATION_MARKERS, 3, "ReportedMarkerPos");
//	std::cout << "----" << std::endl;
//
//
//	if (!localmin && (sqrt(vNormSq<EQNDIMENSION>(error3)) < 0.02)) {
//		std::cout << "Test passed for Inserted Length =" << InsertedLength << std::endl;
//		std::cout << "###" << std::endl;
//	}
//	else {
//		std::cout << "TEST FAIL: " << "Results do not match!" << std::endl;
//		std::cout << "###" << std::endl;
//		FailFlag++;
//	}
//
//
//	//
//	//  Test code for CRMShootingMethodBVP under tip constraint
//	//
//
//	std::cout << std::endl << "Constrained Tip Deflection Test: " << std::endl << std::endl;
//
//	//double ustarlistnew[Nc][3] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
//	double ustarlistnew[Nc][3] = { 0.000128, -0.000173, 0.0, 0.00567, -0.00822, 0.0, -0.00250, 0.00188, 0.0 }; //{-0.0025, 0.0019, 0.0, 0.0057, -0.0082, 0.0, 0.00013, -0.00017, 0.0};
//	// Let's apply a tip force and calculate where the catheter tip ends up.
//	double TipForce_2[3] = { 0.001, 0.0, 0.0 };
//	double TipForce_2_Real[3]; for (int i = 0; i < 3; i++) TipForce_2_Real[i] = TipForce_2[i];
//	ContactMode = ContactModeType::FREE_TIP;
//	double u0_initialguess_2[3] = { 0.0, 0.0, 0.0 };//{ 0.001, 0.001, 0.001 };
//	double ftip_initialguess_2[3] = { 0.0, 0.0, 0.0 };//{ 0.009, 0.001, 0.001 };
//	CRMShootingMethodBVP_Prep(B0, gravity, &(x0[12])/*p0*/, &(x0[3])/*R0*/,
//		SegmentLengths, MarkerLoc, iRlist, oRlist, YoungModlist, ShearModlist, ustarlistnew,
//		CoilAlignmentAngles, CoilTurnAreaMat, rho, ActMass,
//		CathParams, CathConfig);
//	CRMConstructShootingMethodParamSet<double>(CathParams, CathConfig, InsertedLength, ActuationCurrents, ContactMode, TipConstraintPoint, TipForce_2, IntegrationStepSize, BVPParams);
//	CRMShootingMethodBVP(BVPParams, u0_initialguess_2, ftip_initialguess_2, u0_calc, ftip_calc, localmin);
//	for (int i = 0; i < 3; i++) expectedu0[i] = x0Real[i] = u0_calc[i];
//	CRMSolverIVP(BVPParams, x0Real, TipForce_2_Real, false, xfReal, residualReal, ReportedMarkerPos);
//	printMatrix(&(xf[12]), 1, 3, "Expected tip location (calculated in free mode)");
//	printMatrix(u0_calc, 1, 3, "u0 - calculated");
//	printMatrix(residualReal, 1, 3, "residual");
//	std::cout << "localmin " << localmin << std::endl;
//	std::cout << std::endl;
//	// Let's use this point as the constraint point and let's calculate what the resulting contact force will be
//	ContactMode = ContactModeType::FIXED_TIP;
//	mCopy_AB<3>(&(xf[12]), TipConstraintPoint);
//	CRMConstructShootingMethodParamSet<double>(CathParams, CathConfig, InsertedLength, ActuationCurrents, ContactMode, TipConstraintPoint, TipForce_2, IntegrationStepSize, BVPParams);
//	CRMShootingMethodBVP(BVPParams, u0_initialguess_2, ftip_initialguess_2, u0_calc, ftip_calc, localmin);
//	//CRMShootingMethodBVP(BVPParams, expectedu0, TipForce_2, u0_calc, ftip_calc, &localmin, &errorflag);
//	for (int i = 0; i < 3; i++) x0Real[i] = u0_calc[i];
//	CRMSolverIVP(BVPParams, x0Real, ftip_calc, false, xfReal, residualReal, ReportedMarkerPos);
//	printMatrix(&(xfReal[12]), 1, 3, "Resulting tip location (calculated in fixed tip mode)");
//	printMatrix(residualReal, 1, 3, "residual");
//	std::cout << "localmin  " << localmin << std::endl;
//	double error4[3], error5[3];
//	mSub_AB<3, 1>(u0_calc, expectedu0, error4);
//	mSub_AB<3, 1>(ftip_calc, TipForce_2, error5);
//	printMatrix(u0_calc, 1, 3, "Calculated u0");
//	printMatrix(expectedu0, 1, 3, "Expected u0");
//	printMatrix(error4, 1, 3, "error");
//	std::cout << "Error Norm:" << sqrt(vNormSq<3>(error4)) << std::endl;
//	printMatrix(ftip_calc, 1, 3, "Calculated Tip Force");
//	printMatrix(TipForce_2, 1, 3, "Expected Tip Force");
//	printMatrix(error5, 1, 3, "error");
//	std::cout << "Error Norm:" << sqrt(vNormSq<3>(error5)) << std::endl;
//	std::cout << "localmin " << localmin << std::endl;
//	std::cout << "----" << std::endl;
//	if (!localmin && (sqrt(vNormSq<EQNDIMENSION>(error5)) < 0.02)) {
//		std::cout << "Test passed for Inserted Length =" << InsertedLength << std::endl;
//		std::cout << "###" << std::endl;
//	}
//	else {
//		std::cout << "TEST FAIL: " << "Results do not match!" << std::endl;
//		std::cout << "###" << std::endl;
//		FailFlag++;
//	}
//
//
//	if (FailFlag==0) {
//		std::cout << std::endl << std::endl << "###" << std::endl;
//		std::cout << "ALL TESTS PASSED!..." << std::endl;
//		std::cout << "###" << std::endl;
//	}
//	else {
//		std::cout << std::endl << std::endl << "###" << std::endl;
//		std::cout << FailFlag << " TEST(S) FAILED!..." << std::endl;
//		std::cout << "###" << std::endl;
//	}
//	return (FailFlag);
//
//
//}
//
//
