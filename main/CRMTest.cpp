#include <iostream>
#include <cmath>
#include <chrono>

#include "CRM.hpp"
using namespace std::chrono;

#ifdef ANALYTICAL_SE3_STEP
constexpr int NUM_INTEGRATION_STATES = 3;   // u[0..2]
#else
constexpr int NUM_INTEGRATION_STATES = 15;   // u[0..2],R[0..9],p[0..2]
#endif



template <typename T>
void printMatrix(T *p, int D1, int D2, const char *text) {
    std::cout << text << " --" << std::endl;
    for (int i=0; i<D1; i++) {
        for (int j=0; j<D2; j++) {
            std::cout << *(p+i*D2+j) << " ";
        }
        std::cout << std::endl;
    }
    //std::cout << "----" << std::endl;
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


int RunExample(void);
int RunFKExample(void);
int RunFKDynamicsExample(void);
int RunTests(void);
int ConstrainedTipTest(void);

int main(int argc, char** argv) {

//	int FailFlag = 0;

    //RunExample();
//	RunFKExample();
//	RunFKExample();

    RunFKDynamicsExample();

//	std::cout << std::endl << std::endl << "### Running Tests ..." << std::endl << std::endl;
//	FailFlag += RunTests();
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

//
// Example showing how to use the CRM_ForwardKinematics functions
//		to calculate the catheter forward kinematics for the constrained tip case
// This is the preferred method
//
int RunFKDynamicsExample(void) {

    std::cout << "### RunFKContactExample() --- CRM Dynamics Examples... " << std::endl;
    std::cout << std::endl << "Constrained Tip Deflection Example: " << std::endl << std::endl;

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
    double ustarlist[NUM_FLEX_SEG][3] = { 0.000148626102272827, 0.000094448815853795, 0, 0.0007399849479553773, -0.0002292201697658481, 0};//0.00567, -0.00822, 0.0};

    double area_ = M_PI *(oRlist[0] * oRlist[0] -  iRlist[0] * iRlist[0]);
    double Area[NUM_SEGMENTS] = { area_, area_, area_}; // Area of the tubing

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
    ContactModeType ContactMode = ContactModeType::FREE_TIP;// ContactModeType::FIXED_TIP;
    // External point force (in spatial coordinates) applied at the tip of the catheter (\lambda = 0)  - unit: ??
    //   (this will be used when ContactMode == ContactModeType::FREE_TIP)
    //   not used in the solution since we are in FIXED_TIP mode
    double TipForce[3] = { 0.0, 0.0, 0.0 };
    // The spatial coordinates of the point where the catheter tip is constrained to be
    //   (this will be used when ContactMode == ContactModeType::FIXED_TIP)
    double TipConstraintPoint[3] = { 6.12, 40.57, 93.96 };


    // *** Control Inputs
    // Inserted Length of the catheter (length of the catheter that is inside the heart chamber) - unit: mm
    double InsertedLength  = 0.0;
    for (int i = 0; i < NUM_SEGMENTS; ++i) {
        InsertedLength += SegmentLengths[i];
    }
    double ActuationCurrents[NUM_ACT_SET][3] = {0.0, 0.0, 0.0};// { {0.100, 0.100, 0.100}};

    // *** Numerical Computation Parameters
    // Stepsize used in numerical integration along the length of the catheter during IVP - unit: mm
    double IntegrationStepSize = 0.2;

    /**
 * These are hard coded, need to revise these later
 */
    // we are adding the tubing mass of the coil section to the total mass of actuator:
    double tubing_mass = rho[0] * SegmentLengths[1];
    double ActInertia[NUM_ACT_SET][9];
    for (int i = 0; i < NUM_ACT_SET; ++i)
    {
        double I_zz = 0.5 * (ActMass[i]) * (oRlist[0] * oRlist[0] + iRlist[0] * iRlist[0]);
        double I_xx = 0.25 * (ActMass[i]) * (oRlist[0] * oRlist[0] + iRlist[0] * iRlist[0]) + 1.0 / 12 * (ActMass[i]) * SegmentLengths[2*i+1] * SegmentLengths[2*i+1];
//        std::cout << "Izz " << I_zz << std::endl;
//        std::cout << "I_xx " << I_xx << std::endl;
//        I_xx += 0.001;
        ActInertia[i][0] = I_xx; ActInertia[i][1] = 0.0; ActInertia[i][2] = 0.0;
        ActInertia[i][3] = 0.0; ActInertia[i][4] = I_xx; ActInertia[i][5] = 0.0;
        ActInertia[i][6] = 0.0; ActInertia[i][7] = 0.0; ActInertia[i][8] = I_zz;
    }
    double v_L_pre[3] = { 0.0, 0.0, 0.0 };
    double w_L_pre[3] = { 0.0, 0.0, 0.0 };

    // *** Storage for storing localization marker positions
    double ReportedMarkerPos[NUM_LOCALIZATION_MARKERS][3];

    // Define variables of convenience that will be use to package catheter physical parameters and spatial configuration parameters
    CRMCatheterModelParams CathParams{};		// Catheter physical parameters will be packaged into this structure
    CatheterConfiguration CathConfig{};		// Catheter spatial configuration parameters will be packages into this structure


    // Package catheter physical parameters and spatial configuration parameters
    //   this step would typically needs to be executed only once
    CRMShootingMethodBVP_Prep(B0, gravity, p0, R0,
                              SegmentLengths, MarkerLoc, iRlist, oRlist, YoungModlist, ShearModlist, ustarlist,
                              CoilAlignmentAngles, CoilTurnAreaMat, rho, ActMass,
                              CathParams, CathConfig);

    // Define initial guesses to be used when solving boundary value problem
    double u0_initialguess[3] = { 0.000709849479553773, -0.000302201697658481, 0};
    // initial guess for the contstraint force at the catheter tip (this will be used when ContactMode == ContactModeType::FIXED_TIP)
    double ftip_initialguess[3] = { 0.0, 0.0, 0.0 };
    // Define initial guesses to be used when solving boundary value problem
    double nL_initialguess[3] = { 0.0, 0.0, 0.0 };

    double mL_initialguess[3] = { 0.0, 0.0, 0.0 };
    // numerical nonlinear equation solver diagnostic outputs
    int localmin;

    double damping[6] = {0, 0, 0, 0, 0, 0};
    double DELTA_T = 0.01;

    /**
     * get the p, R at coil
     */
    // Declare output variables
    // catheter shape state at the entry point of the catheter
    //   states are packed u[0..2], R[0..8], p[0..2](R: 3x3 matrix stored in row major order R11 R12 R13 R21 R22 R23 R31 R32 R33)
    // moment residual at the catheter tip - this should converge to {0,0,0} if the catheter is at its equilibrium configuration
    double residual[NUM_RESIDUAL];
    // spatial coordinates of the localization markers
    double xf[NUM_STATES];
    double u0_calc[3], ftip_calc[3];

    CRMShootingMethod_DYNParams<double> BVPParams{};

//     Populate the variable of convenience for conveniently passing lots of arguments to the BVP Solver
//	   this step would typicall need to be executed every time BVP needs to be calculated as actuation variables would change at every time step
    CRMConstructShootingMethodParamSet<double>(CathParams, CathConfig, InsertedLength, ActuationCurrents, ContactMode,
                                               TipConstraintPoint, TipForce, IntegrationStepSize, ActInertia,
                                               v_L_pre, w_L_pre, p0,  R0, damping, DELTA_T,
                                               BVPParams);

    // Cosserat Rod Model - Solve the Boundary Value Problem to calculate the equilibrium configuration of the catheter
    CRMShootingMethodBVP(BVPParams, u0_initialguess, ftip_initialguess, u0_calc, ftip_calc, localmin);

//    // output for time advance, not used in BVP, just placeholders
    double pL[3], RL[9], TBcoil[3];
//
//    // Cosserat Rod Model - Solve the Initial Value Problem to calculate the shape of the catheter
    double in_mL[3], in_nL[3];
//    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
//        for (int j = 0; j < 3; ++j) {
//            in_u0[i][j] = 0;
//            in_mL[i][j] = 0;
//            in_nL[i][j] = 0;
//        }
//    }
//
//    for (int i = 0; i < 3; ++i) {in_u0[0][i] = u0_calc[i];}
    for (int i = 0; i < 3; ++i) {
        in_mL[i] = 0;
        in_nL[i] = 0;
    }
    CRMSolverIVP(BVPParams, u0_calc, in_mL, in_nL, ftip_calc, true, xf, residual,TBcoil, pL, RL, ReportedMarkerPos);

    //
    // We will use the new Forward Kinematics Function
    //
    // The Forward Kinematics Function uses a different parameter structure than BVP
    CRMDYNFKData<double> FKParams{};
    FKParams.CathConfig = &CathConfig;
    FKParams.CathParams = &CathParams;
    FKParams.ContactMode = ContactMode;
    FKParams.FinalValueOnly = true;  // we want the localization coil locations, too
    FKParams.ReportedMarkerPos = &ReportedMarkerPos;
    for (int i = 0; i < 3; i++) FKParams.TipConstraintPoint[i] = TipConstraintPoint[i];
    for (int i = 0; i < 3; i++) FKParams.TipForce[i] = TipForce[i];
    for (int i = 0; i < 3; i++) FKParams.u0_initialguess[i] = u0_initialguess[i];
    for (int i = 0; i < 3; i++) FKParams.mL_initialguess[i] = mL_initialguess[i];
    for (int i = 0; i < 3; i++) FKParams.nL_initialguess[i] = nL_initialguess[i];
    for (int i = 0; i < 3; i++) FKParams.ftip_initialguess[i] = ftip_initialguess[i];
    FKParams.IntegrationStepSize = IntegrationStepSize;
    for (int i = 0; i < 3; i++) FKParams.v_L_pre[i] = v_L_pre[i];
    for (int i = 0; i < 3; i++) FKParams.w_L_pre[i] = w_L_pre[i];
    for (int i = 0; i < 3; i++) FKParams.pL_pre[i] = pL[i];
    for (int i = 0; i < 9; i++) FKParams.RL_pre[i] = RL[i];
    for (int i = 0; i < 6; i++) FKParams.damping[i] = damping[i];
    FKParams.DELTA_T = DELTA_T;

    for (int i = 0; i < NUM_ACT_SET; ++i) {
        for (int j = 0; j < 9; ++j) FKParams.actInertia[i][j] = ActInertia[i][j];
    }

    // inputs
    ActuationCurrents[0][2] = 0.01;

    VectorXd control_inputs(NUM_ACT_SET * 3 + 1); // actuator currents (distal to proximal) followed by inserted length
    for (int i = 0; i < NUM_ACT_SET; i++)  for (int j = 0; j < 3; j++) control_inputs(i * 3 + j) = ActuationCurrents[i][j];
    control_inputs(NUM_ACT_SET * 3) = InsertedLength;
    // and outputs
    VectorXd FKsolution(out_Dim); // p[0..2],R[0..8],u0[0..2](,ftip[0..2])  R: in row major order

    // Multiple repetitions to more reliably measure time
    int REPS = 100;
    // Get starting timepoint
    auto start = high_resolution_clock::now();

    // Cosserat Rod Model - Solve the Forward Kinematics for the Constrained Tip Case
//	for (int cnt = 0; cnt < REPS; cnt++)
    FKsolution = CRM_DynamicsFK_FreeSpace<double>(control_inputs, FKParams, localmin);

    // Get ending timepoint
    auto stop = high_resolution_clock::now();
    // Get duration. Substart timepoints to
    // get duration. To cast it to proper unit
    // use duration cast method
    auto duration = duration_cast<microseconds>(stop - start);
//	std::cout << std::endl << "Average time taken by Constrained Tip FK Solution in " << REPS << " repetitions: " << duration.count() / REPS << " microseconds" << std::endl << std::endl;

    //display the results,
    std::cout << "FK Output -- v,w,p,R: \n" << FKsolution.transpose() << "\n";
    std::cout << "localmin:" << localmin << std::endl;
    std::cout << "----" << std::endl;

    return (localmin);
}

