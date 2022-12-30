#include <iostream>
#include <cmath>
#include <chrono>
#include <cstring>

#include<fstream>

#include "CRMTest.h"
#include "CRMDYN.hpp"
#include <vector>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

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



void read_input_file(const char *path_to_input, double input_currents[][3], double time_step_[]){
    FILE *file = fopen(path_to_input, "r");
    if (file == NULL) {
        printf("Impossible to load input currents! Are you in the right path ?\n");
    }
    std::vector< float > input_vec, time_vec;
    float x;
    std::ifstream File;
    File.open(path_to_input);
    while (File >> x) {
        input_vec.push_back(x);
        if (File.peek() == ',')
            File.ignore();
    }

    int SIZE_CUR_SET = input_vec.size() / 9;
    for (int i = 0; i < SIZE_CUR_SET; ++i) {
        for (int j = 0; j < 3; ++j) {
            input_currents[i][j] = 0.001 * input_vec[3 * (3 * i + j) + 1];
        }
        time_step_[i] = 0.001 * input_vec[9*(i+1)-1];
    }

//    for (int i = 0; i < size_input; ++i) {
//        std::cout << "input_currents: " << input_currents[i][0] << " " << input_currents[i][1]  << " " << input_currents[i][2] << std::endl;
//        std::cout << time_step_[i] << std::endl;
//    }
//    std::cout << "size_input: " << size_input << std::endl;
}

void display_points(std::vector<std::vector<double> > &input_pts){

    int size_ = input_pts.size();
    std::vector<double> v1, v2, v3;
    for (int id_ = 0; id_ < size_; ++id_) {
        v1.push_back(input_pts[id_][0]);
        v2.push_back(input_pts[id_][1]);
        v3.push_back(input_pts[id_][2]);
    }

    plt::scatter(v1, v2, v3,3);
    plt::xlabel("x");
    plt::ylabel("y");
    plt::set_zlabel("z"); // set_zlabel rather than just zlabel, in accordance with the Axes3D method
    plt::show();
    plt::pause(0);

}


int RunTests(void);
int RunExample(void);
//int RunDynamicsExample(int SIZE_CUR_SET, double currents_set[][3], double time_step[], std::vector<std::vector<double> > &output_tip_pts) ;

int main(int argc, char** argv) {

//    const char *path_to_input = "/home/ranhao/Documents/CRM_Dynamics/input_file/circle_01_new.play";
//
//    FILE *file = fopen(path_to_input, "r");
//    if (file == NULL) {
//        printf("Impossible to load input currents! Are you in the right path ?\n");
//    }
//    std::vector< float > input_vec;
//    float x;
//    std::ifstream File;
//    File.open(path_to_input);
//    while (File >> x) {
//        input_vec.push_back(x);
//        if (File.peek() == ',')
//            File.ignore();
//    }
//
//    int SIZE_CUR_SET = input_vec.size() / 9;
//    double input_currents[SIZE_CUR_SET][3], time_step_[SIZE_CUR_SET];
//
//    for (int i = 0; i < SIZE_CUR_SET; ++i) {
//        for (int j = 0; j < 3; ++j) {
//            input_currents[i][j] = 0.001 * input_vec[3 * (3 * i + j) + 1];
//        }
//        time_step_[i] = 0.001 * input_vec[9*(i+1)-1];
//    }

//    double sum = 0;
//    for (int i = 2; i < SIZE_CUR_SET; ++i) {
//        sum += time_step_[i];
//    }

//    for (int i = 0; i < SIZE_CUR_SET; ++i) {
//        std::cout << "input_currents: " << input_currents[i][0] << " " << input_currents[i][1]  << " " << input_currents[i][2] << std::endl;
//        std::cout << time_step_[i] << std::endl;
//    }
//    std::cout << "size_input: " << SIZE_CUR_SET << std::endl;

//    std::vector< std::vector< double > > tip_pts;
//    tip_pts.resize(SIZE_CUR_SET);
//    RunDynamicsExample(SIZE_CUR_SET, input_currents, time_step_, tip_pts);
//
//    for (int i = 0; i < SIZE_CUR_SET; ++i) {
//        std::cout << "input_currents: " << input_currents[i][0] << " " << input_currents[i][1]  << " " << input_currents[i][2] << std::endl;
//    }
//
//    int output_length = tip_pts.size();
//    for (int i = 0; i < output_length; ++i) {
//        std::cout << "output: " << tip_pts[i][0] << " " << tip_pts[i][1]  << " " << tip_pts[i][2] << std::endl;
//
//    }
//    display_points(tip_pts);

	RunExample();



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
    double YoungModlist[NUM_FLEX_SEG] =  { 5.3948, 5.3948 };// { 4.9105, 4.9105 };
    // Shear moduli of each of the flexbile segments - unit: ??
    double ShearModlist[NUM_FLEX_SEG] = { 2.3881, 2.3881 };//{ 1.6545, 1.6545 };
    // Alignment angles for the side coils, for each of the actuator sets - unit: radians
    double CoilAlignmentAngles[NUM_ACT_SET][2] = { {0.0, 0.0}  };// { {-0.0871, -0.3934}  };
    // Coil turn area matrices for each of the actuator sets  - unit: mm2
    double CoilTurnAreaMat[NUM_ACT_SET][9] = { { 1.3851,  0.0, 0.0, 0.0, 1.44, 0.0, 0.0, 0.0, 1.60 } };
    // Lengths of each of the catheter segments - unit: mm
    double SegmentLengths[NUM_SEGMENTS] = { 19.85,  18.3, 59.40 };

    // Mass of each of the actuator sets - unit: ??
    double ActMass[NUM_ACT_SET] = { 5.7736e-5};
    // Array of lambda values for marker locations  (distance from the tip to each of the markers) - unit: mm
    double MarkerLoc[NUM_LOCALIZATION_MARKERS] = { 2.18, 18.79, 50.06, 101.52, 104.02 };
    // Length density for each of the catheter segments
    double rho[NUM_SEGMENTS] = { 2.8814e-7, 2.8814e-7, 2.8814e-7 };
    // Local curvature in unloaded configuration for each of the flexible segments; NUM_FLEX_SEG*3 long array corresponding to NUM_FLEX_SEG many 3x1 vectors
    double ustarlist[NUM_FLEX_SEG][3] = { 0.000148626102272827, 0.00094448815853795, 0, 0.000799849479553773, -0.0001892201697658481, 0};
//    double ustarlist[NUM_FLEX_SEG][3] = { 0.000148626102272827, 0.00094448815853795, 0,  0.000128, -0.000173, 0.0};//0.00567, -0.00822, 0.0};

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
    // Actuation currents for each of the coils for each of the coil sets - unit: A
    double ActuationCurrents[NUM_ACT_SET][3] = {0.0, 0.0, 0.0};// { {0.100, 0.100, 0.100}};


	// *** Numerical Computation Parameters
	// Stepsize used in numerical integration along the length of the catehter during IVP - unit: mm
	double IntegrationStepSize = 0.2;

//    double tubing_mass = rho[0] * SegmentLengths[1];
//    ActMass[0] += tubing_mass;
    double ActInertia[NUM_ACT_SET][9];
    for (int i = 0; i < NUM_ACT_SET; ++i)
    {
        double I_zz = 0.5 * ( ActMass[i]) * (oRlist[0] * oRlist[0] + iRlist[0] * iRlist[0]);
        double I_xx = 0.25 * ( ActMass[i]) * (oRlist[0] * oRlist[0] + iRlist[0] * iRlist[0]) + 1.0 / 12 * ActMass[i] * SegmentLengths[2*i+1] * SegmentLengths[2*i+1] ;
        ActInertia[i][0] = I_xx; ActInertia[i][1] = 0.0; ActInertia[i][2] = 0.0;
        ActInertia[i][3] = 0.0; ActInertia[i][4] = I_xx; ActInertia[i][5] = 0.0;
        ActInertia[i][6] = 0.0; ActInertia[i][7] = 0.0; ActInertia[i][8] = I_zz;
    }
    double v_L_pre[3] = { 0.0, 0.0, 0.0 };
    double w_L_pre[3] = { 0.0, 0.0, 0.0 };

    double v0[3] = {0.0,0.0,0.0};              //initial linear velocity in local frame of entry point
    double w0[3] = {0.0,0.0,0.0};

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

    paramHistory u_history[NUM_FLEX_SEG], v_history[NUM_FLEX_SEG], w_history[NUM_FLEX_SEG];
    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
        int length = ( SegSteps[i] + 1 ) * 3 ;
        if (length < NUM_HISTORY_LENGTH){
            u_history[i].length = length;
            for (int j = 0; j < ( SegSteps[i] + 1 ) ; ++j) {
                u_history[i].data[j*3] = ustarlist[NUM_FLEX_SEG -1 - i][0];
                u_history[i].data[j*3+1] = ustarlist[NUM_FLEX_SEG -1 -i][1];
                u_history[i].data[j*3+2] = ustarlist[NUM_FLEX_SEG -1 -i][2];
            }

            v_history[i].length = length;
            v_history[i].data[0] = v0[0];
            v_history[i].data[1] = v0[1];
            v_history[i].data[2] = v0[2];
            for (int j = 1; j < ( SegSteps[i] + 1 ) ; ++j) {
                v_history[i].data[j*3] = 0.0;
                v_history[i].data[j*3+1] = 0.0;
                v_history[i].data[j*3+2] = 0.0;
            }

            w_history[i].length = length;
            w_history[i].data[0] = w0[0];
            w_history[i].data[1] = w0[1];
            w_history[i].data[2] = w0[2];
            for (int j = 1; j < ( SegSteps[i] + 1 ) ; ++j) {
                w_history[i].data[j*3] = 0.0;
                w_history[i].data[j*3+1] = 0.0;
                w_history[i].data[j*3+2] = 0.0;
            }
        }else{std::cout <<" Maximum history length exceeded!";
            exit(3);
        }


    }


//
//    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
//        for (int j = 0; j < ( SegSteps[i] + 1 ) ; ++j) {
//            std::cout << " u_history " << u_history[i].data[j*3] << " " << u_history[i].data[j*3+1] << " " << u_history[i].data[j*3+2] << std::endl;
//            std::cout <<"end of first round: " << j*3+2 << std::endl;
//        }
//        std::cout << "length_of u : " <<  u_history[i].length << std::endl;
//
//    }

//    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
//        for (int j = 0; j < ( SegSteps[i] + 1 ) ; ++j) {
//            std::cout << " w_history " << w_history[i].data[j*3] << " " << w_history[i].data[j*3+1] << " " << w_history[i].data[j*3+2] << std::endl;
//            std::cout <<"end of first round: " << j*3+2 << std::endl;
//        }
//        std::cout << "length_of w: " <<  w_history[i].length << std::endl;
//
//    }



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

    // Define variables of convenience that will be use to package catheter physical parameters and spatial configuration parameters
    CRMCatheterModelParams CathParams{};		// Catheter physical parameters will be packaged into this structure
    CatheterConfiguration CathConfig{};		// Catheter spatial configuration parameters will be packages into this structure


    // Package catheter physical parameters and spatial configuration parameters
	//   this step would typically needs to be executed only once
    CRMShootingMethodBVP_Prep(B0,gravity, p0, R0, v0, w0,
                              SegmentLengths, MarkerLoc, iRlist, oRlist, YoungModlist, ShearModlist, ustarlist,
                              CoilAlignmentAngles, CoilTurnAreaMat, rho, ActMass, ActInertia, tubingInertia,
                              CathParams, CathConfig);

	// Define initial guesses to be used when solving boundary value problem
	// initial guess for the curvature at the catheter base
	double u0_initialguess[3] = { 0.0007399849479553773, -0.00001892201697658481, 0};// NEED TO START IN A CLOSE NEIGHBORHOOD
	// initial guess for the contstraint force at the catheter tip (this will be used when ContactMode == ContactModeType::FIXED_TIP)
	double ftip_initialguess[3] = { 0.0, 0.0, 0.0 };

    double n0_initialguess[3] = { 0.0, 0.0, 0.0 };

    double nL_initialguess[3] = { 0.0, 0.0, 0.0 };

    double uL_initialguess[3] = { 0.0000148626102272827, 0.00094448815853795, 0 };

    double initial_guesses[NUM_RESIDUAL];
    for (int i = 0; i < 3; ++i) {
        initial_guesses[i] = u0_initialguess[i];
        initial_guesses[i+3] = n0_initialguess[i];
    }
    if (NUM_RESIDUAL > 6){
        for (int i = 0; i < 3; ++i) {
            initial_guesses[i+6] = uL_initialguess[i];
            initial_guesses[i+9] = nL_initialguess[i];
        }
    }

    // Declare the output variables for BVP
    // calculated curvature at the catheter base
    double out_calc[NUM_RESIDUAL];
    // calculated contstraint force at the catheter tip (this will be used when ContactMode == ContactModeType::FIXED_TIP)
    double ftip_calc[3];
    // numerical nonlinear equation solver diagnostic outputs
    int localmin;

    double DELTA_T = 0.05;
    double damping_tubing[3] = {0,0,0};
    double damping_coil[6] = {0,0,0,0,0,0};
	// Define the variable convenience used to package the arguments passed to BVP Solver
	CRMShootingMethodParams<double> BVPParams{};
	// Populate the variable of convenience for conveniently passing lots of arguments to the BVP Solver
	//   this step would typicall need to be executed every time BVP needs to be calculated as actuation variables would change at every time step
	CRMConstructShootingMethodParamSet<double>(CathParams, CathConfig, InsertedLength, ActuationCurrents, ContactMode, TipConstraintPoint, TipForce,
                                               IntegrationStepSize, v_L_pre, w_L_pre, u_history, v_history, w_history, DELTA_T, damping_tubing, damping_coil, BVPParams);

    // Cosserat Rod Model - Solve the Boundary Value Problem to calculate the equilibrium configuration of the catheter
	CRMShootingMethodBVP(BVPParams, initial_guesses, ftip_initialguess,
                         out_calc, ftip_calc, localmin);


    paramHistory u_history_update[NUM_FLEX_SEG], v_history_update[NUM_FLEX_SEG], w_history_update[NUM_FLEX_SEG];
    double v_L_update[3], w_L_update[3], pL_update[3], RL_update[9], h0_new[NUM_FLEX_SEG];
	// Declare output variables
	// catheter shape state at the entry point of the catheter
	//   states are packed u[0..2], R[0..8], p[0..2](R: 3x3 matrix stored in row major order R11 R12 R13 R21 R22 R23 R31 R32 R33)
	double xf[NUM_STATES];
	// moment residual at the catheter tip - this should converge to {0,0,0} if the catheter is at its equilibrium configuration
	double residual[NUM_RESIDUAL];
	// spatial coordinates of the localization markers
	double ReportedMarkerPos[NUM_LOCALIZATION_MARKERS][3];

	// Cosserat Rod Model - Solve the Initial Value Problem to calculate the shape of the catheter
    CRMSolverIVP(BVPParams, out_calc, ftip_calc, true, xf, residual, ReportedMarkerPos, u_history_update, v_history_update, w_history_update,
                 v_L_update, w_L_update, pL_update, RL_update, h0_new);

	// Print outputs
	std::cout << "localmin  " << localmin << std::endl;
	printMatrix(&(out_calc[0]), 1, 3, "Calculated curvature at base");
    printMatrix(&(out_calc[3]), 1, 3, "Calculated internal force at base");
    if (NUM_RESIDUAL > 6){
        printMatrix(&(out_calc[6]), 1, 3, "Calculated curvature at coil top");
        printMatrix(&(out_calc[9]), 1, 3, "Calculated force at coil top");
    }


    printMatrix(ftip_calc, 1, 3, "Calculated Tip Force");
	printMatrix(&(xf[12]), 1, 3, "Catheter Tip Position");
	printMatrix(&(xf[3]), 1, 9, "Catheter Tip Orientation");

    printMatrix(&(xf[15]), 1, 3, "v at tip");
    printMatrix(&(xf[18]), 1, 3, "w at tip");

    printMatrix(v_L_update, 1, 3, "v at coil");
    printMatrix(w_L_update, 1, 3, "w at coil");

    printMatrix(pL_update, 1, 3, "p at coil");
    printMatrix(RL_update, 1, 9, "R at coil");

    printMatrix(residual, 1, 3, "Force residual at coil");
    printMatrix(&(residual[3]), 1, 3, "Moment residual at coil");



//    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
//        for (int j = 0; j < ( SegSteps[i] + 1 ) ; ++j) {
//            std::cout << " v_history_update " << v_history_update[i].data[j*3] << " " << v_history_update[i].data[j*3+1] << " " << v_history_update[i].data[j*3+2] << std::endl;
//            std::cout <<"end of first round: " << j*3+2 << std::endl;
//        }
//        std::cout << "length_of v_history_update : " <<  v_history_update[i].length << std::endl;
//
//    }


    std::cout << "----" << std::endl;


	return (localmin);


}
//
//int RunDynamicsExample(int SIZE_CUR_SET, double currents_set[][3], double time_step[], std::vector<std::vector<double> > &output_tip_pts) {
//
//    std::cout << "############################ " << std::endl;
//    std::cout << "### CRM Dynamics Examples... " << std::endl;
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
//    // Young's moduli of each of the flexible segments - unit: ??
//    double YoungModlist[NUM_FLEX_SEG] = { 5.3948, 5.3948 };
//    // Shear moduli of each of the flexbile segments - unit: ??
//    double ShearModlist[NUM_FLEX_SEG] = { 2.3881, 2.3881 };
//    // Alignment angles for the side coils, for each of the actuator sets - unit: radians
//    double CoilAlignmentAngles[NUM_ACT_SET][2] = { {0.0, 0.0}  };
//    // Coil turn area matrices for each of the actuator sets  - unit: mm2
//    double CoilTurnAreaMat[NUM_ACT_SET][9] = { { 1.44, 0.0, 0.0, 0.0, 1.3851, 0.0, 0.0, 0.0, 1.60 } };
//    // Lengths of each of the catheter segments : Tip to Base Ordering  - unit: mm
//    double SegmentLengths[NUM_SEGMENTS] = { 19.85,  18.3, 60.35 };
////    double SegmentLengths[NUM_SEGMENTS] = { 10.72, 16.14, 15.395 };
//
//    // Mass of each of the actuator sets - unit: ??  does this include the tubing inside??????
//    double ActMass[NUM_ACT_SET] = { 7.7736e-5};
//    // Array of lambda values for marker locations  (distance from the tip to each of the markers) - unit: mm
//    double MarkerLoc[NUM_LOCALIZATION_MARKERS] = { 2.18, 18.79, 50.06, 60.52, 80.02 };
//    // Length density for each of the catheter segments
//    double rho[NUM_SEGMENTS] = { 2.8814e-7, 2.8814e-7, 2.8814e-7 };
//    // Local curvature in unloaded configuration for each of the flexible segments; NUM_FLEX_SEG*3 long array corresponding to NUM_FLEX_SEG many 3x1 vectors
//    double ustarlist[NUM_FLEX_SEG][3] = { 0.000128, -0.000173, 0.0, 0.000128, -0.000173, 0.0}; //0.00567, -0.00822, 0.0};
//
//    double area_ = M_PI *(oRlist[0] * oRlist[0] -  iRlist[0] * iRlist[0]);
//    double Area[NUM_SEGMENTS] = { area_, area_, area_}; // Area of the tubing
//
//
////    // *** Physical Description of the Catheter
////    // IMPORTANT NOTE: For now, most proximal segment is assumed to be always flexible
////    //    and the flexible and rigid segments are assumed to be alternating
////    //    most distal segment can be flexible or rigid
////    // Items are listed in distal-to-proximal order (starting from the tip of the catheter towards the base)
////    // Outer radii of each of the flexible segments - unit: mm
////    double oRlist[NUM_FLEX_SEG] = { 1.5875  };
////    // Inner radii of each of the flexible segments - unit: mm
////    double iRlist[NUM_FLEX_SEG] = { 0.9906 };
////    // Young's moduli of each of the flexible segments - unit: ??
////    double YoungModlist[NUM_FLEX_SEG] = { 5.3948 };
////    // Shear moduli of each of the flexbile segments - unit: ??
////    double ShearModlist[NUM_FLEX_SEG] = { 2.3881 };
////    // Alignment angles for the side coils, for each of the actuator sets - unit: radians
////    double CoilAlignmentAngles[NUM_ACT_SET][2] = { {0.0, 0.0}  };
////    // Coil turn area matrices for each of the actuator sets  - unit: mm2
////    double CoilTurnAreaMat[NUM_ACT_SET][9] = { { 1.44, 0.0, 0.0, 0.0, 1.3851, 0.0, 0.0, 0.0, 1.60 } };
////    // Lengths of each of the catheter segments : Tip to Base Ordering  - unit: mm
////    double SegmentLengths[NUM_SEGMENTS] = {  18.3, 60.35 };
//////    double SegmentLengths[NUM_SEGMENTS] = { 10.72, 16.14, 15.395 };
////
////    // Mass of each of the actuator sets - unit: ??  does this include the tubing inside??????
////    double ActMass[NUM_ACT_SET] = { 7.7736e-5};
////    // Array of lambda values for marker locations  (distance from the tip to each of the markers) - unit: mm
////    double MarkerLoc[NUM_LOCALIZATION_MARKERS] = { 2.18, 18.79, 50.06, 60.52, 80.02 };
////    // Length density for each of the catheter segments
////    double rho[NUM_SEGMENTS] = { 2.8814e-7, 2.8814e-7};
////    // Local curvature in unloaded configuration for each of the flexible segments; NUM_FLEX_SEG*3 long array corresponding to NUM_FLEX_SEG many 3x1 vectors
////    double ustarlist[NUM_FLEX_SEG][3] = { 0.000128, -0.000173, 0.0}; //0.00567, -0.00822, 0.0};
////
////    double area_ = M_PI *(oRlist[0] * oRlist[0] -  iRlist[0] * iRlist[0]);
////    double Area[NUM_SEGMENTS] = { area_, area_}; // Area of the tubing
//
//
//
//    // *** Catheter Configuration in spatial coordinates
//    // B0 field vector of the MRI scanner (in spatial coordinates) - unit: Tesla
//    double B0[3] = { 0.0, 3.0, 0.0 };
//    // Gravity vector - unit: ??
//    double gravity[3] = { 0.0, 0.0, 9.81 };
//    // Catheter entry point coordinates (in spatial coordinates) - unit: mm
//    double p0[3] = { 0.0, 0.0, 0.0 };
//    // Catheter entry point orientation (3x3 rotation matrix describing catheter entry point frame orientation relative to the spatial frame stored in row major order) - unit: unitless
//    double R0[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
//
//    // *** Other External variables
//    // specify if catheter is in free space or if the catheter tip is constrained to a contact point
//    ContactModeType ContactMode = ContactModeType::FREE_TIP;
//    // External point force (in spatial coordinates) applied at the tip of the catheter (\lambda = 0)  - unit: ??
//    //   (this will be used when ContactMode == ContactModeType::FREE_TIP)
//    double TipForce[3] = { 0.0, 0.0, 0.0 };
//    // The spatial coordinates of the point where the catheter tip is constrained to be
//    //   (this will be used when ContactMode == ContactModeType::FIXED_TIP)
//    double TipConstraintPoint[3] = { 0.0, 0.0, 0.0 };
//
//    // *** Numerical Computation Parameters
//    // Stepsize used in numerical integration along the length of the catehter during IVP - unit: mm
//    double IntegrationStepSize = 0.2;
//
////	double tubing_mass = M_PI * ( oRlist[0] * oRlist[0] -  iRlist[0] * iRlist[0]) * rho[0] * SegmentLengths[1];
//
//    double tubing_mass = rho[0] * SegmentLengths[1];
//    ActMass[0] += tubing_mass;
//    double ActInertia[NUM_ACT_SET][9];
//    for (int i = 0; i < NUM_ACT_SET; ++i)
//    {
//        double I_zz = 0.5 * ( ActMass[i]) * (oRlist[0] * oRlist[0] + iRlist[0] * iRlist[0]);
//        double I_xx = 0.25 * ( ActMass[i]) * (oRlist[0] * oRlist[0] + iRlist[0] * iRlist[0]) + 1.0 / 12 * ActMass[i] * SegmentLengths[2*i+1] * SegmentLengths[2*i+1] ;
//        ActInertia[i][0] = I_xx; ActInertia[i][1] = 0.0; ActInertia[i][2] = 0.0;
//        ActInertia[i][3] = 0.0; ActInertia[i][4] = I_xx; ActInertia[i][5] = 0.0;
//        ActInertia[i][6] = 0.0; ActInertia[i][7] = 0.0; ActInertia[i][8] = I_zz;
//    }
//    double v_L_pre[3] = { 0.0, 0.0, 0.0 };
//    double w_L_pre[3] = { 0.0, 0.0, 0.0 };
//
//    double v0[3] = {0.0,0.0,0.0};              //initial linear velocity in local frame of entry point
//    double w0[3] = {0.0,0.0,0.0};
//
//
//    /**
//     * initial curvature is assumed constant curvature, we are initializing the full length u_history,
//     * Prepare for interpolation of dynamic insertion length
//     */
//    double SegEnds[NUM_SEGMENTS+1];
//    SegEnds[0] = 0.0;
//    for (int i = 1; i < NUM_SEGMENTS+1; ++i) {
//        SegEnds[i] = SegEnds[i-1] + SegmentLengths[NUM_SEGMENTS - i]; // SegmentLengths is ordered from the distal
//    }
//
//    int SegSteps[NUM_FLEX_SEG];
//    double h0[NUM_FLEX_SEG];
//    for (int i=0; i<NUM_FLEX_SEG; i++) { 		// flexible catheter segment
//        // Calculate the number of integration steps based on the given IntegrationStepSize
//        SegSteps[i]=int(ceil( (SegEnds[2*i+1]-SegEnds[2*i]) / IntegrationStepSize ) );
//        h0[i] = dVal((SegEnds[2*i+1]-SegEnds[2*i]))/(SegSteps[i]*1.0);
//    }
//
//    paramHistory u_history[NUM_FLEX_SEG], v_history[NUM_FLEX_SEG], w_history[NUM_FLEX_SEG];
//    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
//        int length = ( SegSteps[i] + 1 ) * 3 ;
//        u_history[i].length = length;
//        u_history[i].data = new double[length];
//        for (int j = 0; j < ( SegSteps[i] + 1 ) ; ++j) {
//            u_history[i].data[j*3] = ustarlist[NUM_FLEX_SEG -1 - i][0];
//            u_history[i].data[j*3+1] = ustarlist[NUM_FLEX_SEG -1 -i][1];
//            u_history[i].data[j*3+2] = ustarlist[NUM_FLEX_SEG -1 -i][2];
//        }
//
//        v_history[i].length = length;
//        v_history[i].data = new double[length];
//        v_history[i].data[0] = v0[0];
//        v_history[i].data[1] = v0[1];
//        v_history[i].data[2] = v0[2];
//        for (int j = 1; j < ( SegSteps[i] + 1 ) ; ++j) {
//            v_history[i].data[j*3] = 0.0;
//            v_history[i].data[j*3+1] = 0.0;
//            v_history[i].data[j*3+2] = 0.0;
//        }
//
//        w_history[i].length = length;
//        w_history[i].data = new double[length];
//        w_history[i].data[0] = w0[0];
//        w_history[i].data[1] = w0[1];
//        w_history[i].data[2] = w0[2];
//        for (int j = 1; j < ( SegSteps[i] + 1 ) ; ++j) {
//            w_history[i].data[j*3] = 0.0;
//            w_history[i].data[j*3+1] = 0.0;
//            w_history[i].data[j*3+2] = 0.0;
//        }
//    }
//
//    // Area Moment of Inertia for a hollow cylindrical section
//    // https://en.wikipedia.org/wiki/Second_moment_of_area
//    // https://www.engineeringtoolbox.com/area-moment-inertia-d_1328.html
//    double tubingInertia[NUM_FLEX_SEG][9];
//
//    for (int i = 0; i < NUM_FLEX_SEG; ++i)
//    {
//        double I_xx =  M_PI * (oRlist[i] * oRlist[i] * oRlist[i] * oRlist[i] - iRlist[i] * iRlist[i] * iRlist[i] * iRlist[i] )  * 0.25;
//        double I_zz = 2 * I_xx;
//        tubingInertia[i][0] = I_xx; tubingInertia[i][1] = 0.0; tubingInertia[i][2] = 0.0;
//        tubingInertia[i][3] = 0.0; tubingInertia[i][4] = I_xx; tubingInertia[i][5] = 0.0;
//        tubingInertia[i][6] = 0.0; tubingInertia[i][7] = 0.0; tubingInertia[i][8] = I_zz;
//    }
//
//
//
//
//    // Define variables of convenience that will be use to package catheter physical parameters and spatial configuration parameters
//    CRMCatheterModelParams CathParams{};		// Catheter physical parameters will be packaged into this structure
//    CatheterConfiguration CathConfig{};		// Catheter spatial configuration parameters will be packages into this structure
//    CRMShootingMethodParams<double> BVPParams{};
//    CRMDynamicsParams<double> DynamicsParams{};
//
//
//    CRMShootingMethodBVP_Prep(B0,gravity, p0, R0, v0, w0,
//                              SegmentLengths, MarkerLoc, iRlist, oRlist, YoungModlist, ShearModlist, ustarlist,
//                              CoilAlignmentAngles, CoilTurnAreaMat, rho, ActMass, ActInertia, tubingInertia,
//                              CathParams, CathConfig);
//
//    // Define initial guesses to be used when solving boundary value problem
//    // initial guess for the curvature at the catheter base
//    double u0_initialguess[3] = { 0.000128, -0.000173, 0.0} ; //{ 0.000128, -0.000173, 0.0 }
//    // initial guess for the internal force at the catheter base
//    double n0_initialguess[3] = { 0.0, 0.0, 0.0 };
//    // initial guess for the contstraint force at the catheter tip (this will be used when ContactMode == ContactModeType::FIXED_TIP)
//    double ftip_initialguess[3] = { 0.0, 0.0, 0.0 };
//
////    // initial guess for the internal force applied to the next flexible segment at upper side of the coil (L)
////    double nL_initialguess[3] = { 0.0, 0.0, 0.0 };
////    // initial guess for the internal moment applied to the next flexible segment at upper side of the coil (L)
////    double mL_initialguess[3] = { 0.0, 0.0, 0.0 };
//
//    double ReportedMarkerPos[NUM_LOCALIZATION_MARKERS][3];
//
//    DynamicsParams.CathParams = &CathParams;
//    DynamicsParams.ContactMode = ContactMode;
//    DynamicsParams.FinalValueOnly = false;  // we want the localization coil locations, too
//    DynamicsParams.ReportedMarkerPos = &ReportedMarkerPos; //report final shape of the catheter
//
//    for (int i = 0; i < 3; i++) DynamicsParams.TipConstraintPoint[i] = TipConstraintPoint[i];
//    for (int i = 0; i < 3; i++) DynamicsParams.TipForce[i] = TipForce[i];
//    for (int i = 0; i < 3; ++i) {
//        DynamicsParams.initial_guess[i] = u0_initialguess[i];
//        DynamicsParams.initial_guess[i+3] = n0_initialguess[i];
//    }
//    if (NUM_RESIDUAL>6){
//        DynamicsParams.initial_guess[6] = 0.000128; DynamicsParams.initial_guess[7] = -0.000173;DynamicsParams.initial_guess[8] =  0.0; //
//        DynamicsParams.initial_guess[9] =  0.0; DynamicsParams.initial_guess[10] = 0.0;DynamicsParams.initial_guess[11] =  0.0;
//    }
//    for (int i = 0; i < 3; i++) DynamicsParams.ftip_initialguess[i] = ftip_initialguess[i];
//    DynamicsParams.IntegrationStepSize = IntegrationStepSize;
//
//    //*** Control Inputs
//    // Inserted Length of the catheter (length of the catheter that is inside the heart chamber) - unit: mm
//    double InsertedLength  = 0.0;
//    for (int i = 0; i < NUM_SEGMENTS; ++i) {
//        InsertedLength += SegmentLengths[i];
//    }
//
//    double x_t[NUM_STATES], control_t[NUM_CONTROL];
//    for (int i = 0; i < NUM_STATES; i++) {
//        if (i < 3) {
//            x_t[i] = u0_initialguess[i];
//        }
//        else if (i < 12) {
//            x_t[i] = R0[i - 3];
//        }
//        else if (i < 15) {
//            x_t[i] = p0[i - 12];
//        }
//        else if (i < 18) {
//            x_t[i] = v0[i - 15];
//        }
//        else if (i < 21) {
//            x_t[i] = w0[i - 18];
//        }
//        else if (i < 24) {
//            x_t[i] = n0_initialguess[i-21];
//        }
//        else {x_t[i] = 0.0;} //m0 will be calculated inside with curvature
//    }
//
//
//    // Declare the output variables for BVP
//    double out_calc[NUM_RESIDUAL];
//    // calculated contstraint force at the catheter tip (this will be used when ContactMode == ContactModeType::FIXED_TIP)
//    double ftip_calc[3];
//    // numerical nonlinear equation solver diagnostic outputs
//    int localmin;
//    double residual[NUM_RESIDUAL];
//
//    double x_t_tip[NUM_STATES], v_L_new[3], w_L_new[3], h0_pre_new[NUM_FLEX_SEG];
//    paramHistory u_history_new[NUM_FLEX_SEG], v_history_new[NUM_FLEX_SEG], w_history_new[NUM_FLEX_SEG];
//
//    control_t[NUM_CONTROL-1] = InsertedLength;
//
//
//    //initial inputs
//    double ActuationCurrents[NUM_ACT_SET][3] = {0.0, 0.0, 0.0};// { {0.100, 0.100, 0.100}};
//
//    CRMConstructShootingMethodParamSet<double>(CathParams, CathConfig, InsertedLength, ActuationCurrents, ContactMode, TipConstraintPoint, TipForce,
//                                               IntegrationStepSize, v_L_pre, w_L_pre, u_history, v_history, w_history, h0,  BVPParams);
//
//
//    for (int i = 0; i < SIZE_CUR_SET; ++i) {
//        for (int j = 0; j < 3; ++j) {
//            control_t[j] = currents_set[i][j];
//        }
//        control_t[0] = control_t[1] = 0.0;
//        control_t[2] = 0.0;
//        int duration = (int) (time_step[i] / DELTA_T) ;
//
//        for (int t = 0; t < duration; ++t) {
//            //        std::cout << "v_L_pre: " << v_L_pre[0] << " " << v_L_pre[1] << " " << v_L_pre[2] << std::endl;
//            //        std::cout << "w_L_pre: " << w_L_pre[0] << " " << w_L_pre[1] << " " << w_L_pre[2] << std::endl;
//            //        std::cout << "h0: " << h0[0]  << std::endl;
//            std::cout << "control_t : " << control_t[0] << " " << control_t[1] << " " << control_t[2] << std::endl;
//
//            CRM_Dynamics(x_t, control_t, BVPParams, u_history, v_history, w_history, v_L_pre, w_L_pre, h0,  DynamicsParams,
//                         x_t_tip, out_calc, ReportedMarkerPos, u_history_new, v_history_new, w_history_new, v_L_new, w_L_new, h0_pre_new, residual);
//
//            printMatrix(&(x_t_tip[12]), 1, 3, "Catheter Tip Position");
//            printMatrix(&(x_t_tip[15]), 1, 3, "v at tip");
//            printMatrix(&(x_t_tip[18]), 1, 3, "w at tip");
//            std::cout << "v_L_new: " << v_L_new[0] << " " << v_L_new[1] << " " << v_L_new[2] << std::endl;
//            std::cout << "w_L_new: " << w_L_new[0] << " " << w_L_new[1] << " " << w_L_new[2] << std::endl;
//            std::cout << "h0_pre_new: " << h0_pre_new[0]  << std::endl;
//            std::cout << " Residual force in core: " << residual[0] << " " << residual[1] << " " << residual[2] << std::endl;
//            std::cout << " Residual moment in core: " << residual[3] << " " << residual[4] << " " << residual[5] << std::endl;
//
//            printMatrix(&(out_calc[0]), 1, 3, "Calculated curvature at base");
//            printMatrix(&(out_calc[3]), 1, 3, "Calculated internal force at base");
//
//            for (int j = 0; j < NUM_FLEX_SEG; ++j) { // first segment
//                u_history[j].length = u_history_new[j].length;
//                u_history[j].data = new double[ u_history[j].length ];
//                for (int k = 0; k < u_history[j].length; ++k) {
//                    u_history[j].data[k] = u_history_new[j].data[k];
//                }
//
//                v_history[j].length = v_history_new[j].length;
//                v_history[j].data = new double[ v_history[j].length ];
//                for (int k = 0; k < v_history[j].length; ++k) {
//                    v_history[j].data[k] = v_history_new[j].data[k];
//                }
//
//                w_history[j].length = w_history_new[j].length;
//                w_history[j].data = new double[ w_history[j].length ];
//                for (int k = 0; k < w_history[j].length; ++k) {
//                    w_history[j].data[k] = w_history_new[j].data[k];
//                }
//
//                h0[j] = h0_pre_new[j];
//            }
//
//            //accumulated small errors
//            if (vNormSq<3>(v_L_new) > EPS){
//                mCopy_AB<3>(v_L_new, v_L_pre);
//            }else{ v_L_pre[0] = v_L_pre[1] = v_L_pre[2] = 0.0;}
//
//            if (vNormSq<3>(w_L_new) > EPS){
//                mCopy_AB<3>(w_L_new, w_L_pre);
//            }else{ w_L_pre[0] = w_L_pre[1] = w_L_pre[2] = 0.0;}
//
//            // prepare for the next round
//            auto & initial_guess = DynamicsParams.initial_guess;
//            mCopy_AB<NUM_RESIDUAL>(out_calc, initial_guess);
//
////            for (int k = 0; k < NUM_RESIDUAL; ++k) {
////                std::cout << "out_calc initial guess: " << DynamicsParams.initial_guess[k] << std::endl;
////            }
//
//
//        }
//
//        output_tip_pts[i].resize(3);
//        for (int j = 0; j < 3; ++j) {
//            output_tip_pts[i][j] = x_t_tip[j+12];
//        }
//
//
//    }
//
//
//    // Print outputs
////    std::cout << "localmin  " << localmin << std::endl;
//    printMatrix(&(out_calc[0]), 1, 3, "Calculated curvature at base");
//    printMatrix(&(out_calc[3]), 1, 3, "Calculated internal force at base");
//    if(NUM_RESIDUAL> 6) {
//        printMatrix(&(out_calc[6]), 1, 3, "Calculated internal force at upper side coil");
//        printMatrix(&(out_calc[9]), 1, 3, "Calculated internal torque at upper side coil");
//    }
////    printMatrix(ftip_calc, 1, 3, "Calculated Tip Force");
//    printMatrix(&(x_t_tip[12]), 1, 3, "Catheter Tip Position");
//    printMatrix(&(x_t_tip[3]), 1, 9, "Catheter Tip Orientation");
//
//    printMatrix(&(x_t_tip[15]), 1, 3, "v at tip");
//    printMatrix(&(x_t_tip[18]), 1, 3, "w at tip");
//
//    printMatrix(v_L_new, 1, 3, "v_L at coil");
//    printMatrix(w_L_new, 1, 3, "w_L at coil");
//    std::cout << " Residual force in core: " << residual[0] << " " << residual[1] << " " << residual[2] << std::endl;
//    std::cout << " Residual moment in core: " << residual[3] << " " << residual[4] << " " << residual[5] << std::endl;
//////    for (int i = 0; i < 1; ++i) { //check the dynamic length
//////        for (int j = 0; j < u_history_new[i].length / 3; ++j) {
//////            std::cout << j << " th u_history update: " << u_history_new[i].data[j*3] <<" " << u_history_new[i].data[j*3+1] <<" " << u_history_new[i].data[j*3+2] << std::endl;
//////        }
//////    }
//////    std::cout << "updated u_history length: " << u_history_new[0].length << std::endl;
////
//////    printMatrix(residual, 1, 3, "Force residual at coil");
//////    printMatrix(&(residual[3]), 1, 3, "Moment residual at coil");
////
////
////    for (int i = 0; i < NUM_LOCALIZATION_MARKERS; ++i) {
////        std::cout << "out_ReportedMarkerPos: " << ReportedMarkerPos[i][0] << " " << ReportedMarkerPos[i][1] << " " << ReportedMarkerPos[i][2] << std::endl;
////    }
////
//    std::cout << "----" << std::endl;
//
//    return 0;
//}
