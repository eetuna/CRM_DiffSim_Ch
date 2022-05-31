#include "CRMDYN.hpp"

#include <fstream>

using namespace std;

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



void ivp_test(){

    cout << "### CRMSolverIVP Test: " << endl;

//    // *** Physical Description of the Catheter
//    // IMPORTANT NOTE: For now, most proximal segment is assumed to be always flexible
//    //    and the flexible and rigid segments are assumed to be alternating
//    //    most distal segment can be flexible or rigid
//    // Items are listed in distal-to-proximal order (starting from the tip of the catheter towards the base)
//    // Outer radii of each of the flexible segments - unit: mm
    double oRlist[NUM_FLEX_SEG] = { 1.5875, 1.5875 };
    // Inner radii of each of the flexible segments - unit: mm
    double iRlist[NUM_FLEX_SEG] = { 0.9906, 0.9906 };
//    // Young's moduli of each of the flexible segments - unit: ??
    double YoungModlist[NUM_FLEX_SEG] = { 5.3948, 5.3948 };
//    // Shear moduli of each of the flexbile segments - unit: ??
    double ShearModlist[NUM_FLEX_SEG] = { 2.3881, 2.3881 };
//    // Alignment angles for the side coils, for each of the actuator sets - unit: radians
//    double CoilAlignmentAngles[NUM_ACT_SET][2] = { {0.0, 0.0}};
//    // Coil turn area matrices for each of the actuator sets  - unit: mm2
//    double CoilTurnAreaMat[NUM_ACT_SET][9] = { { 1.44, 0.0, 0.0, 0.0, 1.3851, 0.0, 0.0, 0.0, 1.60 }  };
//    // Lengths of each of the catheter segments - unit: mm
    double SegmentLengths[NUM_SEGMENTS] = { 10.72, 16.14 };
//    // Mass of each of the actuator sets - unit: ??
    double ActMass[NUM_ACT_SET] = { 7.7736e-5 };

//    // set up the test problem inputs
//    double IntegrationStepSize = 0.2;
//    double InsertedLength = 104.02f;
//    double B0[3] = { 0.0f, 3.0f, 0.0f };			// B0 field vector of the MRI scanner (in spatial coordinates)
//    double g[3] = { 0.0, 0.0, 9.81 };
//    double x0[NUM_STATES] = { 0.0f,0.0f,0.0f,   1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f,   0.0f, 0.0f, 0.0f,   0.0f,0.0f,0.0f,   0.0f,0.0f,0.0f,};
//    int FinalStepOnly = false;
//    const int Na = NUM_ACT_SET; 				// Number of actuator sets
//    const int Nc = NUM_FLEX_SEG;				// Number of flexible catheter segments
//    const int Ns = NUM_SEGMENTS;  				// NumSegments: Ns = 2 * Na + 1   // derived quantity
//    const int Nm = NUM_LOCALIZATION_MARKERS;	// Number of Markers
//    //	For all parameters below, segments and actuator units are numbered/ordered from the tip of the catheter towards the base
//    // Array of lambda values for segment endpoints; Ns long array
//    double SegEnds[Ns] = { 10.72f, 26.86f, 42.255f};//, 57.865f, 104.02f };
//    // Array of lambda values for marker locations
//    double MarkerLoc[Nm] = { 2.18f, 18.79f, 50.06f, 101.52f, 104.02f };
//    // Catheter Rigidity Matrix; (Na+1)*9 long array, (Na+1) 3x3 matrices stored in row major order
//    double Klist[Nc][9] = { {22.8304161859862f, 0.0f, 0.0f, 0.0f, 22.8304161859862f, 0.0f, 0.0f, 0.0f, 20.2125442625320f} };
//    // Inverses of K matrices; (Na+1)*9 long array, (Na+1) 3x3 matrices stored in row major order  // derived quantity
//    double Kinvlist[Nc][9] = { {0.0438012164059375f, 0.0f, 0.0f, 0.0f, 0.0438012164059375f, 0.0f, 0.0f, 0.0f, 0.0494742268470230f} };
//    // Local curvature in unloaded configuration for each of the flexible segments; (Na+1)*3 long array, (Na+1) 3x1 vectors
//    double ustarlist[Nc][3] = { 0.000127537459505030f, -0.000173238845806410f, 0.0f, 0.00566574659169662f, -0.00822266255995812f, 0.0f }; //{-0.0025f, 0.0019f, 0.0f, 0.0057f, -0.0082f, 0.0f, 0.00013f, -0.00017f, 0.0f};
//    // Actuator magnetization moments; Na*3 long array, Na 3x1 vectors; MagMoment = CoilAlignMat * CoilTurnAreaMat * ActuationCurrentVector
//    double MagMomentlist[Na][3] = { 0.0f, 0.0f, 0.1600f };
//    // Reciprocal of \Delta s (= \Delta \lambda) used in discretizing fcum  ( deltasinv = 1 / (L/N) = N/L)
//    double dlambdainv = 1.0f;
//    double FullLength = SegEnds[Ns - 1];
//    // 3*(N+1) by 1 array (grouped by 3 floats) storing cumulative external force integrated from \lambda = index * \Delta\lambda to the catheter tip (\lambda=0)
//    double fcumIVP[NUM_FCUM_LAMBDA + 1][3] = {};
//
//




    // set up the test problem inputs
    double IntegrationStepSize = 0.2f;
    double InsertedLength = 104.02f;
    double B0[3] = { 0.0f, 3.0f, 0.0f };			// B0 field vector of the MRI scanner (in spatial coordinates)
    double g[3] = { 0.0, 0.0, 9.81 };
    double x0[NUM_STATES] = { 0.0f,0.0f,0.0f,   1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f,   0.0f, 0.0f, 0.0f };
    int FinalStepOnly = false;
    const int Na = NUM_ACT_SET; 				// Number of actuator sets
    const int Nc = NUM_FLEX_SEG;				// Number of flexible catheter segments
    const int Ns = NUM_SEGMENTS;  				// NumSegments: Ns = 2 * Na + 1   // derived quantity
    const int Nm = NUM_LOCALIZATION_MARKERS;	// Number of Markers
    //	For all parameters below, segments and actuator units are numbered/ordered from the tip of the catheter towards the base
    // Array of lambda values for segment endpoints; Ns long array
    double SegEnds[Ns] = { 10.72f, 26.86f, 42.255f, 57.865f, 104.02f };
    // Array of lambda values for marker locations
    double MarkerLoc[Nm] = { 2.18f, 18.79f, 50.06f, 101.52f, 104.02f };
    // Catheter Rigidity Matrix; (Na+1)*9 long array, (Na+1) 3x3 matrices stored in row major order
    double Klist[Nc][9] = { {22.8304161859862f, 0.0f, 0.0f, 0.0f, 22.8304161859862f, 0.0f, 0.0f, 0.0f, 20.2125442625320f}, {22.8304161859862f, 0.0f, 0.0f, 0.0f, 22.8304161859862f, 0.0f, 0.0f, 0.0f, 20.2125442625320f}, {22.8304161859862f, 0.0f, 0.0f, 0.0f, 22.8304161859862f, 0.0f, 0.0f, 0.0f, 20.2125442625320f} };
    // Inverses of K matrices; (Na+1)*9 long array, (Na+1) 3x3 matrices stored in row major order  // derived quantity
    double Kinvlist[Nc][9] = { {0.0438012164059375f, 0.0f, 0.0f, 0.0f, 0.0438012164059375f, 0.0f, 0.0f, 0.0f, 0.0494742268470230f}, {0.0438012164059375f, 0.0f, 0.0f, 0.0f, 0.0438012164059375f, 0.0f, 0.0f, 0.0f, 0.0494742268470230f}, {0.0438012164059375f, 0.0f, 0.0f, 0.0f, 0.0438012164059375f, 0.0f, 0.0f, 0.0f, 0.0494742268470230f} };
    // Local curvature in unloaded configuration for each of the flexible segments; (Na+1)*3 long array, (Na+1) 3x1 vectors
    double ustarlist[Nc][3] = { 0.000127537459505030f, -0.000173238845806410f, 0.0f, 0.00566574659169662f, -0.00822266255995812f, 0.0f, -0.00250030480901699f, 0.00188434148765547f, 0.0f }; //{-0.0025f, 0.0019f, 0.0f, 0.0057f, -0.0082f, 0.0f, 0.00013f, -0.00017f, 0.0f};
    // Actuator magnetization moments; Na*3 long array, Na 3x1 vectors; MagMoment = CoilAlignMat * CoilTurnAreaMat * ActuationCurrentVector
    double MagMomentlist[Na][3] = { 0.0f, 0.0f, 0.1600f, 0.0f, 0.0f, 0.0f };
    // Reciprocal of \Delta s (= \Delta \lambda) used in discretizing fcum  ( deltasinv = 1 / (L/N) = N/L)
    double dlambdainv = 1.0f;
    double FullLength = SegEnds[Ns - 1];
    // 3*(N+1) by 1 array (grouped by 3 floats) storing cumulative external force integrated from \lambda = index * \Delta\lambda to the catheter tip (\lambda=0)
    double fcumIVP[NUM_FCUM_LAMBDA + 1][3] = {};





    double ActInertia[NUM_ACT_SET][9];
    for (int i = 0; i < NUM_ACT_SET; ++i)
    {
        double I_zz = 0.5 * ActMass[i] * (oRlist[0] * oRlist[0] + iRlist[0] * iRlist[0]);
        double I_xx = 0.25 * ActMass[i] * (oRlist[0] * oRlist[0] + iRlist[0] * iRlist[0]) + 1 / 12 * ActMass[i] * SegmentLengths[2*i+1];
        ActInertia[i][0] = I_xx; ActInertia[i][1] = 0.0; ActInertia[i][2] = 0.0;
        ActInertia[i][3] = 0.0; ActInertia[i][4] = I_xx; ActInertia[i][5] = 0.0;
        ActInertia[i][6] = 0.0; ActInertia[i][7] = 0.0; ActInertia[i][8] = I_zz;
    }
    double v_L_pre[3] = { 0.0, 0.0, 0.0 };
    double w_L_pre[3] = { 0.0, 0.0, 0.0 };
    double ftip[3] = {0.0, 0.0, 0.0};

    // declare the output variables
    double xf[NUM_STATES];
    double residual[6];
    double ReportedMarkerPos[NUM_LOCALIZATION_MARKERS][3];



    // Cosserat Rod Model - Integrator for Solving the Initial Value Problem
    CRMSolverIVP(x0, IntegrationStepSize, InsertedLength, dlambdainv, SegEnds, MarkerLoc, Klist, Kinvlist, ustarlist, v_L_pre, w_L_pre,
                 ActMass, ActInertia, MagMomentlist, fcumIVP, ftip, B0, g, FinalStepOnly, xf, residual, ReportedMarkerPos);
}


int RunTests (void) {

    int FailFlag = 0;

    std::cout << "######## CRM Tools Test Results " << std::endl;

    //
    //  Test code for CRMSolverIVP  (Solver for Initial Value Problem)
    //
    std::cout << "### CRMSolverIVP Test: " << std::endl;

    double oRlist[NUM_FLEX_SEG] = { 1.5875, 1.5875 };
    // Inner radii of each of the flexible segments - unit: mm
    double iRlist[NUM_FLEX_SEG] = { 0.9906, 0.9906 };

    // set up the test problem inputs
    double IntegrationStepSize = 0.2;
    double InsertedLength = 104.02;
    double B0[3] = { 0.0, 3.0, 0.0 };			// B0 field vector of the MRI scanner (in spatial coordinates)
    double g[3] = { 0.0, 0.0, 9.81 };
    double x0[NUM_STATES] = { 0.0,0.0,0.0,   1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,   0.0, 0.0, 0.0 };
    bool FinalStepOnly = false;
    const int Na = NUM_ACT_SET; 				// Number of actuator sets
    const int Nc = NUM_FLEX_SEG;				// Number of flexible catheter segments
    const int Ns = NUM_SEGMENTS;  				// NumSegments: Ns = 2 * Na + 1   // derived quantity
    const int Nm = NUM_LOCALIZATION_MARKERS;	// Number of Markers
    //	For all parameters below, segments and actuator units are numbered/ordered from the tip of the catheter towards the base
    // Array of lambda values for segment endpoints; Ns long array
    double SegEnds[Ns] = { 10.72, 26.86, 42.255, 57.865, 104.02 };
    // Array of lambda values for marker locations
    double MarkerLoc[Nm] = { 2.18, 18.79, 50.06, 101.52, 104.02 };
    // Catheter Rigidity Matrix; (Na+1)*9 long array, (Na+1) 3x3 matrices stored in row major order
    double Klist[Nc][9] = { {22.8304161859862, 0.0, 0.0, 0.0, 22.8304161859862, 0.0, 0.0, 0.0, 20.2125442625320}, {22.8304161859862, 0.0, 0.0, 0.0, 22.8304161859862, 0.0, 0.0, 0.0, 20.2125442625320}, {22.8304161859862, 0.0, 0.0, 0.0, 22.8304161859862, 0.0, 0.0, 0.0, 20.2125442625320} };
    // Inverses of K matrices; (Na+1)*9 long array, (Na+1) 3x3 matrices stored in row major order  // derived quantity
    double Kinvlist[Nc][9] = { {0.0438012164059375, 0.0, 0.0, 0.0, 0.0438012164059375, 0.0, 0.0, 0.0, 0.0494742268470230}, {0.0438012164059375, 0.0, 0.0, 0.0, 0.0438012164059375, 0.0, 0.0, 0.0, 0.0494742268470230}, {0.0438012164059375, 0.0, 0.0, 0.0, 0.0438012164059375, 0.0, 0.0, 0.0, 0.0494742268470230} };
    // Local curvature in unloaded configuration for each of the flexible segments; (Na+1)*3 long array, (Na+1) 3x1 vectors
    double ustarlist[Nc][3] = { 0.000127537459505030, -0.000173238845806410, 0.0, 0.00566574659169662, -0.00822266255995812, 0.0, -0.00250030480901699, 0.00188434148765547, 0.0 }; //{-0.0025, 0.0019, 0.0, 0.0057, -0.0082, 0.0, 0.00013, -0.00017, 0.0};
    // Actuator magnetization moments; Na*3 long array, Na 3x1 vectors; MagMoment = CoilAlignMat * CoilTurnAreaMat * ActuationCurrentVector
    double MagMomentlist[Na][3] = { 0.0, 0.0, 0.1600, 0.0, 0.0, 0.0 };
    // Reciprocal of \Delta s (= \Delta \lambda) used in discretizing fcum  ( deltasinv = 1 / (L/N) = N/L)
    double dlambdainv = 1.0;
    double FullLength = SegEnds[Ns - 1];
    // 3*(N+1) by 1 array (grouped by 3 doubles) storing cumulative external force integrated from \lambda = index * \Delta\lambda to the catheter tip (\lambda=0)
    double fcumIVP[NUM_FCUM_LAMBDA + 1][3] = {};
    double ftip[3] = { 0.0, 0.0, 0.0 };

    double ActMass[NUM_ACT_SET] = { 7.7736e-5 };

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

    // declare the output variables
    double xf[NUM_STATES];
    double residual[3];
    double ReportedMarkerPos[NUM_LOCALIZATION_MARKERS][3];

    // Cosserat Rod Model - Integrator for Solving the Initial Value Problem
    CRMSolverIVP<double>(x0, IntegrationStepSize, InsertedLength, dlambdainv, SegEnds, MarkerLoc, Klist, Kinvlist, ustarlist, v_L_pre, w_L_pre,
                         ActMass, ActInertia, MagMomentlist, fcumIVP, ftip, B0, g, FinalStepOnly, xf, residual, ReportedMarkerPos);
    // compare the results to the reference values
    double ExpectedStates[NUM_STATES] = { 0.0225, -0.0024, -0.0074, 0.982564773234176, 0.029114762388711, -0.183626784558004, -0.094107711560958, 0.929675871807347, -0.356155178150316, 0.160344017817671, 0.367226228689942, 0.916206686783898, -5.497054669396932, -5.570783408417309, 103.0917629459139 }; //  ???f, ???f, ???f,  ???f, ???f, ???f, ???f, ???f, ???f, ???f, ???f, ???f,   ???f, ???f, ???f };
    double ExpectedMarkerPos[NUM_LOCALIZATION_MARKERS][3] = { {-5.1038, -4.8439, 101.0747}, {-2.4420, -1.9709, 84.9666}, {0.0000, 0.0000, 53.9600}, {0.0000, 0.0000, 2.5000}, {0, 0, 0} };
    double errorv2[NUM_STATES], ssqerrorv2;
    double errorv3[NUM_LOCALIZATION_MARKERS][3], ssqerrorv3;

    ssqerrorv2 = 0.0;
    for (int i = 0; i < NUM_STATES; i++) ssqerrorv2 += pow((errorv2[i] = ExpectedStates[i] - xf[i]), 2);
    ssqerrorv3 = 0.0;
    for (int i = 0; i < NUM_LOCALIZATION_MARKERS; i++)
        for (int j = 0; j < 3; j++)
            ssqerrorv3 += pow((errorv3[i][j] = ExpectedMarkerPos[i][j] - ReportedMarkerPos[i][j]), 2);

    printMatrix(ExpectedStates, 1, NUM_STATES, "ExpectedStates");
    printMatrix(xf, 1, NUM_STATES, "xf");
//    printMatrix(errorv2, 1, NUM_STATES, "errorv2");
//    std::cout << "Error Norm:" << sqrt(ssqerrorv2) << std::endl;
//    std::cout << "----" << std::endl;
//    printMatrix(&(ExpectedMarkerPos[0][0]), NUM_LOCALIZATION_MARKERS, 3, "ExpectedMarkerPos");
//    printMatrix(&(ReportedMarkerPos[0][0]), NUM_LOCALIZATION_MARKERS, 3, "ReportedMarkerPos");
//    printMatrix(&(errorv3[0][0]), NUM_LOCALIZATION_MARKERS, 3, "errorv3");
//    std::cout << "Error Norm:" << sqrt(ssqerrorv3) << std::endl;
//    std::cout << "----" << std::endl;
//
//    if ((sqrt(ssqerrorv2) < 0.2) && (sqrt(ssqerrorv3) < 0.2)) {
//        std::cout << "Test passed for Inserted Length =" << InsertedLength << std::endl;
//        std::cout << "###" << std::endl;
//    }
//    else {
//        std::cout << "TEST FAIL: " << "Results do not match!" << std::endl;
//        std::cout << "###" << std::endl;
//        FailFlag++;
//    }

}

int main(int argc, char** argv){

//    ivp_test();
    RunTests();

}