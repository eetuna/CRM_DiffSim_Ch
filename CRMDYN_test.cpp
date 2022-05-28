#include "CRMDYN.hpp"

#include <fstream>

using namespace std;
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

    // set up the test problem inputs
    double IntegrationStepSize = 0.2;
    double InsertedLength = 104.02f;
    double B0[3] = { 0.0f, 3.0f, 0.0f };			// B0 field vector of the MRI scanner (in spatial coordinates)
    double g[3] = { 0.0, 0.0, 9.81 };
    double x0[NUM_STATES] = { 0.0f,0.0f,0.0f,   1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f,   0.0f, 0.0f, 0.0f,   0.0f,0.0f,0.0f,   0.0f,0.0f,0.0f,};
    int FinalStepOnly = false;
    const int Na = NUM_ACT_SET; 				// Number of actuator sets
    const int Nc = NUM_FLEX_SEG;				// Number of flexible catheter segments
    const int Ns = NUM_SEGMENTS;  				// NumSegments: Ns = 2 * Na + 1   // derived quantity
    const int Nm = NUM_LOCALIZATION_MARKERS;	// Number of Markers
    //	For all parameters below, segments and actuator units are numbered/ordered from the tip of the catheter towards the base
    // Array of lambda values for segment endpoints; Ns long array
    double SegEnds[Ns] = { 10.72f, 26.86f, 42.255f};//, 57.865f, 104.02f };
    // Array of lambda values for marker locations
    double MarkerLoc[Nm] = { 2.18f, 18.79f, 50.06f, 101.52f, 104.02f };
    // Catheter Rigidity Matrix; (Na+1)*9 long array, (Na+1) 3x3 matrices stored in row major order
    double Klist[Nc][9] = { {22.8304161859862f, 0.0f, 0.0f, 0.0f, 22.8304161859862f, 0.0f, 0.0f, 0.0f, 20.2125442625320f} };
    // Inverses of K matrices; (Na+1)*9 long array, (Na+1) 3x3 matrices stored in row major order  // derived quantity
    double Kinvlist[Nc][9] = { {0.0438012164059375f, 0.0f, 0.0f, 0.0f, 0.0438012164059375f, 0.0f, 0.0f, 0.0f, 0.0494742268470230f} };
    // Local curvature in unloaded configuration for each of the flexible segments; (Na+1)*3 long array, (Na+1) 3x1 vectors
    double ustarlist[Nc][3] = { 0.000127537459505030f, -0.000173238845806410f, 0.0f, 0.00566574659169662f, -0.00822266255995812f, 0.0f }; //{-0.0025f, 0.0019f, 0.0f, 0.0057f, -0.0082f, 0.0f, 0.00013f, -0.00017f, 0.0f};
    // Actuator magnetization moments; Na*3 long array, Na 3x1 vectors; MagMoment = CoilAlignMat * CoilTurnAreaMat * ActuationCurrentVector
    double MagMomentlist[Na][3] = { 0.0f, 0.0f, 0.1600f };
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

    // declare the output variables
    double xf[NUM_STATES];
    double residual[6];
    double ReportedMarkerPos[NUM_LOCALIZATION_MARKERS][3];

    double v_L_pre[3] = { 0.0, 0.0, 0.0 };
    double w_L_pre[3] = { 0.0, 0.0, 0.0 };
    double ftip[3] = {0.0, 0.0, 0.0};

    // Cosserat Rod Model - Integrator for Solving the Initial Value Problem
    CRMSolverIVP(x0, IntegrationStepSize, InsertedLength, dlambdainv, SegEnds, MarkerLoc, Klist, Kinvlist, ustarlist, v_L_pre, w_L_pre,
                 ActMass, ActInertia, MagMomentlist, fcumIVP, ftip, B0, g, FinalStepOnly, xf, residual, ReportedMarkerPos);
}


int main(int argc, char** argv){

    ivp_test();


}