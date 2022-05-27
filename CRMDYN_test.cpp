#include "CRMDYN.hpp"

#include <fstream>

using namespace std;
void ivp_test(){
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
//
//    double ActInertia[NUM_ACT_SET][9];
//    for (int i = 0; i < NUM_ACT_SET; ++i)
//    {
//        double I_zz = 0.5 * ActMass[i] * (oRlist[0] * oRlist[0] + iRlist[0] * iRlist[0]);
//        double I_xx = 0.25 * ActMass[i] * (oRlist[0] * oRlist[0] + iRlist[0] * iRlist[0]) + 1 / 12 * ActMass[i] * SegmentLengths[2*i+1];
//        ActInertia[i][0] = I_xx; ActInertia[i][1] = 0.0; ActInertia[i][2] = 0.0;
//        ActInertia[i][3] = 0.0; ActInertia[i][4] = I_xx; ActInertia[i][5] = 0.0;
//        ActInertia[i][6] = 0.0; ActInertia[i][7] = 0.0; ActInertia[i][8] = I_zz;
//    }

//    // Array of lambda values for marker locations  (distance from the tip to each of the markers) - unit: mm
//    double MarkerLoc[NUM_LOCALIZATION_MARKERS] = { 2.18, 18.79, 50.06, 101.52, 104.02 };
//    // Length density for each of the catheter segments
//    double rho[NUM_SEGMENTS] = { 2.8814e-7, 2.8814e-7, 2.8814e-7};
//    // Local curvature in unloaded configuration for each of the flexible segments; NUM_FLEX_SEG*3 long array corresponding to NUM_FLEX_SEG many 3x1 vectors
//    double ustarlist[NUM_FLEX_SEG][3] = { 0.000128, -0.000173, 0.0, 0.00567, -0.00822, 0.0};
//
//    // *** Catheter Configuration in spatial coordinates
//    // B0 field vector of the MRI scanner (in spatial coordinates) - unit: Tesla
//    double B0[3] = { 0.0, 3.0, 0.0 };
//    // Gravity vector - unit: ??
//    double g[3] = { 0.0, 0.0, 9.81 };
//    // Catheter entry point coordinates (in spatial coordinates) - unit: mm
//    double p0[3] = { 0.0, 0.0, 0.0 };
//    // Catheter entry point orientation (3x3 rotation matrix describing catheter entry point frame orientation relative to the spatial frame stored in row major order) - unit: unitless
//    double R0[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
//
//    // *** Other External variables
//    // specify if catheter is in free space or if the catheter tip is constrained to a contact point
////    ContactModeType ContactMode = ContactModeType::FREE_TIP;
//    // External point force (in spatial coordinates) applied at the tip of the catheter (\lambda = 0)  - unit: ??
//    //   (this will be used when ContactMode == ContactModeType::FREE_TIP)
//    double TipForce[3] = { 0.0, 0.0, 0.0 };
//    // The spatial coordinates of the point where the catheter tip is constrained to be
//    //   (this will be used when ContactMode == ContactModeType::FIXED_TIP)
//    double TipConstraintPoint[3] = { 0.0, 0.0, 0.0 };
//
//    // *** Control Inputs
//    // Inserted Length of the catheter (length of the catheter that is inside the heart chamber) - unit: mm
//    double InsertedLength = 104.0;
//    // Actuation currents for each of the coils for each of the coil sets - unit: A
//    double ActuationCurrents[NUM_ACT_SET][3] = { {0.100, 0.100, 0.100} };
//
//    // *** Numerical Computation Parameters
//    // Stepsize used in numerical integration along the length of the catehter during IVP - unit: mm
//    double IntegrationStepSize = 0.2;
//
//
//    double u0_initialguess[3] = { 0.0, 0.0, 0.0 };
//    double v0_initialguess[3] = { 0.0, 0.0, 0.0 };
//    double w0_initialguess[3] = { 0.0, 0.0, 0.0 };
//    double x_0[NUM_STATES];
//    for (int i = 0; i < NUM_STATES; i++) {
//        if (i < 3) {
//            x_0[i] = u0_initialguess[i];
//        }
//        else if (i < 12) {
//            x_0[i] = R0[i - 3];
//        }
//        else if( i< 15) {
//            x_0[i] = p0[i - 12];
//        }
//    }
//    for (int i = 0; i < 3; ++i) {
//        x_0[i+3+9+3] = v0_initialguess[i];
//        x_0[i+3+9+6] = w0_initialguess[i];
//
//    }
//
//    double SegEndLambdas[NUM_SEGMENTS], LocMarkerLambdas[NUM_LOCALIZATION_MARKERS];
//    for (int i = 0; i < NUM_SEGMENTS; ++i) {
//        SegEndLambdas[i] = 1.0; //test numbers
//    }
//    for (int i = 0; i < NUM_LOCALIZATION_MARKERS; ++i) {
//        LocMarkerLambdas[i] = 20.0; // test numbers
//    }
//
//    double oR, iR, mI, pmI, E, G, K[9], Kinv[9], in_K[NUM_FLEX_SEG][9], in_Kinv[NUM_FLEX_SEG][9], in_u_star[NUM_FLEX_SEG][3];
//    for (int i = 0; i < NUM_FLEX_SEG; i++) {
//        oR = oRlist[i];
//        iR = iRlist[i];
//        E = YoungModlist[i];
//        G = ShearModlist[i];
//        mI = 0.25 * M_PI * (POW4(oR) - POW4(iR));	// Area moment of inertia along x - axis
//        pmI = 0.5 * M_PI * (POW4(oR) - POW4(iR));	// Polar moment of inertia of area
//        K[0] = E * mI;		K[1] = 0.0;			K[2] = 0.0;
//        K[3] = 0.0;			K[4] = E * mI;		K[5] = 0.0;
//        K[6] = 0.0;			K[7] = 0.0;			K[8] = G * pmI;
//        Kinv[0] = 1.0 / (E * mI);		Kinv[1] = 0.0;				Kinv[2] = 0.0;
//        Kinv[3] = 0.0;				Kinv[4] = 1.0 / (E * mI);		Kinv[5] = 0.0;
//        Kinv[6] = 0.0;				Kinv[7] = 0.0;				Kinv[8] = 1.0 / (G * pmI);
//        mCopy_AB<9>(K, in_K[i]);
//        mCopy_AB<9>(Kinv, in_Kinv[i]);
//        mCopy_AB<3>(u0_initialguess, in_u_star[i]);
//    }
//
//    double CoilAlignMat[9], c0, s0, c1, s1;
//    double tempf[3], in_MagMoment[NUM_ACT_SET][3];
//    for (int i = 0; i < NUM_ACT_SET; i++) {
//        c0 = cos(CoilAlignmentAngles[i][0]);
//        s0 = sin(CoilAlignmentAngles[i][0]);
//        c1 = cos(CoilAlignmentAngles[i][1]);
//        s1 = sin(CoilAlignmentAngles[i][1]);
//        mMult_AB<3, 3, 1>(CoilTurnAreaMat[i], ActuationCurrents[i], tempf);
//        CoilAlignMat[0] = c0;	CoilAlignMat[1] = -s1;	CoilAlignMat[2] = 0.0;
//        CoilAlignMat[3] = s0;	CoilAlignMat[4] = c1;	CoilAlignMat[5] = 0.0;
//        CoilAlignMat[6] = 0.0;	CoilAlignMat[7] = 0.0;	CoilAlignMat[8] = 1.0;
//        mMult_AB<3, 3, 1>(CoilAlignMat, tempf, in_MagMoment[i]);
//    }
//    double in_fcumlambda[NUM_FCUM_LAMBDA+1][3];
//    in_fcumlambda[0][0] = 0.0;
//    in_fcumlambda[0][1] = 0.0;
//    in_fcumlambda[0][2] = 0.0;
//
//
//    double Li =100.0;
//    double dlambdainv = 1 / 100.0;
//
//    //u_history
//
//    CRMIVPCoreParams<double> testParams;
//    CRMSolverIVP_Prep ( x_0, IntegrationStepSize,
//            Li, dlambdainv,
//            SegEndLambdas, LocMarkerLambdas,
//            in_K, in_Kinv, in_u_star,
//            in_MagMoment, in_fcumlambda,
//            B0, g,
//            true, testParams);
//
//







    cout << "### CRMSolverIVP Test: " << endl;

    // set up the test problem inputs
    float IntegrationStepSize = 0.2f;
    float InsertedLength = 104.02f;
    float B0[3] = { 0.0f, 3.0f, 0.0f };			// B0 field vector of the MRI scanner (in spatial coordinates)
    double g[3] = { 0.0, 0.0, 9.81 };
    float x0[NUM_STATES] = { 0.0f,0.0f,0.0f,   1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f,   0.0f, 0.0f, 0.0f };
    int FinalStepOnly = false;
    const int Na = NUM_ACT_SET; 				// Number of actuator sets
    const int Nc = NUM_FLEX_SEG;				// Number of flexible catheter segments
    const int Ns = NUM_SEGMENTS;  				// NumSegments: Ns = 2 * Na + 1   // derived quantity
    const int Nm = NUM_LOCALIZATION_MARKERS;	// Number of Markers
    //	For all parameters below, segments and actuator units are numbered/ordered from the tip of the catheter towards the base
    // Array of lambda values for segment endpoints; Ns long array
    float SegEnds[Ns] = { 10.72f, 26.86f, 42.255f};//, 57.865f, 104.02f };
    // Array of lambda values for marker locations
    float MarkerLoc[Nm] = { 2.18f, 18.79f, 50.06f, 101.52f, 104.02f };
    // Catheter Rigidity Matrix; (Na+1)*9 long array, (Na+1) 3x3 matrices stored in row major order
    float Klist[Nc][9] = { {22.8304161859862f, 0.0f, 0.0f, 0.0f, 22.8304161859862f, 0.0f, 0.0f, 0.0f, 20.2125442625320f} };
    // Inverses of K matrices; (Na+1)*9 long array, (Na+1) 3x3 matrices stored in row major order  // derived quantity
    float Kinvlist[Nc][9] = { {0.0438012164059375f, 0.0f, 0.0f, 0.0f, 0.0438012164059375f, 0.0f, 0.0f, 0.0f, 0.0494742268470230f} };
    // Local curvature in unloaded configuration for each of the flexible segments; (Na+1)*3 long array, (Na+1) 3x1 vectors
    float ustarlist[Nc][3] = { 0.000127537459505030f, -0.000173238845806410f, 0.0f, 0.00566574659169662f, -0.00822266255995812f, 0.0f }; //{-0.0025f, 0.0019f, 0.0f, 0.0057f, -0.0082f, 0.0f, 0.00013f, -0.00017f, 0.0f};
    // Actuator magnetization moments; Na*3 long array, Na 3x1 vectors; MagMoment = CoilAlignMat * CoilTurnAreaMat * ActuationCurrentVector
    float MagMomentlist[Na][3] = { 0.0f, 0.0f, 0.1600f };
    // Reciprocal of \Delta s (= \Delta \lambda) used in discretizing fcum  ( deltasinv = 1 / (L/N) = N/L)
    float dlambdainv = 1.0f;
    float FullLength = SegEnds[Ns - 1];
    // 3*(N+1) by 1 array (grouped by 3 floats) storing cumulative external force integrated from \lambda = index * \Delta\lambda to the catheter tip (\lambda=0)
    float fcumIVP[NUM_FCUM_LAMBDA + 1][3] = {};

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
    float xf[NUM_STATES];
    float residual[3];
    float ReportedMarkerPos[NUM_LOCALIZATION_MARKERS][3];

    // Cosserat Rod Model - Integrator for Solving the Initial Value Problem
//    CRMSolverIVP(x0, IntegrationStepSize, InsertedLength, dlambdainv, SegEnds, MarkerLoc, Klist, Kinvlist, ustarlist, MagMomentlist, fcumIVP, B0, FinalStepOnly, xf, residual, ReportedMarkerPos);


}


int main(int argc, char** argv){

    ivp_test();


}