#pragma once
#include <cmath>
#include "CRMMatrixOperations.hpp"

#define NUM_STATES 21   // u[0..2],R[0..9],p[0..2],v[0..2],w[0..2]
#define NUM_INTEGRATION_STATES 9   // u[0..2],v[0..2],w[0..2]
#define EQNDIMENSION 9	// domain: u[0..2],v[0..2],w[0..2] range: m_tip[0..2]
#define DELTA_T

#define NUM_ACT_SET 1								// Number of actuator sets
#define NUM_FLEX_SEG 2								// Number of flexible segments
#define NUM_SEGMENTS (NUM_ACT_SET+NUM_FLEX_SEG)		// Total number of segments
#define NUM_LOCALIZATION_MARKERS 10					// Total number of localization markers

#define IVALUE_SCALE_U	1.0 //(0.01)			// the variable used in Nonlinear Solver is multiplied with this scale to calculate u (curvature) that will be used in IVP
#define RESIDUAL_SCALE_M	1.0 //(10.0)			// the residual for tip moment coming out of the IVP will be multiplied with this scale to return to the Nonlinear Solver
#define RESIDUAL_SCALE_P	100.0 //(10.0)			// the residual for tip position error coming out of the IVP will be multiplied with this scale to return to the Nonlinear Solver

#define DELTA_T 0.02
#define TRUSTREGION							// Trust Region Method with the numerical jacobian (Default)


// Structure for defining dynamics model parameters
struct CRMDynamicsModelParams {
    //	For all parameters below, segments and actuator units are numbered/ordered from the tip of the catheter towards the base (distal to proximal)
    double 	SegLengths[NUM_SEGMENTS];				// Array of segment lengths; NUM_SEGMENTS long array
    double	LocMarkers[NUM_LOCALIZATION_MARKERS];	// Array of localization marker locations (in lambda coordinates); NUM_LOCALIZATION_MARKERS long array
    double	InnerRadius[NUM_FLEX_SEG];				// Inner radii of the flexible catheter segments
    double	OuterRadius[NUM_FLEX_SEG];				// Outer radii of the flexible catheter segments
    double	YoungsModulus[NUM_FLEX_SEG];			// Youngs Moduli of the flexible catheter segments
    double	ShearModulus[NUM_FLEX_SEG];				// Shear Moduli of the flexible catheter segments
    double 	ustar[NUM_FLEX_SEG][3];					// Local curvature in unloaded configuration for each of the flexible segments; (NUM_FLEX_SEG)*3 long array, (NUM_FLEX_SEG) 3x1 vectors
    //		Actuator segments are assumed to be straight
    double 	CoilAlignmentAngles[NUM_ACT_SET][2];	// Coil Alignment Angles; NUM_ACT_SET*2 long array, for each actuator set, the angle for the first coil is relative to x axis, and the angle for the second coil is relative to y axis
    double 	CoilTurnAreaMat[NUM_ACT_SET][9]; 		// Coil Turn Area matrices; NUM_ACT_SET*9 long array, NUM_ACT_SET 3x3 matrices stored in row major order
    double 	rho[NUM_SEGMENTS];						// Length density (mass per unit length) of the flexible catheter substrate (tubing); (NUM_SEGMENTS) long array
    double 	ActMass[NUM_ACT_SET];					// Actuator segment masses, does not include the flexible substrate; (NUM_ACT_SET) long array
    double delta_t
};

// Structure for defining catheter configuration parameters
struct CatheterConfiguration {
    double 	B0[3];		//  B0 field vector of the MRI scanner (in spatial coordinates)
    double 	g[3];		//  Gravity vector (in spatial coordinates)
    double 	p0[3];		//  Catheter entry port position (in spatial coordinates)
    double 	R0[9];		//	Catheter orientation at the entry port (relative to the spatial frame); 3x3 matrix stored in row major order
};

// Structure for passing parameters to the CRM Shooting Method Boundary Value Problem Solver CRMShootingMethodBVP
template <typename adType>
struct CRMDYNShootingMethodParams {
    //	For all parameters below, segments and actuator units are numbered/ordered from the tip of the catheter towards the base
    //  DISTAL TO PROXIMAL ORDERING
    adType 	Li;
    double 	dlambdainv;
    double 	SegEndLambdas[NUM_SEGMENTS];
    double	LocMarkerLambdas[NUM_LOCALIZATION_MARKERS];
    double 	K[NUM_FLEX_SEG][9];
    double 	Kinv[NUM_FLEX_SEG][9];
    double 	ustar[NUM_FLEX_SEG][3];
    adType 	MagMoment[NUM_ACT_SET][3];
    double 	fcumlambda[NUM_FCUM_LAMBDA + 1][3];
    double 	B0[3];
    double	IntegrationStepSize;
    double	R0[9];
    double	p0[3];
};

//
// ---------------------------------------------------------
//        Cosserat Rod Dynamics - Integrator for Solving the Initial Value Problem
// ---------------------------------------------------------
//

//  Parameter Set used to call CRMSolverIVP_Core
template <typename adType>
struct CRMIVPCoreParams {
    //  These are the intermediate variables that will be used to call CRMSolverIVP_Core
    //	For all parameters below, segments and actuator units are numbered/ordered from the base of the catheter towards the tip (proximal to distal)
    adType xi[NUM_STATES];  							// Initial value of the state for the next segment to be integrated
    //  States are packed u[0..2],R[0..8],p[0..2]  (R: 3x3 matrix stored in row major order R11 R12 R13 R21 R22 R23 R31 R32 R33)
    adType SegBounds[NUM_SEGMENTS+1];					// The s values as each of the segment boundaries  (For NUM_SEGMENTS segments, there are NUM_SEGMENTS+1 boundaries)
    int   SegSteps[NUM_FLEX_SEG];						// Number of integration steps in each flexible catheter segment

    adType InsertedLength;								// Inserted Length (length of the catheter from the entry point to the tip)
    double dlambdainv; 									// reciprocal of dlambda (lambda stepsize used in discretizing fcumlambda)
    double K[NUM_FLEX_SEG][9]; 							// Catheter Rigidity Matrices;
    double Kinv[NUM_FLEX_SEG][9];						// Inverses of K Matrices
    double ustar[NUM_FLEX_SEG][3];						// Local curvature in unloaded configuration for each of the flexible segments
    double fcumlambda[NUM_FCUM_LAMBDA+1][3];			// 3*(NUM_FCUM_LAMBDA+1) by 1 array (grouped by 3 doubles) storing cumulative external force (excluding tip force) integrated from \lambda = index * \Delta\lambda to the catheter tip (\lambda=0) - in spatial (catheter base frame) coordinates
    bool   FinalValueOnly;								// Flag used to indicate if only final value (xf) is returned (true) or if Marker Locations are returned as well (false)
    double LocMarkers[NUM_LOCALIZATION_MARKERS];		// The s values of each of the localization markers (markers not yet inserted into the catheter would have a negative s value

    double B0[3];										// B0 field vector of the MRI scanner (in spatial coordinates)
    adType MagMoment[NUM_ACT_SET][3]; 					// Actuator magnetization moments in body coordinates; Na*3 long array, Na 3x1 vectors; MagMoment = CoilAlignMat * CoilTurnAreaMat * ActuationCurrentVector

    int   StartSegmentIndex;							// Index of the segment where the integration to solve IVP will start -- the segment located at the entry point; note that segment indices start at 0
    int	  NextLocMarker;								// Next Localization Marker to be computed
    double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]; // positions of markers (ordered proximal to distal)
    double u_history[NUM_FLEX_SEG][SegSteps[NUM_FLEX_SEG]][3]; //History of u
    //   only the entries 0..NextLocMarker-1 are filled
};

// CRMSolverIVP API which exposes all of the individual parameters
template <typename adType>
void CRMSolverIVP (	double in_x_0[NUM_STATES], double in_IntegrationStepSize,
                       double in_Li, double in_dlambdainv,
                       double in_SegEndLambdas[NUM_SEGMENTS], double	in_LocMarkerLambdas[NUM_LOCALIZATION_MARKERS],
                       double in_K[NUM_FLEX_SEG][9], double in_Kinv[NUM_FLEX_SEG][9],
                       double in_ustar[NUM_FLEX_SEG][3], double in_vstar[NUM_FLEX_SEG][3],double in_wstar[NUM_FLEX_SEG][3],
                       double in_MagMoment[NUM_ACT_SET][3], double in_fcumlambda[NUM_FCUM_LAMBDA+1][3], double in_ftp[3],
                       double in_B0[3],
                       bool in_FinalValueOnly,
                       double out_x_N[NUM_STATES], double out_WrenchResidual[6],
                       double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]
);


// Preparation of CRMIVPCoreParams for subsequent call to CRMSolverIVP_Core
template <typename adType>
void CRMSolverIVP_Prep ( adType in_x_0[NUM_STATES], double in_IntegrationStepSize,
                         adType in_Li, double in_dlambdainv,
                         double in_SegEndLambdas[NUM_SEGMENTS], double in_LocMarkerLambdas[NUM_LOCALIZATION_MARKERS],
                         double in_K[NUM_FLEX_SEG][9], double in_Kinv[NUM_FLEX_SEG][9],
                         double in_ustar[NUM_FLEX_SEG][3],double in_vstar[NUM_FLEX_SEG][3],double in_wstar[NUM_FLEX_SEG][3],
                         adType in_MagMoment[NUM_ACT_SET][3], double in_fcumlambda[NUM_FCUM_LAMBDA+1][3],
                         double in_B0[3],
                         bool in_FinalValueOnly,
                         CRMIVPCoreParams<adType> &out_CoreParams	);


// Core Computations used in CRMSolverIVP - Integrator for Solving the Initial Value Problem
//    parameters are prepared using CRMSolverIVP_Prep
//    CRMSolverIVP_Core can be called multiple times by changing only u[0..2] components of xi (in_u[]) and in_ftip[]
//    for a single execution of CRMSolverIVP_Prep - the values in in_params should not be changed by user
template <typename adType>
void CRMSolverIVP_Core ( CRMIVPCoreParams<adType> in_params,
                         adType in_u[3], adType in_v[3], adType in_w[3], adType in_ftip[3],
                         adType out_x_N[NUM_STATES], adType out_WrenchResidual[6],
                         double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3] );


// Function for copying data in device memory to global memory --- used for dataflow pipelining
template <typename adType>
void CRMSolverIVP_Return (  adType in_x_N[NUM_STATES], adType out_WrenchResidual[6],
                            double in_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3],
                            double out_x_N[NUM_STATES], double out_WrenchResidual[6],
                            double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]);

// Cosserat Dynamics Model Integrand
template <typename adType>
void CRMIntegrand (	adType s, adType x[NUM_STATES],
                       adType Li, double dlambdainv,
                       double in_K[9], double in_Kinv[9], double in_l[3], double in_ustar[3], double in_fcumlambda[NUM_FCUM_LAMBDA+1][3],
                       adType in_ftip[3], double u_pre[3],
                       adType xdot[NUM_INTEGRATION_STATES]);
//double 	Li;				// Inserted length of the catheter (from s=0 to the tip)
//double 	dlambdainv;		// Reciprocal of \Delta \lambda (= \Delta s) used in discretizing fcum  ( dlambdainv = 1 / (Lf/NUM_FCUM_LAMBDA) = NUM_FCUM_LAMBDA/Lf ), Lf: functional (full) length of the catheter
//double 	ustar[3];		// Local curvature in unloaded configuration at current s (3x1 array)
// -UNUSED- double ustardot[3];	// Derivative of local curvature in unloaded configuration at current s; we assume ustardot=0.0 since our rest shape model is piecewise constant curvature
//double 	K[9];	  		// Catheter Rigidity Matrix at current s (3x3 matrix stored in row major order)
//double 	Kinv[9];  		// Inverse of K at current s (3x3 matrix stored in row major order)
// -UNUSED- double Kdot[9];   	// Derivative of K at current s (3x3 matrix stored in row major order);  we assume Kdot=0.0
//double 	l[3];       	// External moment at current s (world coordinates)
//double 	fcumlambda[NUM_FCUM_LAMBDA+1][3];	// (NUM_FCUM_LAMBDA+1)x3 array storing cumulative external force (exluding tip force) integrated from \lambda = index * \Delta\lambda to the catheter tip (\lambda=0)
//double		ftip[3];		// External point force (in spatial coordinates) applied at the tip of the catheter (\lambda = 0)
// Returning xdot-- EQ:(9) Rucker 2010, u_pre_: updated u as previous time value
// support function to copy location marker positions
//template <typename adType>
//void LocMarkerUpdate(double p[3], adType xi[NUM_STATES], adType t);

template <typename adType>
double dVal(adType x) { return (x); }

/*
 * Numerical Integration Functions
 */
//
//ABM4: 4th order Adams-Bashforth Prediction and Adams-Moulton Correction Numerical Integration
//    Note: The first three steps are calculated using RK2
template <typename adType>
void ABM4 (	adType in_x_0[NUM_STATES], adType t_0, int N, double h,
               adType Li, double dlambdainv, double in_K[9], double in_Kinv[9], double in_l[3], double in_ustar[3], double u_history[N][3], double in_fcumlambda[NUM_FCUM_LAMBDA+1][3], adType in_ftip[3],
               bool FinalValueOnly, double	in_LocMarkers[NUM_LOCALIZATION_MARKERS], int *inout_NextLocMarkerIdx,
               adType out_x_N[NUM_STATES], double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]	);

//ABM4_step One step of 4th order Adams-Bashforth Prediction and Adams-Moulton Correction
template <typename adType>
void ABM4_step(	adType in_x_n[NUM_STATES], adType t_n, double h,
                   adType in_xdot_nm1[NUM_INTEGRATION_STATES], adType in_xdot_nm2[NUM_INTEGRATION_STATES], adType in_xdot_nm3[NUM_INTEGRATION_STATES],
                   adType in_x_nm1[NUM_STATES], adType in_x_nm2[NUM_STATES], adType in_x_nm3[NUM_STATES],
                   adType Li, double dlambdainv, double in_K[9], double in_Kinv[9], double in_l[3], double in_ustar[3], double in_fcumlambda[NUM_FCUM_LAMBDA + 1][3], adType in_ftip[3],
                   adType out_x_np1[NUM_STATES], adType out_xdot_n[NUM_INTEGRATION_STATES]);

//RK2_step One step of 2nd Order Runge-Kutta Integration
template <typename adType>
void RK2_step(	adType in_x_n[NUM_STATES], adType t_n, double h,
                  adType Li, double dlambdainv, double in_K[9], double in_Kinv[9], double in_l[3], double in_ustar[3], double u_pre[3], double in_fcumlambda[NUM_FCUM_LAMBDA+1][3], adType in_ftip[3],
                  adType out_x_np1[NUM_STATES], adType out_xdot_n[NUM_INTEGRATION_STATES] );


//
// Robotic Kinematics Related Functions
//

// cross product operator
// 1 dim array stored in row major order
template <typename adType>
void wHat(adType in_w[3], adType out_what[9]);

// Rodrigues' Formula for calculating expm(what*theta) for ||w||=1
//    basic formula
template <typename adType>
void RodriguesFormula (double in_w[3], double in_theta, double out_R[9]);
//    expanded form of the formula - presumably more efficient
template <typename adType>
void RodriguesExpanded (adType in_w[3], adType in_theta, adType out_R[9]);

// calculate R_np1 and p_np1 analytically using twist exponential, without numerical integration
//   g_np1 = g_n * expm ( \hat{\xi}^b *h ),  where \xi^b= [ 0 0 1 u_n^T ]^T
//   g = [R p; 0 0 0 1];
template <typename adType>
void SE3_Analytical_Step(adType in_R_n[9], adType in_p_n[3], adType in_u_n[3], double h, adType out_R_np1[9], adType out_p_np1[3]);


#define MAX(a,b) 	( ((b)>(a))?(b):(a) )
#define MIN(a,b) 	( ((b)<(a))?(b):(a) )
#define POW4(a)		( (a)*(a)*(a)*(a) )


