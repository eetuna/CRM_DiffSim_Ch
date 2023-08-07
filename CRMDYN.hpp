#pragma once
#include <cmath>
#include "CRMMatrixOperations.hpp"
#include <math.h>
#include <stdlib.h>

#include <chrono>

using namespace std::chrono;

#define NUM_STATES 15   // u[0..2],R[0..9],p[0..2]
#define ANALYTICAL_SE3_STEP
#ifdef ANALYTICAL_SE3_STEP
	#define NUM_INTEGRATION_STATES 3   // u[0..2]
#else
	#define NUM_INTEGRATION_STATES 15   // u[0..2],R[0..9],p[0..2]
#endif
#define EQNDIMENSION 6	// domain: u[0..2], range: m_tip[0..2]

#define NUM_ACT_SET 1								// Number of actuator sets
#define NUM_FLEX_SEG 2								// Number of flexible segments
#define NUM_SEGMENTS (NUM_ACT_SET+NUM_FLEX_SEG)		// Total number of segments
#define NUM_LOCALIZATION_MARKERS 5					// Total number of localization markers
	// IMPORTANT NOTE: for now most proximal segment is assumed to be always flexible
	//    and the flexible and rigid segments are assumed to be alternating
	//    most distal segment can be flexible or rigid

#define NUM_COIL_STATES (6+3+9) // v[3], w[3], p[3] ,R[9]


#define NUM_RESIDUAL 3
#define NUM_DYN_RESIDUAL (NUM_ACT_SET*6) // (m+n) * NUM_ACT_SET
//#define DELTA_T 0.01
#define NUM_CONTROL (NUM_ACT_SET * 3 + 1) // Control dimension for the dynamic model
#define RESIDUAL_SCALE_F	1.0	// the residual for coil force coming out of the IVP will be multiplied with this scale to return to the Nonlinear Solver
#define RESIDUAL_SCALE_M	1.0 //(10.0)			// the residual for tip moment coming out of the IVP will be multiplied with this scale to return to the Nonlinear Solver


#define NUM_FCUM_LAMBDA 104  	//  Number of steps used in calculating fcumlambda (cumulative forces); number of entries in fcumlambda is (NUM_FCUM+1)

// Regularization scales used for Nonlinear Solver

#define IVALUE_SCALE_M  1.0 //
#define IVALUE_SCALE_N	1.0		// the variable used in Nonlinear Solver is multiplied with this scale to calculate ftip (tip force) that will be used in IVP
#define IVALUE_SCALE_U	1.0 //(0.01)			// the variable used in Nonlinear Solver is multiplied with this scale to calculate u (curvature) that will be used in IVP
#define IVALUE_SCALE_F	0.01//(0.01)			// the variable used in Nonlinear Solver is multiplied with this scale to calculate ftip (tip force) that will be used in IVP
#define RESIDUAL_SCALE_P	1e+1 //(10.0)			// the residual for tip position error coming out of the IVP will be multiplied with this scale to return to the Nonlinear Solver
#define RESIDUAL_SCALE_R	1e+1 //(10.0)			// the residual for tip position error coming out of the IVP will be multiplied with this scale to return to the Nonlinear Solver

// NL Solver method selection
#define TRUSTREGION							// Trust Region Method with the numerical jacobian (Default)

#define MAX(a,b) 	( ((b)>(a))?(b):(a) )
#define MIN(a,b) 	( ((b)<(a))?(b):(a) )
#define POW4(a)		( (a)*(a)*(a)*(a) )


// Uncomment the following if CRMShootingMethodBVP will be executed on FPGA PL as a kernel
//#define BVP_AS_A_KERNEL

// Enumerated type defining catheter contact mode.  
enum class ContactModeType { FREE_TIP, FIXED_TIP };
// FREE_TIP : no contact, FIXED_TIP : catheter tip is constrained at a given point
// Contact Mode: FREE_TIP ---  EQNDIMENSION will be 3	// NLEquation - domain: u[0..2]; range: m_tip[0..2]  (moment at tip)
// Contact Mode: FIXED_TIP --- EQNDIMENSION will be 6	// NLEquation domain: u[0..2], ftip[0..2]; range: m_tip[0..2], delta_p[0..2] (tip position error)


//used to add the extra variables for returning debug data from the kernel
//#define ADD_DEBUG
#ifdef ADD_DEBUG
#define DEBUG_BUFFER_LIMIT 2000 // size of the debug message buffer
#endif

//
// Implementation Note - 10/23/2021 MCC:
//    This version of the code has been implemented to be compatible with autodifferentiation (using autodiff library) using dual numbers.
//    Specifically, all of the CRM code has been templated such that all of the variables that would be differentiated (outputs) and variables that would be 
//    differentiated with respect to (inputs) are defined as the templated "adType".  Similarly, all of the local variables which has adType as an lvalue 
//    is also defined to be the adType.  All of the other i/o arguments of the API and the local variables remaing to be defined as regular double type.
//


////
//// ---------------------------------------------------------
////
//// Cosserat Rod Model Kinematics
////
//// ---------------------------------------------------------
////

// Structure for defining catheter model parameters
struct CRMCatheterModelParams {
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
    double  ActInertia[NUM_ACT_SET][9];             // Inertia matrix of the coils
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
struct CRMShootingMethodParams {
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
	ContactModeType ContactMode;
	double	TipConstraintPoint[3];
	double	TipForce[3];

    double 	g[3];
    double v_L_pre[NUM_ACT_SET][3];                                  // The linear velocity at the coil (L)
    double w_L_pre[NUM_ACT_SET][3];                                  // The angular velocity at the coil (L)
    double p_pre[NUM_ACT_SET][3];
    double R_pre[NUM_ACT_SET][9];
    double actMass[NUM_ACT_SET];
    double actInertia[NUM_ACT_SET][9];
    double damping[NUM_ACT_SET][6];
    double DELTA_T;

};


////
//// ---------------------------------------------------------
////        Cosserat Rod Model - Integrator for Solving the Initial Value Problem
//// ---------------------------------------------------------
////

// CRMSolverIVP API which uses CRMShootinMethodParams data structure
//    note that in_ftip[] will be used, ignoring the value in_Params.TipForce[]
template <typename adType>
void CRMSolverIVP(	CRMShootingMethodParams<adType> in_Params,
                      adType in_u0[3], adType in_ftip[3],
                      bool in_FinalValueOnly,
                      adType out_x_N[NUM_STATES], adType out_MomentResidual[NUM_RESIDUAL], adType pcoil[NUM_ACT_SET][3], adType Rcoil[NUM_ACT_SET][9],
                      double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]);


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
														//   only the entries 0..NextLocMarker-1 are filled

    double 	g[3];		                                //  Gravity vector (in spatial coordinates)
    double v_L_pre[NUM_ACT_SET][3];                                  // The linear velocity at the coil (L) at the previous time period
    double w_L_pre[NUM_ACT_SET][3];                                  // The angular velocity at the coil (L) at the previous time period
    double actMass[NUM_ACT_SET];
    double actInertia[NUM_ACT_SET][9];

    double p_pre[NUM_ACT_SET][3];
    double R_pre[NUM_ACT_SET][9];
    double m_L[NUM_ACT_SET][3];
    double n_L[NUM_ACT_SET][3];
    bool DYN_Core; // true is we are calculating the coer function inside the dynamics function
    double damping[NUM_ACT_SET][6]; //v and w of the coil
    double DELTA_T;
	};


template <typename adType>
void CRMSolverIVP_Prep ( adType in_x_0[NUM_STATES], double in_IntegrationStepSize,
                         adType in_Li, double in_dlambdainv,
                         double in_SegEndLambdas[NUM_SEGMENTS], double in_LocMarkerLambdas[NUM_LOCALIZATION_MARKERS],
                         double in_K[NUM_FLEX_SEG][9], double in_Kinv[NUM_FLEX_SEG][9], double in_ustar[NUM_FLEX_SEG][3],
                         adType in_MagMoment[NUM_ACT_SET][3], double in_fcumlambda[NUM_FCUM_LAMBDA+1][3],
                         double in_B0[3], double in_g[3], const double in_actMass[NUM_ACT_SET], double in_actInertia[NUM_ACT_SET][9],
                         double in_damping[NUM_ACT_SET][6], double in_delta_t,
                         double in_v_L_pre[NUM_ACT_SET][3], double in_w_L_pre[NUM_ACT_SET][3], double in_p_pre[NUM_ACT_SET][3], double in_R_pre[NUM_ACT_SET][9],
                         double in_mL[NUM_ACT_SET][3], double in_nL[NUM_ACT_SET][3],
                         bool in_DYNCORE, bool in_FinalValueOnly, CRMIVPCoreParams<adType> &out_CoreParams) ;

/**
 *  Core Computations used in CRMSolverIVP - Integrator for Solving the Initial Value Problem
 *  parameters are prepared using CRMSolverIVP_Prep
 * @tparam adType
 * @param in_params
 * @param in_u
 * @param in_nL
 * @param in_tauL : only use this in the dynamic core, in quasi-static case, it should be 0
 * @param in_ftip
 * @param out_x_N
 * @param out_Residual
 * @param out_Tbcoil
 * @param out_pcoil
 * @param out_Rcoil
 * @param out_p_atLocMarkers
 */
template <typename adType>
void CRMSolverIVP_Core ( CRMIVPCoreParams<adType> in_params,
                         adType in_u[3], adType in_ftip[3],
                         adType out_x_N[NUM_STATES], adType out_Residual[NUM_RESIDUAL], adType out_pcoil[NUM_ACT_SET][3], adType out_Rcoil[NUM_ACT_SET][9],
                         double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]);


// support function to copy location marker positions
template <typename adType>
void LocMarkerUpdate(double p[3], adType xi[NUM_STATES], adType t);

// support function - needed for autodiff compatibility
template <typename adType>
double dVal(adType x) { return (x); }

// Cosserat Rod Model Integrand
template <typename adType>
void CRMIntegrand (	adType s, adType x[NUM_STATES], adType Li, double dlambdainv,
                       const double in_K[9], const double in_Kinv[9], const double in_l[3], const double in_ustar[3], adType in_nL[3],
                       double in_fcumlambda[NUM_FCUM_LAMBDA+1][3], adType in_ftip[3],
                       adType xdot[NUM_INTEGRATION_STATES]) ;
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


////
//// ---------------------------------------------------------
////        Cosserat Rod Model - Boundary Value Problem Solver
//// ---------------------------------------------------------
////

// CRMShootingMethodBVP API which uses CRMShootinMethodParams data structure
//    note that when in_ContactMode == ContactModeType::FIXED_TIP inParams.TipForce[] will not be used 
//              when in_ContactMode == ContactModeType::FREE_TIP  in_ftip_initialguess[] will not be used, and out_ftip[] will be set to inParams.TipForce[]
template <typename adType>
void CRMShootingMethodBVP(	CRMShootingMethodParams<adType> in_Params,
                              const double in_u0_initialguess[3], const double in_ftip_initialguess[3],
                              adType out_u0[3], adType out_ftip[3], int& out_localmin);
//

// Preparation of NLEqnParams for subsequent call to Nonlinear Equation Solvers
void CRMShootingMethodBVP_Prep (
        double in_B0[3], double in_g[3], double in_p0[3], double in_R0[9],
        double in_SegLengths[NUM_SEGMENTS], double in_LocMarkers[NUM_LOCALIZATION_MARKERS],
        double in_InnerRadius[NUM_FLEX_SEG],	double in_OuterRadius[NUM_FLEX_SEG],
        double in_YoungsModulus[NUM_FLEX_SEG], double in_ShearModulus[NUM_FLEX_SEG],
        double in_ustar[NUM_FLEX_SEG][3],
        double in_CoilAlignmentAngles[NUM_ACT_SET][2], double in_CoilTurnAreaMat[NUM_ACT_SET][9],
        double in_rho[NUM_SEGMENTS],	double in_ActMass[NUM_ACT_SET], double in_ActInertia[NUM_ACT_SET][9],
        CRMCatheterModelParams &CathParams, CatheterConfiguration &CathConfig );

// wrapper for equation to be solved -- needed for CRMNonlinearSolver
template <typename adType>
struct NLEqnParams : CRMIVPCoreParams<adType> {
	ContactModeType ContactMode;			// Enumerated type defining catheter contact mode.  ContactMode == FREE_TIP if the catheter is not in contact with a surface, FIXED_TIP if catheter tip is constrained to be at TipContraintPoint
	double			TipConstraintPoint[3];	// The spatial coordinates of the point where the catheter tip is constrained to be (used if ContactMode == FIXED_TIP)
	double			TipForce[3];			// External point force (in spatial coordinates) applied at the tip of the catheter (\lambda = 0) (used if ContactMode == FREE_TIP)

    // terms used in the dynamics
//    double	u0_initialguess[3];
    double	ftip_initialguess[3];
    double xf[NUM_STATES]; //Tip state, if applicable
};


template <typename adType>
void NLEquation(adType in_x[], adType out_y[], NLEqnParams<adType> Params);

// Construct Shooting Method Parameter Set from Catheter Model and Configuration Parameters
template <typename adType>
void CRMConstructShootingMethodParamSet(	CRMCatheterModelParams CathParams, CatheterConfiguration CathConfig,
                                            adType InsertionLength, adType ActuationCurrents[NUM_ACT_SET][3],
                                            ContactModeType ContactMode,
                                            double TipConstraintPoint[3], double TipForce[3],
                                            double IntegrationStepSize,
                                            double in_v_L_pre[NUM_ACT_SET][3], double in_w_L_pre[NUM_ACT_SET][3], double in_p_pre[NUM_ACT_SET][3], double in_R_pre[NUM_ACT_SET][9], double in_damping[NUM_ACT_SET][6], double in_DELTA_T,
                                            CRMShootingMethodParams<adType> &ShootingParams);
////
//// ---------------------------------------------------------
////        Cosserat Rod Model - Dynamics Functions
//// ---------------------------------------------------------
////

template <typename adType>
struct CRMDynamicsParams {
    CRMCatheterModelParams* CathParams;
    ContactModeType ContactMode;
    double	TipConstraintPoint[3];
    double	TipForce[3];
    double	u0_initialguess[3];
    double	nL_initialguess[3];
    double	mL_initialguess[3];
    double	ftip_initialguess[3];
    double	IntegrationStepSize;
    bool	FinalValueOnly;				// Flag used to indicate if only final value is returned (true) or if Marker Locations are returned as well (false)
    double	(*ReportedMarkerPos)[NUM_LOCALIZATION_MARKERS][3];  // output
};


template <typename adType>
void CRMShootingMethodBVP_DYN(	NLEqnParams<adType> in_Params, adType out_u0[3], adType out_ftip[3], int& out_localmin);

template <typename adType>
void CoilDynamics( adType in_coil_state[NUM_COIL_STATES], adType in_n[3], adType g[3],
                   adType actMass, adType actInertia[9], adType damping[6], double DELTA_T, adType in_B0[3], adType in_muhat[9], adType in_mL[3], adType out_coil_state[NUM_COIL_STATES], adType out_xdot_n[6]);

template <typename adType>
void CoilIntegrad(adType in_twist[6], adType in_n[3], adType g[3], adType R[9], adType actMass, adType actInertia[9], adType damping[6],adType in_B0[3], adType in_muhat[9], adType in_mL[3], adType twistdot[6]);

template <typename adType>
void SE3_TimeSpace(adType in_R_n[9], adType in_p_n[3], adType h, adType in_twist_n[6], adType out_R_np1[9], adType out_p_np1[3]) ;

template <typename adType>
void RK2_coildyn(adType in_x_n[NUM_COIL_STATES], adType in_n[3], adType g[3],  adType actMass, adType actInertia[9], adType damping[6], adType in_B0[3], adType in_muhat[9], adType in_mL[3],
                 adType out_x_np1[NUM_COIL_STATES], adType out_xdot_n[6] ) ;

template <typename adType>
void ABM4_coildyn(	adType in_x_n[NUM_COIL_STATES],adType in_xdot_nm1[6], adType in_xdot_nm2[6], adType in_xdot_nm3[6],
                      adType in_x_nm1[NUM_COIL_STATES], adType in_x_nm2[NUM_COIL_STATES], adType in_x_nm3[NUM_COIL_STATES],
                      adType in_n[3], adType g[3],  adType actMass, adType actInertia[9], adType damping[6], adType in_B0[3], adType in_muhat[9], adType in_mL[3],
                      adType out_x_np1[NUM_COIL_STATES], adType out_xdot_n[6]);

template <typename adType>
void DYNNLEquation(adType in_x[], adType out_y[], NLEqnParams<adType> Params, adType out_u0[3]);

template <typename adType>
void DynamicsBVP(CRMShootingMethodParams<adType> in_Params, const double xf[NUM_STATES],
                 double in_mL_initialguess[NUM_ACT_SET][3], double in_nL_initialguess[NUM_ACT_SET][3], double in_ftip_initialguess[3],
                 adType out_u0[3], adType out_mL[NUM_ACT_SET][3], adType out_nL[NUM_ACT_SET][3], adType out_ftip[3], int& out_localmin) ;

template <typename adType>
void DYNSolverIVP(	CRMShootingMethodParams<adType> in_Params, adType in_u0[3],
                      adType in_mL[NUM_ACT_SET][3], adType in_nL[NUM_ACT_SET][3], adType in_ftip[3],
                      bool in_FinalValueOnly,
                      adType out_x_N[NUM_STATES], adType out_coil_state[NUM_ACT_SET][NUM_COIL_STATES],
                      double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]);

template <typename adType>
void CRMFlexible_IVP_Back ( int SegmentIndex, const double in_p[3], const double in_R[9],  CRMIVPCoreParams<adType> in_params,
                            const double in_u[3], const double in_n_L[3],
                            adType out_u[3], double out_p[3], double out_R[9]);

/**
 *  The update function: runs inside the Dynamic functions to pass dynamic data (parameters) that are updated during the loop to the shooting method
 * @tparam adType
 * @param CathParams
 * @param in_p0
 * @param in_R0
 * @param in_v0
 * @param in_w0
 * @param InsertionLength
 * @param ActuationCurrents
 * @param ContactMode
 * @param TipConstraintPoint
 * @param TipForce
 * @param IntegrationStepSize
 * @param in_v_L_pre
 * @param in_w_L_pre
 * @param ShootingParams
 */
template <typename adType>
void CRMDYNCoreDataUpdate(CRMCatheterModelParams CathParams, double in_p0[3], double in_R0[9],
                          adType InsertionLength, adType ActuationCurrents[NUM_ACT_SET][3],
                          ContactModeType ContactMode, double TipConstraintPoint[3], double TipForce[3], double IntegrationStepSize,
                          double in_v_L_pre[3], double in_w_L_pre[3], double pL_pre[3], double RL_pre[9],
                          CRMShootingMethodParams<adType> &ShootingParams);


////
//// Numerical Integration Functions
////

//ABM4: 4th order Adams-Bashforth Prediction and Adams-Moulton Correction Numerical Integration
//    Note: The first three steps are calculated using RK2
template <typename adType>
void ABM4 (adType in_x_0[NUM_STATES], adType t_0, int N, double h,
           adType Li, double dlambdainv, double in_K[9], double in_Kinv[9], double in_l[3], double in_ustar[3], adType in_nL[3],
           double in_fcumlambda[NUM_FCUM_LAMBDA+1][3], adType in_ftip[3], bool FinalValueOnly,
           double	in_LocMarkers[NUM_LOCALIZATION_MARKERS], int *inout_NextLocMarkerIdx,
           adType out_x_N[NUM_STATES], double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3] );

//ABM4_step One step of 4th order Adams-Bashforth Prediction and Adams-Moulton Correction
template <typename adType>
void ABM4_step(	adType in_x_n[NUM_STATES], adType t_n, double h,
				adType in_xdot_nm1[NUM_INTEGRATION_STATES], adType in_xdot_nm2[NUM_INTEGRATION_STATES], adType in_xdot_nm3[NUM_INTEGRATION_STATES],
				adType in_x_nm1[NUM_STATES], adType in_x_nm2[NUM_STATES], adType in_x_nm3[NUM_STATES],
				adType Li, double dlambdainv, double in_K[9], double in_Kinv[9], double in_l[3], double in_ustar[3],  double in_nL[3],  double in_fcumlambda[NUM_FCUM_LAMBDA + 1][3], adType in_ftip[3],
                adType out_x_np1[NUM_STATES], adType out_xdot_n[NUM_INTEGRATION_STATES]);

//RK2_step One step of 2nd Order Runge-Kutta Integration
template <typename adType>
void RK2_step(	adType in_x_n[NUM_STATES], adType t_n, double h,
				adType Li, double dlambdainv, double in_K[9], double in_Kinv[9], double in_l[3], double in_ustar[3],  double in_nL[3],
                  double in_fcumlambda[NUM_FCUM_LAMBDA+1][3], adType in_ftip[3],
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

template <typename adType>
void rotationMatrixToEulerAngles(adType in_R[9], adType out_v[3]);


#include "CRMIVP_Defs.hpp"
#include "CRMIVP_Defs _backwards.hpp"
#include "CRMBVP_Defs.hpp"
#include "CRMDynamics.hpp"
#include "CoilDynamics_Defs.hpp"
