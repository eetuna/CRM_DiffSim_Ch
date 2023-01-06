#pragma once

//
// ---------------------------------------------------------
// This file contains the definition for the internally used 
//    components and interfaces for the 
//    Cosserot Rod Model of the Catheter
// ---------------------------------------------------------
//

//
// Implementation Note - 10/23/2021 MCC:
//    This version of the code has been implemented to be compatible with autodifferentiation (using autodiff library) using dual numbers.
//    Specifically, all of the CRM code has been templated such that all of the variables that would be differentiated (outputs) and variables that would be 
//    differentiated with respect to (inputs) are defined as the templated "adType".  Similarly, all of the local variables which has adType as an lvalue 
//    is also defined to be the adType.  All of the other i/o arguments of the API and the local variables remaing to be defined as regular double type.
//

// Implementation Note - 7/17/2022 MCC:
//    In the descriptions below, two different parametrizations are refered to when parametrizing locations along the length of the catheter, namely,
//    \lambda parameters, and s parameters:
//      \lambda refers to distances along the length of the cathteter measured from the tip towards the base (distal to proximal)
//      s refers to distances along the length of the cathteter measured from the entry point towards the top (proximal to distal)
//    It is important to note that s-parameters of specific locations/features changes as a function of the Insertion Length of the catheter whereas 
//      \lambda parameters are fixed.  
//   

// we will use c++ type traits to translate State Vector type to the correspnding State Derivative Vector type
template <typename T, typename adType> struct StDerivativeVect_trait { typedef void type; };
template<typename adType> struct StDerivativeVect_trait<StateVector<adType>,adType> { typedef StateDerivativeVector<adType> type; };
template<typename adType> struct StDerivativeVect_trait<AugmentedStateVector<adType, IVPJacobiansMini>, adType> { typedef AugmentedStateDerivativeVector<adType, IVPJacobiansMini> type; };
template<typename adType> struct StDerivativeVect_trait<AugmentedStateVector<adType,IVPJacobiansFull>,adType> { typedef AugmentedStateDerivativeVector<adType, IVPJacobiansFull> type; };
template<typename T, typename adType> using StDerivativeVectType = typename StDerivativeVect_trait<T,adType>::type;

// helper to be able to match state vector types
template <class T> using expr_type = std::remove_cv_t<std::remove_reference_t<T>>;

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
	adType	TipForce[3];
	double  CoilAlignmentTurnAreaMatrix[NUM_ACT_SET][9];
};

// Construct Shooting Method Parameter Set from Catheter Model and Configuration Parameters
template <typename adType>
void CRMConstructShootingMethodParamSet	(	CRMCatheterModelParams CathParams, CatheterConfiguration CathConfig,
											adType InsertionLength, adType ActuationCurrents[NUM_ACT_SET][3], 
											ContactModeType ContactMode,
											double TipConstraintPoint[3], adType TipForce[3],
											double IntegrationStepSize,
											CRMShootingMethodParams<adType> &ShootingParams);


//
// ---------------------------------------------------------
//        Cosserat Rod Model - Integrator for Solving the Initial Value Problem
// ---------------------------------------------------------
//

// CRMSolverIVP API which uses CRMShootinMethodParams data structure
//    note that in_ftip[] will be used, ignoring the value in_Params.TipForce[]
template <typename adType>
void CRMSolverIVP(	CRMShootingMethodParams<adType> in_Params,
					adType in_u0[3], adType in_ftip[3],
					bool in_FinalValueOnly,
					adType out_x_N[NUM_STATES], adType out_MomentResidual[3],
					double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]			
					);

// CRMSolverIVP API which exposes all of the individual parameters
template <typename adType>
void CRMSolverIVP (	double in_x_0[NUM_STATES], double in_IntegrationStepSize,
					double in_Li, double in_dlambdainv,
					double in_SegEndLambdas[NUM_SEGMENTS], double	in_LocMarkerLambdas[NUM_LOCALIZATION_MARKERS],
					double in_K[NUM_FLEX_SEG][9], double in_Kinv[NUM_FLEX_SEG][9], double in_ustar[NUM_FLEX_SEG][3],
					double in_CoilAlignmentTurnAreaMatrix[NUM_ACT_SET][9],
					adType in_MagMoment[NUM_ACT_SET][3], double in_fcumlambda[NUM_FCUM_LAMBDA+1][3], adType in_ftp[3],
					double in_B0[3],
					bool in_FinalValueOnly,
					double out_x_N[NUM_STATES], double out_MomentResidual[3],
					double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]
					);

// INPUT PARAMETERS
//	In all input parameters, segments, markers, and actuator units etc are numbered/ordered from the tip of the catheter towards the base (distal to proximal)
//double 	x_0[NUM_STATES];							// Initial value of the state at the entry point of the catheter
//														//  States are packed p[0..2],R[0..8],u[0..2]  (R: 3x3 matrix stored in row major order R11 R12 R13 R21 R22 R23 R31 R32 R33)
//double 	Li;											// Inserted Length (length of the catheter from the entry point to the tip)
//double 	dlambdainv;									// Reciprocal of \Delta \lambda (= \Delta s) used in discretizing fcum  // derived quantity ( dlambdainv = 1 / (Lf/NUM_FCUM) = NUM_FCUM/Lf ), Lf: functional length of the catheter
//double 	SegEndLambdas[NUM_SEGMENTS];				// Array of lambda values for segment endpoints; NUM_SEGMENTS long array
//double	LocMarkerLambdas[NUM_LOCALIZATION_MARKERS];	// Array of lambda values for localization markers; NUM_LOCALIZATION_MARKERS long array
//double 	K[NUM_FLEX_SEG][9];		  					// Catheter Rigidity Matrix; (NUM_FLEX_SEG)x9 long array, (NUM_FLEX_SEG) 3x3 matrices stored in row major order
														//		Actuator segments are assumed to be rigid
														//		In the future, we may want to just pass the parameters and construct the matrix inside
//double 	Kinv[NUM_FLEX_SEG][9];  					// Inverses of K matrices; (NUM_FLEX_SEG)x9 long array, (NUM_FLEX_SEG) 3x3 matrices stored in row major order  // derived quantity
// -UNUSED- double 	Kdot[NUM_FLEX_SEG][9]; 					// Derivative of K;  we assume Kdot=0.0
//double 	ustar[NUM_FLEX_SEG][3];						// Local curvature in unloaded configuration for each of the flexible segments; (NUM_FLEX_SEG)x3 long array, (NUM_FLEX_SEG) 3x1 vectors
														//		Actuator segments are assumed to be straight
// -UNUSED- double 	ustardot[NUM_FLEX_SEG][3];				// Derivative of local curvature in unloaded configuration; we assume ustardot=0.0 since our rest shape model is piecewise constant curvature
//double 	MagMoment[NUM_ACT_SET][3];					// Actuator magnetization moments; NUM_ACT_SETx3 long array, NUM_ACT_SET 3x1 vectors; MagMoment = CoilAlignMat * CoilTurnAreaMat * ActuationCurrentVector
//double 	fcumlambda[NUM_FCUM_LAMBDA+1][3];			// (NUM_FCUM_LAMBDA+1)x3 array (grouped by 3 doubles) storing cumulative external force (excluding tip force) integrated from \lambda = index * \Delta\lambda 
														//    to the catheter tip (\lambda=0) - in spatial (catheter base, i.e., entry port, frame) coordinates
														//    fcumlambda is parameterized by \lambda, distance from the tip, not s, distance from the entry point
														//    \lambda = Li - s
//double	ftip[3];									// External point force (in spatial coordinates) applied at the tip of the catheter (\lambda = 0)
//double 	B0[3];										//  B0 field vector of the MRI scanner (in spatial coordinates)
//bool  	FinalValueOnly;								// Flag used to indicate if only final value (xf) is returned (true) or if Marker Locations are returned as well (false)// OUTPUT VALUES
//	All output values are reported as numbered/ordered from the tip of the catheter towards the base (distal to proximal)
//double 	x_N[NUM_STATES];							// Final value of the state at the distal tip of the catheter
//														//    States are packed u[0..2],R[0..8],p[0..2]  (R: 3x3 matrix stored in row major order R11 R12 R13 R21 R22 R23 R31 R32 R33)
//double 	out_MomentResidual[3];						// Residual moment at the tip of the catheter, for use in boundary value problem solution
//double 	p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3];// Positions of markers (ordered distal to proximal) at lambda coordinates given in LocMarkerLambdas  --- filled if FinalValueOnly is false



//  Parameter Set used to call CRMSolverIVP_Core
template <typename adType>
struct CRMIVPCoreParams {
	//  These are the intermediate variables that will be used to call CRMSolverIVP_Core
	//	For all parameters below, segments and actuator units are numbered/ordered from the base of the catheter towards the tip (proximal to distal)
	adType xi[NUM_STATES];  							// Initial value of the state for the next segment to be integrated
														//  States are packed p[0..2],R[0..8],u[0..2]  (R: 3x3 matrix stored in row major order R11 R12 R13 R21 R22 R23 R31 R32 R33)
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
	double CoilAlignmentTurnAreaMatrix[NUM_ACT_SET][9]; // The product CoilAlignMat * CoilTurnAreaMat for each of the actuators; Na long array of 3x3 matrices (stored in row major order)

	double B0[3];										// B0 field vector of the MRI scanner (in spatial coordinates)
	adType MagMoment[NUM_ACT_SET][3]; 					// Actuator magnetization moments in body coordinates; Na*3 long array, Na 3x1 vectors; MagMoment = CoilAlignMat * CoilTurnAreaMat * ActuationCurrentVector

	int   StartSegmentIndex;							// Index of the segment where the integration to solve IVP will start -- the segment located at the entry point; note that segment indices start at 0
	int	  NextLocMarker;								// Next Localization Marker to be computed
	double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]; // positions of markers (ordered proximal to distal)
														//   only the entries 0..NextLocMarker-1 are filled
	};


// Preparation of CRMIVPCoreParams for subsequent call to CRMSolverIVP_Core
template <typename adType>
void CRMSolverIVP_Prep ( adType in_x_0[NUM_STATES], double in_IntegrationStepSize,
						 adType in_Li, double in_dlambdainv,
						 double in_SegEndLambdas[NUM_SEGMENTS], double in_LocMarkerLambdas[NUM_LOCALIZATION_MARKERS],
						 double in_K[NUM_FLEX_SEG][9], double in_Kinv[NUM_FLEX_SEG][9], double in_ustar[NUM_FLEX_SEG][3],
						 double in_CoilAlignmentTurnAreaMatrix[NUM_ACT_SET][9],
						 adType in_MagMoment[NUM_ACT_SET][3], double in_fcumlambda[NUM_FCUM_LAMBDA+1][3], 
						 double in_B0[3],
						 bool in_FinalValueOnly,
						 CRMIVPCoreParams<adType> &out_CoreParams	);


// Core Computations used in CRMSolverIVP - Integrator for Solving the Initial Value Problem
//    parameters are prepared using CRMSolverIVP_Prep
//    CRMSolverIVP_Core can be called multiple times by changing only u[0..2] components of xi (in_u[]) and in_ftip[]
//    for a single execution of CRMSolverIVP_Prep - the values in in_params should not be changed by user
template <typename adType>
void CRMSolverIVP_Core ( const CRMIVPCoreParams<adType>& in_params,
						 adType in_u[3], adType in_ftip[3],
						 StateVector<adType>& out_x_N, adType out_MomentResidual[3],
						 double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3] );


// Function for copying data in device memory to global memory --- used for dataflow pipelining
template <typename adType>
void CRMSolverIVP_Return (  const StateVector<adType>& in_x_N, adType in_MomentResidual[3],
							double in_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3],
							double out_x_N[NUM_STATES], double out_MomentResidual[3],
							double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]);


// support function to Propagate Boundary Condition through a Rigid Link - used by CRMSolverIVP_Core
template<typename adType>
void CRMSolverIVP_PropagateBCThroughRigidLink(StateVector<adType>& xi_ip1, adType Residual_ip1[3],
	const adType RigidSegmentLength, const unsigned int ActNo, const adType MagMoment[3], const double CoilAlignmentTurnAreaMatrix[9], 
	const double B0[3], const double ustar_i[3], const double K_i[9], const double ustar_ip1[3], const double Kinv_ip1[9],
	const StateVector<adType>& xf_i, const adType Residual_i[3]);
	//StateVector<adType> xi_ip1;	//  Initial value (xi) of the state for the next segment (i+1) to be integrated
	//adType Residual_ip1[3];		//  Residual at the end of the rigid link -- will be returned
	//adType RigidSegmentLength;	//  Length of the Rigid Segment
	//adType MagMoment[3]; 			//  Actuator magnetization moments in body coordinates for the actuator on the rigid link 
									//     MagMoment = CoilAlignMat * CoilTurnAreaMat * ActuationCurrentVector
	//double B0[3];					//  B0 field vector of the MRI scanner (in spatial coordinates)
	//double ustar_i[3];			//  Local curvature in unloaded configuration for the flexible segment before (immediately proximal to) the rigid link
	//double K_i[9]; 				//  Catheter Rigidity Matrix flexible segment before (immediately proximal to) the rigid link
	//double ustar_ip1[3];			//  Local curvature in unloaded configuration for flexible segment after (immediately distal to) the rigid link
	//double Kinv_ip1[9];			//  Inverse of K Matrix for the flexible segment after (immediately distal to) the rigid link
	//StateVector<adType> xf_i;		//  Final value (xf) of the state for the last segment (i) integrated
	//adType Residual_i[3];			//  Residual at the end of the last segment (i)


// support function for calculating location markers
template <typename adType>
void CalculateLocMarkers (  int& NextLocMarker, const StateVector<adType>& xnext, 
							const double LocMarkers[NUM_LOCALIZATION_MARKERS], adType SegBounds_ip1, 
							double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3] );
	//int& NextLocMarker; 									// next localization marker that needs to have an assigned value (updated value is returned)
	//double LocMarkers[NUM_LOCALIZATION_MARKERS];			// The s values of each of the localization markers (markers not yet inserted into the catheter would have a negative s value
	//adType SegBounds_ip1;									// The s value for the end of current segment = beginning of next (distal) segment boundary
	//StateVector<adType> xnext; 							// State vector for the endpoint of current segment
	//double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]; 	// output: positions of markers (ordered proximal to distal)

// support function to copy location marker positions
template <typename adType>
void LocMarkerUpdate(double p[3], const StateVector<adType>& xi, adType t);

template <typename adType>
void LocMarkerUpdate(double p[3], adType xi[NUM_STATES], adType t);


// support function - needed for autodiff compatibility
template <typename adType>
double dVal(adType x) { return (x); }

// Cosserat Rod Model Integrand
template <typename adType>
void CRMIntegrand (	adType s, const StateVector<adType>& in_x, const adType in_Li, const double in_dlambdainv,
					const double in_K[9], const double in_Kinv[9], const double in_l[3], const double in_ustar[3], const double in_fcumlambda[NUM_FCUM_LAMBDA+1][3],
					const adType in_nL[3],
					StateDerivativeVector<adType>& out_xdot);
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


//
// ---------------------------------------------------------
//        Cosserat Rod Model - Boundary Value Problem Solver
// ---------------------------------------------------------
//



// CRMShootingMethodBVP API which uses CRMShootinMethodParams data structure
//    note that when in_ContactMode == ContactModeType::FIXED_TIP inParams.TipForce[] will not be used 
//              when in_ContactMode == ContactModeType::FREE_TIP  in_ftip_initialguess[] will not be used, and out_ftip[] will be set to inParams.TipForce[]
template <typename adType>
void CRMShootingMethodBVP(CRMShootingMethodParams<adType> in_Params,
	double in_u0_initialguess[3], double in_ftip_initialguess[3],
	adType out_u0[3], adType out_ftip[3], int& out_localmin);

// CRMShootingMethodBVP API which exposes all of the individual parameters
template <typename adType>
void CRMShootingMethodBVP(
	ContactModeType in_ContactMode,					// Enumerated type defining catheter contact mode.  in_ContactMode == FREE_TIP if the catheter is not in contact with a surface, FIXED_TIP if catheter tip is constrained to be at TipContraintPoint
	double 	in_u0_initialguess[3],					// Initial guess for the local curvature vector at the entry point
	double 	in_ftip_initialguess[3],				// Initial guess for the tip force (\lambda=0), used when in_ContactMode == FIXED_TIP
	double	in_InsertedLength, 						// Inserted Length (length of the catheter from the entry point to the tip)
	adType 	in_ActuationCurrents[NUM_ACT_SET][3], 	// Actuation current for each of the actuation coils
													//    actuator units are numbered/ordered from the tip of the catheter towards the base (distal to proximal)
	double	in_TipConstraintPoint[3],				// The spatial coordinates of the point where the catheter tip is constrained to be (used if in_ContactMode == FIXED_TIP)
	adType	in_TipForce[3],							// External point force (in spatial coordinates) applied at the tip of the catheter (\lambda = 0) - this value will not be used if in_ContactMode == FIXED_TIP 
	double	in_IntegrationStepSize,					// Stepsize used in numerical integration along the length of the catheter
	// Catheter Configuration Parameters
	double 	in_B0[3],								// B0 field vector of the MRI scanner (in spatial coordinates)
	double 	in_g[3],								// Gravity vector (in spatial coordinates)
	double 	in_p0[3],								// Catheter entry port position (in spatial coordinates)
	double 	in_R0[9],								// Catheter orientation at the entry port (relative to the spatial frame); 3x3 matrix stored in row major order
	// Catheter Model Parameters
	//	For all parameters below, segments and actuator units are numbered/ordered from the tip of the catheter towards the base (distal to proximal)
	double 	in_SegLengths[NUM_SEGMENTS],			// Array of segment lengths; NUM_SEGMENTS long array
	double	in_LocMarkers[NUM_LOCALIZATION_MARKERS],// Array of localization marker locations (in lambda coordinates); NUM_LOCALIZATION_MARKERS long array
	double	in_InnerRadius[NUM_FLEX_SEG],			// Inner radii of the flexible catheter segments
	double	in_OuterRadius[NUM_FLEX_SEG],			// Outer radii of the flexible catheter segments
	double	in_YoungsModulus[NUM_FLEX_SEG],			// Youngs Moduli of the flexible catheter segments
	double	in_ShearModulus[NUM_FLEX_SEG],			// Shear Moduli of the flexible catheter segments
	double 	in_ustar[NUM_FLEX_SEG][3],				// Local curvature in unloaded configuration for each of the flexible segments; (NUM_FLEX_SEG)*3 long array, (NUM_FLEX_SEG) 3x1 vectors
													//		Actuator segments are assumed to be straight
	double 	in_CoilAlignmentAngles[NUM_ACT_SET][2],	// Coil Alignment Angles; NUM_ACT_SET*2 long array, for each actuator set, the angle for the first coil is relative to x axis, and the angle for the second coil is relative to y axis
	double 	in_CoilTurnAreaMat[NUM_ACT_SET][9], 	// Coil Turn Area matrices; NUM_ACT_SET*9 long array, NUM_ACT_SET 3x3 matrices stored in row major order
	double 	in_rho[NUM_SEGMENTS],					// Length density (mass per unit length) of the flexible catheter substrate (tubing); (NUM_SEGMENTS) long array
	double 	in_ActMass[NUM_ACT_SET],				// Actuator segment masses, does not include the flexible substrate; (NUM_ACT_SET) long array
	// Outputs
	adType 	out_u0[3], 								// Local curvature vector at the entry point calculated through the solution of BVP
	adType 	out_ftip[3],							// Tip force calculated through the solution of BVP (\lambda=0) --- used when in_ContactMode == FIXED_TIP 
	int& out_localmin								// out_localmin!=0 if algorithms is stuck at a local minimum, or cannot make further progress
);

// Preparation of NLEqnParams for subsequent call to Nonlinear Equation Solvers
void CRMShootingMethodBVP_Prep(
	double in_B0[3], double in_g[3], double in_p0[3], double in_R0[9],
	double in_SegLengths[NUM_SEGMENTS], double in_LocMarkers[NUM_LOCALIZATION_MARKERS],
	double in_InnerRadius[NUM_FLEX_SEG], double in_OuterRadius[NUM_FLEX_SEG],
	double in_YoungsModulus[NUM_FLEX_SEG], double in_ShearModulus[NUM_FLEX_SEG],
	double in_ustar[NUM_FLEX_SEG][3],
	double in_CoilAlignmentAngles[NUM_ACT_SET][2], double in_CoilTurnAreaMat[NUM_ACT_SET][9],
	double in_rho[NUM_SEGMENTS], double in_ActMass[NUM_ACT_SET],
	CRMCatheterModelParams& CathParams, CatheterConfiguration& CathConfig);

// wrapper for equation to be solved -- needed for CRMNonlinearSolver
template <typename adType>
struct NLEqnParams : CRMIVPCoreParams<adType> {
	ContactModeType ContactMode;			// Enumerated type defining catheter contact mode.  ContactMode == FREE_TIP if the catheter is not in contact with a surface, FIXED_TIP if catheter tip is constrained to be at TipContraintPoint
	double			TipConstraintPoint[3];	// The spatial coordinates of the point where the catheter tip is constrained to be (used if ContactMode == FIXED_TIP)
	adType			TipForce[3];			// External point force (in spatial coordinates) applied at the tip of the catheter (\lambda = 0) (used if ContactMode == FREE_TIP)
};

template <typename adType>
void CRM_NLEquation(adType in_x[], adType out_y[], NLEqnParams<adType> Params);

template<typename adType>
void CRM_NLEquation_AnalyticalJac(adType in_x[], adType out_y[], adType out_fjac[], NLEqnParams<adType> Params);

//
// Numerical Integration Support Functions
//

// In place projection of the State Vector to the appropriate manifold
//   this function is called by ABM4 after every iteration step
template<typename adType>
void Project_State_to_Manifold(StateVector<adType>& State);

//
// Numerical Integration Functions
//

//ABM4: 4th order Adams-Bashforth Prediction and Adams-Moulton Correction Numerical Integration
//    Note: The first three steps are calculated using RK2
template <typename adType, typename StVecType>
void ABM4 (	const StVecType& in_x_0, const adType t_0, const int N, const adType h,
			const adType Li, const double dlambdainv, const double in_K[9], const double in_Kinv[9], const double in_l[3], 
			const double in_ustar[3], const double in_fcumlambda[NUM_FCUM_LAMBDA+1][3], const adType in_nL[3],
			const bool FinalValueOnly, const double in_LocMarkers[NUM_LOCALIZATION_MARKERS], int &inout_NextLocMarkerIdx,
			StVecType& out_x_N, double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]	);

//ABM4_step One step of 4th order Adams-Bashforth Prediction and Adams-Moulton Correction
template <typename adType, typename StVecType>
void ABM4_step(	const StVecType& in_x_n, adType t_n, adType h,
				const StDerivativeVectType<StVecType, adType>& in_xdot_nm1, const StDerivativeVectType<StVecType, adType>& in_xdot_nm2, const StDerivativeVectType<StVecType, adType>& in_xdot_nm3,
				const StVecType& in_x_nm1, const StVecType& in_x_nm2, const StVecType& in_x_nm3,
				const adType Li, const double dlambdainv, const double in_K[9], const double in_Kinv[9], const double in_l[3], 
				const double in_ustar[3], const double in_fcumlambda[NUM_FCUM_LAMBDA + 1][3], const adType in_nL[3],
				StVecType& out_x_np1, StDerivativeVectType<StVecType, adType>& out_xdot_n);

//RK2_step One step of 2nd Order Runge-Kutta Integration
template <typename adType, typename StVecType>
void RK2_step(	const StVecType& in_x_n, const adType t_n, const adType h,
				const adType Li, const double dlambdainv, const double in_K[9], const double in_Kinv[9], const double in_l[3], 
				const double in_ustar[3], const double in_fcumlambda[NUM_FCUM_LAMBDA+1][3], const adType in_nL[3],
				StVecType& out_x_np1, StDerivativeVectType<StVecType, adType>& out_xdot_n );


//
// Robotic Kinematics Related Functions
//

// cross product operator
// 1 dim array stored in row major order
template <typename adType>
void wHat(const adType in_w[3], adType out_what[9]);
// 1 dim array stored in row major order, stride: pointer increment between two consecutive entries 
template <typename adType>
void wHat(const adType in_w[3], adType out_what[9], unsigned int stride);

// cross product operator
// extract the 3dim vector from the so(3) skew symmetric matrix stored in row major order
template <typename adType>
void vee_from_so3(const adType in_what[9], adType out_w[3]);

// Rodrigues' Formula for calculating expm(what*theta) for ||w||=1
//    basic formula
template <typename adType>
void RodriguesFormula (adType in_w[3], adType in_theta, adType out_R[9]);
//    expanded form of the formula - presumably more efficient
template <typename adType>
void RodriguesExpanded (adType in_w[3], adType in_theta, adType out_R[9]);

// calculate R_np1 and p_np1 analytically using twist exponential, without numerical integration
//   g_np1 = g_n * expm ( \hat{\xi}^b *h ),  where \xi^b= [ 0 0 1 u_n^T ]^T
//   g = [R p; 0 0 0 1];
template <typename adType>
void SE3_Analytical_Step(adType in_R_n[9], adType in_p_n[3], adType in_u_n[3], adType h, adType out_R_np1[9], adType out_p_np1[3]);

