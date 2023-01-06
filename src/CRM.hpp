#pragma once
#include <cmath>
#include <tuple>
#include "CRM_MatrixOperations.hpp"
//  define compiler directive DO_NOT_USE_EIGEN  if you to use the Eigen-free version of CRM Forward Kinematics
//      #define DO_NOT_USE_EIGEN
#ifndef DO_NOT_USE_EIGEN
#include <eigen3/Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Matrix;
using Eigen::RowMajor;
#endif

#define ANALYTICAL_SE3_STEP

#define NUM_STATES 15   //
#define NUM_ACT_SET 1								// Number of actuator sets
#define NUM_FLEX_SEG 2								// Number of flexible segments
#define NUM_SEGMENTS (NUM_ACT_SET+NUM_FLEX_SEG)		// Total number of segments
#define NUM_LOCALIZATION_MARKERS 5					// Total number of localization markers
	// IMPORTANT NOTE: for now most proximal segment is assumed to be always flexible
	//    and the flexible and rigid segments are assumed to be alternating
	//    most distal segment can be flexible or rigid

#define NUM_FCUM_LAMBDA 104  	//  Number of steps used in calculating fcumlambda (cumulative forces); number of entries in fcumlambda is (NUM_FCUM+1)

// Regularization scales used for Nonlinear Solver
constexpr double IVALUE_SCALE_U = 1.0;		//(1.0)			// the variable used in Nonlinear Solver is multiplied with this scale to calculate u (curvature) that will be used in IVP
constexpr double IVALUE_SCALE_F = 1.0;		//(100.0)		// the variable used in Nonlinear Solver is multiplied with this scale to calculate ftip (tip force) that will be used in IVP
constexpr double RESIDUAL_SCALE_M = 1.0;	//(1.0)			// the residual for tip moment coming out of the IVP will be multiplied with this scale to return to the Nonlinear Solver
constexpr double RESIDUAL_SCALE_P = 1.0;	//(0.01)		// the residual for tip position error coming out of the IVP will be multiplied with this scale to return to the Nonlinear Solver

// NL Solver method selection
#define FK_TRUSTREGION								// Forward kinematics (inner loop - i.e., free space): Trust Region Method with the numerical jacobian
#define FK_TRUSTREGION_ANALYTICALJAC				//    Use the analytical jacobin in forward kinematics trustregion solution 
#define CONTACT_TRUSTREGION							// Forward kinematics (contact model): Trust Region Method with the numerical jacobian
#define CONTACT_TRUSTREGION_ANALYTICALJAC			//    Use the analytical jacobian in contact trustregion solution 
constexpr double TRUSTREGION_TOLERANCE = 1e-5;		// 

// Numerical Jacobian Calculation Stepsize - two orders of magnitude smaller than typical change scale
constexpr double NUM_JACOBIAN_CURRENT_STEPSIZE = 1e-5;
constexpr double NUM_JACOBIAN_INSERTIONLENGTH_STEPSIZE = 1e-2; 
constexpr double NUM_JACOBIAN_TIPFORCE_STEPSIZE = 1e-5;
constexpr double NUM_JACOBIAN_BASEPOSITION_STEPSIZE = 1e-3;
constexpr double NUM_JACOBIAN_BASEROTATION_STEPSIZE = 1e-3;
constexpr double NUM_JACOBIAN_BASECURVATURE_STEPSIZE = 1e-6;

// Uncomment the following if CRMShootingMethodBVP will be executed on FPGA PL as a kernel
//#define BVP_AS_A_KERNEL

// Enumerated type defining catheter contact mode.  
enum class ContactModeType { FREE_TIP, FIXED_TIP };
// FREE_TIP : no contact, FIXED_TIP : catheter tip is constrained at a given point
// Contact Mode: FREE_TIP ---  Dimension of NLEquation will be 3	// NLEquation - domain: u[0..2]; range: m_tip[0..2]  (moment at tip)
// Contact Mode: FIXED_TIP --- Dimension of NLEquation will be 6	// NLEquation domain: u[0..2], ftip[0..2]; range: m_tip[0..2], delta_p[0..2] (tip position error)


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


//
// ---------------------------------------------------------
//
// Cosserat Rod Model Forward Kinematics
//
// ---------------------------------------------------------
//

// Structure for defining catheter model parameters
struct CRMCatheterModelParams {
	//	For all parameters below, segments and actuator units are numbered/ordered from the tip of the catheter towards the base (distal to proximal)
	double 	SegLengths[NUM_SEGMENTS];				// Array of segment lengths; NUM_SEGMENTS long array
	double	LocMarkers[NUM_LOCALIZATION_MARKERS];	// Array of localization marker locations (in lambda coordinates, distances measured from the tip); NUM_LOCALIZATION_MARKERS long array
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
};

// Structure for defining catheter configuration parameters
struct CatheterConfiguration {
	double 	B0[3];		//  B0 field vector of the MRI scanner (in spatial coordinates)
	double 	g[3];		//  Gravity vector (in spatial coordinates)
	double 	p0[3];		//  Catheter entry port position (in spatial coordinates)
	double 	R0[9];		//	Catheter orientation at the entry port (relative to the spatial frame); 3x3 matrix stored in row major order
};

// Container for parameter data to be sent to Forward Kinematics
template <typename adType>
struct CRMForwardKinematicsData {
	CRMCatheterModelParams* CathParams;		// physical parameters of the catheter
	CatheterConfiguration*	CathConfig;		// catheter configuration in space
	ContactModeType ContactMode;			// contact mode
	double	TipConstraintPoint[3];			// spatial coordinates of the point that the catheter tip is constrained to - only used if ContactMode == ContactModeType::FIXED_TIP
	adType	TipForce[3];					// force applied at the tip of the catheter in spatial coordinates - only used if ContactMode == ContactModeType::FREE_TIP
	double	u0_initialguess[3];				// initial guess for the curvature at the base of the catheter
	double	ftip_initialguess[3];			// initial guess for the tip force in spatial coordinates - only used if ContactMode == ContactModeType::FIXED_TIP
	double	IntegrationStepSize;			// length stepsize used in numerical integration performed as part of IVP calculations (part of BVP)
	bool	FinalValueOnly;					// Flag used to indicate if only final value is returned (true) or if Marker Locations are returned as well (false)
	double	(*ReportedMarkerPos)[NUM_LOCALIZATION_MARKERS][3];  // output which returns the calculated marker locations
};

//useful macros
#define MAX(a,b) 	( ((b)>(a))?(b):(a) )
#define MIN(a,b) 	( ((b)<(a))?(b):(a) )
#define POW4(a)		( (a)*(a)*(a)*(a) )

// Include the internal definitions and the BVP/IVP declarations
#include "CRM_StateVector_Definitions.hpp"
#include "CRM_BVPIVP_APIDeclarations.hpp"
#include "CRM_IVP_Definitions.hpp"
#include "CRM_IVPJacobian_Declarations.hpp"
#include "CRM_IVPJacobian_Definitions.hpp"
#include "CRM_BVP_Definitions.hpp"

//
// Cosserat Rod Model Forward Kinematics API
//

// This is the C++ array version of the main API of the Cosserat Rod Model Forward Kinematics 
//		This is a wrapper for a sequence of CRMShootingMethodBVP + CRMSolverIVP calls
template <typename adType>
int CRM_ForwardKinematics(   // return value ==0 successful; !=0 if algorithms is stuck at a local minimum, or cannot make further progress
	adType in_x[],		// array of length (3 x NUM_ACT_SET + 1)  -- the +1 is for the inserted length of the catheter
						// Actuation currents for each of the actuation coil currents [NUM_ACT_SET][3] stored in row-major order 
						//    followed by the Inserted Length (length of the catheter from the entry point to the tip)
						//    actuator units are numbered/ordered from the tip of the catheter towards the base (distal to proximal)
	adType out_y[],		// output vector
						//    (3+9+3)x1 vector if Params.ContactMode == ContactModeType::FREE_TIP
						//    (3+9+3+3)x1 vector if Params.ContactMode == ContactModeType::FIXED_TIP		
						//    outputs are packed p[0..2],R[0..8],u0[0..2](,ftip[0..2])
						//		p:		tip position in spatial coordinates
						//      R:		tip orientation matrix --- 3x3 matrix stored in row major order (R11 R12 R13 R21 R22 R23 R31 R32 R33)
						//		u0:		curvature at the base of the catheter
						//		ftip:	tip force in spatial coordinates
	CRMForwardKinematicsData<adType> Params);

#ifndef DO_NOT_USE_EIGEN
// we will use c++ type traits to translate adType to corresponding Eigen vector type
template <typename T> struct adTypeVector_trait { typedef void type; };
template<> struct adTypeVector_trait<double> { typedef VectorXd type; };
//template<> struct adTypeVector_trait<real> { typedef autodiff::VectorXreal type; };
template<typename adType> using adTypeVector = typename adTypeVector_trait<adType>::type;

// This is the Eigen vectorized version of the main API of the Cosserat Rod Forward Kinematics 
//		This is a wrapper for a sequence of CRMShootingMethodBVP + CRMSolverIVP calls
template <typename adType>
adTypeVector<adType> CRM_ForwardKinematics(const adTypeVector<adType> in_x, CRMForwardKinematicsData<adType> Params, int& localmin);
//						in_x		input vector (vector version of in_x --- see above)
//						Params		Container for parameter data to be sent to Forward Kinematics
//						localmin	returns 0 if successful; !=0 if algorithms is stuck at a local minimum, or cannot make further progress
//						return		output vector (vector version of out_x --- see above)

//
// This is the API for the analytical calculation of the Forward Kinematics Jacobian 
//		( uses CRMSolverIVPJacobian )
//      (note that this function returns the hybrid manipulator Jacobian + extras, not the Jacobian of the FK map)
//
template <typename adType>
MatrixXd CRM_FKJacobian_Analytical(
	const adTypeVector<adType> in_x,				// Forward Kinematics input vector (see in_x of CRM_ForwardKinematics above)
	const adTypeVector<adType>& in_FKouty,			// Forward Kinematics output vector, corresponding to an equilibrium configuration of the catheter, passed as input (see out_y of CRM_ForwardKinematics above)
	CRMForwardKinematicsData<adType> in_Params);	// Container of the parameter data sent to Forward Kinematics
	// output matrix:
	//  if Params.ContactMode == ContactModeType::FREE_TIP
	//    output matrix is (3+9+3) x (3 x NUM_ACT_SET + 1 +3), packed dp/dz, dp/dft; ws_dz, ws_dft 
	//		dpdz	:	\partial p / \partial in_x			 :	partial derivatives of the tip position in spatial coordinates w.r.t. inputs (should be ~0 for ContactModeType::FIXED_TIP)
	//		dpdft	:	\partial p / \partial f_tip			 :	partial derivatives of the tip position in spatial coordinates w.r.t. tip force (should be ~0 for ContactModeType::FIXED_TIP)
	//      ws_dz	:	( \partial R / \partial in_x ) * R^T :	term corresponding to the spatial angular velocity of the tip
	//      ws_dft	:	( \partial R / \partial f_tip ) * R^T:	term corresponding to the spatial angular velocity of the tip
	//	  the first 3+3=6 rows of the output matrix gives the Hybrid Manipulator Jacobian
	//          see Murray, Li, Sastry Intro to Robotics Book for the definition of the Manipulator Jacobian, and the Hybrid Manipulator Jacobian
	//              and how it is different from the Jacobian of the FK map
	//  if Params.ContactMode == ContactModeType::FIXED_TIP		
	//    outputs matrix is (3) x (3 x NUM_ACT_SET + 1) matrix packed dftip/dz[0..2]
	//		dftip/dz:	\partial ftip / partial in_x		 :	partial derivative of the tip force in spatial coordinates w.r.t. inputs
	//

//
// This is the API for the numerical calculation of the Forward Kinematics Jacobian 
//      (note that this function returns the hybrid manipulator Jacobian + extras, not the Jacobian of the FK map)
//
template <typename adType>
MatrixXd CRM_FKJacobian_Numerical(
	const adTypeVector<adType> in_x,	// Forward Kinematics input vector (see in_x of CRM_ForwardKinematics above)
	adTypeVector<adType>&	out_y,		// Forward Kinematics output vector (see out_y of CRM_ForwardKinematics above)
	CRMForwardKinematicsData<adType> in_Params,	// Container for parameter data to be sent to Forward Kinematics
	int& localmin);						// returns 0 if successful; !=0 if algorithms is stuck at a local minimum, or cannot make further progress 
										//		(returned for forward kinematics calculation at in_x)
	// output matrix:
	//  if Params.ContactMode == ContactModeType::FREE_TIP
	//    output matrix is (3+9+3) x (3 x NUM_ACT_SET + 1), packed dp/di[0..2]; ws[0..2]; du0/di[0..2] 
	//		dpdz	:	\partial p / \partial in_x			 :	partial derivatives of the tip position in spatial coordinates w.r.t. inputs (should be ~0 for ContactModeType::FIXED_TIP)
	//      ws		:	( \partial R / \partial in_x ) * R^T :	term corresponding to the spatial angular velocity of the tip
	//		du0/dz	:	\partial u0 / \partial in_x			 :	partial derivative of curvature at the base of the catheter w.r.t. inputs
	//	  the first 3+3=6 rows of the output matrix gives the Hybrid Manipulator Jacobian
	//          see Murray, Li, Sastry Intro to Robotics Book for the definition of the Manipulator Jacobian, and the Hybrid Manipulator Jacobian
	//              and how it is different from the Jacobian of the FK map
	//  if Params.ContactMode == ContactModeType::FIXED_TIP		
	//    outputs matrix is (3) x (3 x NUM_ACT_SET + 1) matrix packed dftip/di[0..2]
	//		dftip/dz:	\partial ftip / partial in_x		 :	partial derivative of the tip force in spatial coordinates w.r.t. inputs
	//

#endif


#include "CRM_FK_InternalAPIDeclarations.hpp"
#include "CRM_FK_Definitions.hpp"

#include "CRMDYNIVPBVP_Declarations.hpp"
#include "CRMDYNIVPBVP_Defs.hpp"

#include "CoilDynamics_Defs.hpp"
