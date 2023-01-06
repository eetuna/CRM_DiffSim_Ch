#pragma once

//
// ---------------------------------------------------------
// This file contains the definition for the 
//    analytic Jacobian of the IVP solver of the
//    Cosserot Rod Model of the Catheter
// ---------------------------------------------------------
//

//
// Implementation Note - 8/11/2022 MCC:
//    The CRMJacobianIntegrand functions assumes that:
//	    - the distributed moment load density l is assumed to be zero
//      - the K and Kinv matrices are assumed to be diagonal
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


//
// ---------------------------------------------------------
//        Cosserat Rod Model - Integrator for Calculatinf the Jacobian of the Initial Value Problem Solver
// ---------------------------------------------------------
//

// CRMSolverIVPJacobian API which uses CRMShootinMethodParams data structure
//    note that in_ftip[] will be used, ignoring the value in_Params.TipForce[]
template <typename adType>
std::tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd> CRMSolverIVPJacobian(  
		CRMShootingMethodParams<adType> in_Params,
		adType in_u0[3], adType in_ftip[3],
		bool in_FinalValueOnly,
		adType out_x_N[NUM_STATES], adType out_MomentResidual[3]);


// Core Computations used for CRMSolverIVP Jacobian Calculation
//    parameters are prepared using CRMSolverIVP_Prep
//    the values in in_u and in_ftp should be those coming from an equilibrium calculation
//    if needed. CRMSolverIVP_CorewithJacobian can be called multiple times by changing only u[0..2] components of xi (in_u[]) and in_ftip[]
//    for a single execution of CRMSolverIVP_Prep - the values in in_params should not be changed by user
template <typename adType, template<typename> typename IVPJacobians>
void CRMSolverIVP_CoreWithJacobian( CRMIVPCoreParams<adType> in_params,
									adType in_u[3], adType in_ftip[3],
									AugmentedStateVector<adType, IVPJacobians>& out_x_N, adType out_MomentResidual[3]);

template<typename adType, template<typename> typename IVPJacobians>
void CRMSolverIVP_PropagateBCThroughRigidLink( AugmentedStateVector<adType, IVPJacobians>& xi_ip1, adType Residual_ip1[3],
	const adType RigidSegmentLength, const unsigned int ActNo, const adType MagMoment[3], const double CoilAlignmentTurnAreaMatrix[9],
	const double B0[3], const double ustar_i[3], const double K_i[9], const double ustar_ip1[3], const double Kinv_ip1[9],
	const AugmentedStateVector<adType, IVPJacobians>& xf_i, const adType Residual_i[3]);


// Cosserat Rod Model IVP Integrand for Jacobian Calculation
template <typename adType, template<typename> typename IVPJacobians>
void CRMIntegrand (	const adType s, const AugmentedStateVector<adType, IVPJacobians>& in_x, const adType in_Li, const double in_dlambdainv,
					const double in_K[9], const double in_Kinv[9], const double in_l[3], const double in_ustar[3], const double in_fcumlambda[NUM_FCUM_LAMBDA+1][3],
					const adType in_ftip[3],
					AugmentedStateDerivativeVector<adType, IVPJacobians>& out_xdot);


// helper functions
// each column k of  Pnew_k = P_k + (w_k^ R_k)_3 * d, where (w_k^ * R_k)_3 is the third column of the (w_k^ * R_k) matrix ( 3x3 matrices packed in row major order )
template<unsigned int dim, typename adType>
inline void PpluswhatR3timesD(const adType P[], const adType w[], const adType R[], const adType d, adType Pnew[]);

