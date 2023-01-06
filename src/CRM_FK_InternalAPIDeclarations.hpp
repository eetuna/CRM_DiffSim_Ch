#pragma once

// This is the free space model - called by CRM_ForwardKinematics()
//    this function should not be called independently
template <typename adType>
int CRM_ForwardKinematics_FreeSpace(adType in_x[], adType out_y[], CRMForwardKinematicsData<adType> Params);

// This is the contact model - called by CRM_ForwardKinematics()
//    this function should not be called independently
template <typename adType>
int CRM_ForwardKinematics_Contact(adType in_x[], adType out_y[], CRMForwardKinematicsData<adType> Params);

#ifndef DO_NOT_USE_EIGEN
// This is the free space model - called by CRM_ForwardKinematics()
//    this function should not be called independently
template <typename adType>
adTypeVector<adType> CRM_ForwardKinematics_FreeSpace(const adTypeVector<adType> in_x, CRMForwardKinematicsData<adType> Params, int& localmin);

// This is the contact model - called by CRM_ForwardKinematics()
//    this function should not be called independently
template <typename adType>
adTypeVector<adType> CRM_ForwardKinematics_Contact(const adTypeVector<adType> in_x, CRMForwardKinematicsData<adType> Params, int& localmin);


//
// This is the API for the Forward Kinematics Jacobian for the Free Space Model - called by CRM_FKJacobian_Numerical()
//    this function should not be called independently
//       (note that this function returns the hybrid manipulator Jacobian + extras, not the Jacobian of the FK map)
//
// Enumerated type defining catheter contact mode.  
enum class FKFreeJacobianType { ACTUATION_ONLY, ACTUATION_AND_TIP_FORCE };
// ACTUATION_ONLY: Jacobian is calculated only w.r.t. actuation inputs
// ACTUATION_AND_TIP_FORCE: Jacobian is calculated only w.r.t. actuation inputs and tip force
template <typename adType>
MatrixXd CRM_FKJacobian_FreeSpace(
	const adTypeVector<adType> in_x,	// Forward Kinematics input vector (see in_x of CRM_ForwardKinematics)
	adTypeVector<adType>& out_y,		// Forward Kinematics output vector (see out_y of CRM_ForwardKinematics)
	CRMForwardKinematicsData<adType> in_Params,	// Container for parameter data to be sent to Forward Kinematics
	FKFreeJacobianType mode,			// flag to indicate if the Jacobian needs to be calculate w.r.t, actuation inputs only, or alos the tip force
	int& localmin);						// returns 0 if successful; !=0 if algorithms is stuck at a local minimum, or cannot make further progress 
										//		(returned for forward kinematics calculation at in_x)
	// output matrix rows are packed dp[0..2]; ws[0..2]; du0[0..2]
	// output matrix columns are packed (3 x NUM_ACT_SET + 1), actuation currents + inserted length; ( + 3 ), tip force (in spatial coordinates) - only for ACTUATION_AND_TIP_FORCE case
	//	the first 3+3=6 rows of the output matrix gives the Hybrid Manipulator Jacobian
	//          see Murray, Li, Sastry Intro to Robotics Book for the definition of the Manipulator Jacobian, and the Hybrid Manipulator Jacobian
	//              and how it is different from the Jacobian of the FK map
	//
	// if mode == FKFreeJacobianType::ACTUATION_ONLY
	// output matrix:
	//    (3 + 9 + 3) x (3 x NUM_ACT_SET + 1) 
	//		dp/di	:	\partial p / \partial in_x					:	partial derivatives of the tip position in spatial coordinates w.r.t. inputs
	//      ws		:	( ( \partial R / \partial in_x ) * R^T )v	:	term corresponding to the spatial angular velocity of the tip
	//		du0/di	:	\partial u0 / \partial in_x					:	partial derivative of curvature at the base of the catheter w.r.t. inputs
	// if mode == FKFreeJacobianType::ACTUATION_AND_TIP_FORCE
	// output matrix:
	//    (3 + 9 + 3) x (3 x NUM_ACT_SET + 1 + 3) 
	//		[	dp/di	dp/dftip
	//			ws/di	ws/dftip
	//			du0/di	du0/dftip   ]
	//

//
// This is the API for the brute force numerical calculation of the Forward Kinematics Jacobian 
//    i.e., the contact Jacobian in calculated by numerical differentiation of the contact FK model, 
//       not by employing the psedoinverse of the Jacobian of the free space FK
//    this function is provided for debugging purposes - should not be called independently as it is exremely slow
//    it works for bothFREE_TIP and FIXED_TIP modes, but does not calculate jacobian with respect to ftip in FREE_TIP mode
//      (note that this function returns the hybrid manipulator Jacobian + extras, not the Jacobian of the FK map)
//
template <typename adType>
MatrixXd CRM_FKJacobian_BruteForce(
	const adTypeVector<adType> in_x,	// Forward Kinematics input vector (see in_x of CRM_ForwardKinematics above)
	adTypeVector<adType>& out_y,		// Forward Kinematics output vector (see out_y of CRM_ForwardKinematics above)
	CRMForwardKinematicsData<adType> in_Params,	// Container for parameter data to be sent to Forward Kinematics
	int& localmin);						// returns 0 if successful; !=0 if algorithms is stuck at a local minimum, or cannot make further progress 
										//		(returned for forward kinematics calculation at in_x)
	// output matrix:
	//    (3+9+3) x (3 x NUM_ACT_SET + 1) matrix if Params.ContactMode == ContactModeType::FREE_TIP
	//    (3+9+3+3) x (3 x NUM_ACT_SET + 1) matrix if Params.ContactMode == ContactModeType::FIXED_TIP		
	//    outputs are packed dp/di[0..2]; ws[0..2]; du0/di[0..2] (; dftip/di[0..2])
	//		dpdi	:	\partial p / \partial in_x			 :	partial derivatives of the tip position in spatial coordinates w.r.t. inputs
	//      ws		:	( \partial R / \partial in_x ) * R^T :	term corresponding to the spatial angular velocity of the tip
	//		du0/di	:	\partial u0 / \partial in_x			 :	partial derivative of curvature at the base of the catheter w.r.t. inputs
	//		dftip/di:	\partial ftip / partial in_x		 :	partial derivative of the tip force in spatial coordinates w.r.t. inputs (if Params.ContactMode == ContactModeType::FIXED_TIP)
	//	  the first 3+3=6 rows of the output matrix gives the Hybrid Manipulator Jacobian
	//          see Murray, Li, Sastry Intro to Robotics Book for the definition of the Manipulator Jacobian, and the Hybrid Manipulator Jacobian
	//              and how it is different from the Jacobian of the FK map

//
// This function calculates the Forward Kinematics Jacobian by numerically taking derivatives of IVP
//    this function is provided purely for debugging purposes - should not be used or called independently
//
template <typename adType>
MatrixXd CRM_FK_Jacobian_From_IVP_Numerical(const adTypeVector<adType> in_x, const adTypeVector<adType>& in_FKouty, CRMForwardKinematicsData<adType> in_Params);
// helper function for CRM_FK_Jacobian_From_IVP_Numerical
template<typename adType>
VectorXd CRMSolverIVPWrapper(VectorXd in, CRMShootingMethodParams<adType> in_Params);


#endif

//
// these are internally called by TrustRegionDogleg within CRM_ForwardKinematics_Contact()
//
template <typename adType>
struct  CRMContactEquationParams : CRMForwardKinematicsData<adType> {
	adType	Actuation[3 * NUM_ACT_SET + 1];
};

template<typename adType>
void CRM_Contact_Equation(adType in_x[], adType out_y[], CRMContactEquationParams<adType> Params);

template<typename adType>
void CRM_Contact_Equation_AnalyticalJac(adType in_x[], adType out_y[], adType out_fjac[], CRMContactEquationParams<adType> Params);