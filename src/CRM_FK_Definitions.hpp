#pragma once
#include <cmath>

#ifdef DO_NOT_USE_EIGEN
// EIGEN FREE VERSION OF CRM_Forward_Kinematics

template <typename adType>
int CRM_ForwardKinematics(adType in_x[], adType out_y[], CRMForwardKinematicsData<adType> Params) {

	if (Params.ContactMode == ContactModeType::FREE_TIP)
		return CRM_ForwardKinematics_FreeSpace<adType>(in_x, out_y, Params);
	else // (Params.ContactMode == ContactModeType::FIXED_TIP)
		return CRM_ForwardKinematics_Contact<adType>(in_x, out_y, Params);

}


template <typename adType>
int CRM_ForwardKinematics_FreeSpace(adType in_x[], adType out_y[], CRMForwardKinematicsData<adType> Params) {

	ContactModeType ContactMode = Params.ContactMode;
	int X_Dim = NUM_ACT_SET * 3 + 1;	// Dimension of the Input 3 x NUM_ACT_SET + 1 for inserted length
	int Y_Dim;							// Dimension of the Output
	if (ContactMode == ContactModeType::FREE_TIP) {
		Y_Dim = 3 + 9 + 3;					// tip position + R + u_0
	}
	else { // FIXED_TIP
		Y_Dim = 3 + 9 + 3 + 3;				//   ... + tip force
	}

	adType ActuationCurrents[NUM_ACT_SET][3];
	adType InsertedLength;
	adType u0_calc[3];
	adType ftip_calc[3];
	adType x0[NUM_STATES];
	adType xf[NUM_STATES];
	adType residual[3];
	CRMShootingMethodParams<adType> BVPParams;
	int localmin;

	for (int i = 0; i < NUM_ACT_SET; i++)	for (int j = 0; j < 3; j++)	ActuationCurrents[i][j] = in_x[i * 3 + j];
	InsertedLength = in_x[NUM_ACT_SET * 3];

	CRMConstructShootingMethodParamSet(*(Params.CathParams), *(Params.CathConfig), InsertedLength, ActuationCurrents, ContactMode, Params.TipConstraintPoint, Params.TipForce, Params.IntegrationStepSize, BVPParams);

	CRMShootingMethodBVP(BVPParams, Params.u0_initialguess, Params.ftip_initialguess, u0_calc, ftip_calc, localmin);

	for (int i = 0; i < 3; i++) x0[i] = u0_calc[i];
	for (int i = 0; i < 9; i++) x0[i + 3] = Params.CathConfig->R0[i];
	for (int i = 0; i < 3; i++) x0[i + 3 + 9] = Params.CathConfig->p0[i];

	// Cosserat Rod Model - Solve the Initial Value Problem to calculate the shape of the catheter 
	CRMSolverIVP(BVPParams, x0, ftip_calc, Params.FinalValueOnly, xf, residual, *(Params.ReportedMarkerPos));

	for (int i = 0; i < 3; i++) {
		out_y[i] = xf[i];						// p
		out_y[3 + 9 + i] = u0_calc[i];			// u0
		if (ContactMode == ContactModeType::FIXED_TIP) out_y[3 + 9 + 3 + i] = ftip_calc[i];
	}
	for (int i = 0; i < 9; i++) {
		out_y[3 + i] = xf[3 + i];				// R
	}

	return (localmin);
}


template <typename adType>
int CRM_ForwardKinematics_Contact(adType in_x[], adType out_y[], CRMForwardKinematicsData<adType> Params) {

	int localmin;
	int X_Dim = NUM_ACT_SET * 3 + 1;	// Dimension of the Input 3 x NUM_ACT_SET + 1 insertion
	int Y_Dim = 3 + 9 + 3 + 3;			// Dimension of the Output: // tip position + R + u_0 + tip force

	//
	// We will solve the nonlinear equations constraining the tip to the given tip position by varying the tip force
	//

	CRMContactEquationParams<adType> NLEParams;
	NLEParams.CRMForwardKinematicsData::operator= (Params);				// CRMForwardKinematicsData portion of NLEParams comes from Params
	for (int i = 0; i < X_Dim; i++) NLEParams.Actuation[i] = in_x[i];	// Actutation part of the NLEParams comes from in_x
	NLEParams.ContactMode = ContactModeType::FREE_TIP;					// We will use this to solve FREE_TIP model

	const double FSCALE_INV = 1.0 / IVALUE_SCALE_F;
	double initialguessscaled[3];
	adType returnedparamscaled[3];
	for (int i = 0; i < 3; i++) initialguessscaled[i] = FSCALE_INV * Params.ftip_initialguess[i];

#if defined( CONTACT_TRUSTREGION )
	adType x[3]; // we will create a new variable here and not use initialguessscaled since truss-region-dogleg algorithm uses the same variable for both input and output
	adType residual[3];
	int info;
	double tol = TRUSTREGION_TOLERANCE;
	for (int i = 0; i < 3; i++) x[i] = initialguessscaled[i];

	TrustRegionDogleg(CRM_Contact_Equation<adType>, 3, x, residual, tol, info, NLEParams);

	for (int i = 0; i < 3; i++) returnedparamscaled[i] = x[i];
	localmin = (info == 1) ? 0 : (info - 1);
#else // undefined
	exit(1);
#endif

	//
	// Once we know the contact force, we can run the FK again in FREE_TIP mode using the contact force as the tip force
	//
	for (int i = 0; i < 3; i++) NLEParams.TipForce[i] = IVALUE_SCALE_F * returnedparamscaled[i];
	CRM_ForwardKinematics<adType>(in_x, out_y, NLEParams);
	for (int i = 0; i < 3; i++) out_y[15 + i] = IVALUE_SCALE_F * returnedparamscaled[i];  // tip force for output comes from the solution above

	return localmin;  // we are returning localmin that is given by the original NL equations solution
}


#else // ifdef DO_NOT_USE_EIGEN
// Version of CRM Forward Kinematics that uses Eigen

template <typename adType>
int CRM_ForwardKinematics(adType in_x[], adType out_y[], CRMForwardKinematicsData<adType> Params) {

	int localmin;
	int X_Dim = NUM_ACT_SET * 3 + 1;	// Dimension of the Input 3 x NUM_ACT_SET + 1 insertion
	int Y_Dim;						// Dimension of the Output
	if (Params.ContactMode == ContactModeType::FREE_TIP) {
		Y_Dim = 3 + 9 + 3;					// tip position + R + u_0
	}
	else { // FIXED_TIP
		Y_Dim = 3 + 9 + 3 + 3;				//   ... + tip force
	}
	adTypeVector<adType> x(X_Dim), y(Y_Dim);

	for (int i = 0; i < X_Dim; i++) x(i) = in_x[i];

	y = CRM_ForwardKinematics<adType>(x, Params, localmin);

	for (int i = 0; i < Y_Dim; i++) out_y[i] = y(i);

	return localmin;
}


template <typename adType>
adTypeVector<adType> CRM_ForwardKinematics(const adTypeVector<adType> in_x, CRMForwardKinematicsData<adType> Params, int& localmin) {

	if (Params.ContactMode == ContactModeType::FREE_TIP)
		return CRM_ForwardKinematics_FreeSpace<adType>(in_x, Params, localmin);
	else // (Params.ContactMode == ContactModeType::FIXED_TIP)
		return CRM_ForwardKinematics_Contact<adType>(in_x, Params, localmin);

}


template <typename adType>
int CRM_ForwardKinematics_FreeSpace(adType in_x[], adType out_y[], CRMForwardKinematicsData<adType> Params) {

	int localmin;
	int X_Dim = NUM_ACT_SET * 3 + 1;	// Dimension of the Input 3 x NUM_ACT_SET + 1 insertion
	int Y_Dim;						// Dimension of the Output
	if (Params.ContactMode == ContactModeType::FREE_TIP) {
		Y_Dim = 3 + 9 + 3;					// tip position + R + u_0
	}
	else { // FIXED_TIP
		Y_Dim = 3 + 9 + 3 + 3;				//   ... + tip force
	}
	adTypeVector<adType> x(X_Dim), y(Y_Dim);

	for (int i = 0; i < X_Dim; i++) x(i) = in_x[i];

	y = CRM_ForwardKinematics_FreeSpace<adType>(x, Params, localmin);

	for (int i = 0; i < Y_Dim; i++) out_y[i] = y(i);

	return localmin;
}


template <typename adType>
int CRM_ForwardKinematics_Contact(adType in_x[], adType out_y[], CRMForwardKinematicsData<adType> Params) {

	int localmin;
	int X_Dim = NUM_ACT_SET * 3 + 1;	// Dimension of the Input 3 x NUM_ACT_SET + 1 insertion
	int Y_Dim = 3 + 9 + 3 + 3;			// Dimension of the Output // tip position + R + u_0 + tip force
	adTypeVector<adType> x(X_Dim), y(Y_Dim);

	for (int i = 0; i < X_Dim; i++) x(i) = in_x[i];

	y = CRM_ForwardKinematics_Contact<adType>(x, Params, localmin);

	for (int i = 0; i < Y_Dim; i++) out_y[i] = y(i);

	return localmin;
}



template <typename adType>
adTypeVector<adType> CRM_ForwardKinematics_FreeSpace(const adTypeVector<adType> in_x, CRMForwardKinematicsData<adType> Params, int& localmin) {
	
	int X_Dim = NUM_ACT_SET * 3 + 1;	// Dimension of the Input 3 x NUM_ACT_SET + 1 for inserted length
	int Y_Dim;							// Dimension of the Output
	if (Params.ContactMode == ContactModeType::FREE_TIP) {
		Y_Dim = 3 + 9 + 3;					// tip position + R + u_0
	}
	else { // FIXED_TIP
		Y_Dim = 3 + 9 + 3 + 3;				//   ... + tip force
	}
	adTypeVector<adType> out_y(Y_Dim);
	
	adType ActuationCurrents[NUM_ACT_SET][3];
	adType InsertedLength;
	adType u0_calc[3];
	adType ftip_calc[3];
	adType x0[NUM_STATES];
	adType xf[NUM_STATES];
	adType residual[3];
	CRMShootingMethodParams<adType> BVPParams;
	auto& CathParams = *(Params.CathParams);
	auto& CathConfig = *(Params.CathConfig);

	for (int i = 0; i < NUM_ACT_SET; i++)	for (int j = 0; j < 3; j++)	ActuationCurrents[i][j] = in_x(i * 3 + j);
	InsertedLength = in_x(NUM_ACT_SET * 3);

	CRMConstructShootingMethodParamSet(CathParams, CathConfig, InsertedLength, ActuationCurrents, Params.ContactMode, Params.TipConstraintPoint, Params.TipForce, Params.IntegrationStepSize, BVPParams);

	CRMShootingMethodBVP(BVPParams, Params.u0_initialguess, Params.ftip_initialguess, u0_calc, ftip_calc, localmin);

	for (int i = 0; i < 3; i++) x0[i] = u0_calc[i];
	for (int i = 0; i < 9; i++) x0[i + 3] = Params.CathConfig->R0[i];
	for (int i = 0; i < 3; i++) x0[i + 3 + 9] = Params.CathConfig->p0[i];

	// Cosserat Rod Model - Solve the Initial Value Problem to calculate the shape of the catheter 
	CRMSolverIVP(BVPParams, x0, ftip_calc, Params.FinalValueOnly, xf, residual, *(Params.ReportedMarkerPos));

	for (int i = 0; i < 3; i++) {
		out_y(i) = xf[i];						// p
		out_y(3 + 9 + i) = u0_calc[i];			// u0
		if (Params.ContactMode == ContactModeType::FIXED_TIP) out_y(3 + 9 + 3 + i) = ftip_calc[i];
	}
	for (int i = 0; i < 9; i++) {
		out_y(3 + i) = xf[3 + i];				// R
	}

	return (out_y);
}


template <typename adType>
adTypeVector<adType> CRM_ForwardKinematics_Contact(const adTypeVector<adType> in_x, CRMForwardKinematicsData<adType> Params, int& localmin) {

	int X_Dim = NUM_ACT_SET * 3 + 1;	// Dimension of the Input 3 x NUM_ACT_SET + 1 insertion
	int Y_Dim = 3 + 9 + 3 + 3;			// Dimension of the Output: // tip position + R + u_0 + tip force
	adTypeVector<adType> out_y(Y_Dim);

	//
	// We will solve the nonlinear equations constraining the tip to the given tip position by varying the tip force
	//

	CRMContactEquationParams<adType> NLEParams;
	NLEParams.CRMForwardKinematicsData<adType>::operator= (Params);				// CRMForwardKinematicsData portion of NLEParams comes from Params
	for (int i = 0; i < X_Dim; i++) NLEParams.Actuation[i] = in_x(i);	// Actutation part of the NLEParams comes from in_x
	NLEParams.ContactMode = ContactModeType::FREE_TIP;					// We will use this to solve FREE_TIP model

	const double FSCALE_INV = 1.0 / IVALUE_SCALE_F;
	double initialguessscaled[3];
	adType returnedparamscaled[3];
	for (int i = 0; i < 3; i++) initialguessscaled[i] = FSCALE_INV * Params.ftip_initialguess[i];

#if defined( CONTACT_TRUSTREGION )
	adType x[3]; // we will create a new variable here and not use initialguessscaled since truss-region-dogleg algorithm uses the same variable for both input and output
	adType residual[3];
	int info;
	double tol = TRUSTREGION_TOLERANCE;
	for (int i = 0; i < 3; i++) x[i] = initialguessscaled[i];

#if defined  (CONTACT_TRUSTREGION_ANALYTICALJAC)
	TrustRegionDogleg_GivenJacobian<adType, CRMContactEquationParams<adType>>(CRM_Contact_Equation<adType>, CRM_Contact_Equation_AnalyticalJac<adType>, 3, x, residual, tol, info, NLEParams);
#else 	
	TrustRegionDogleg(CRM_Contact_Equation<adType>, 3, x, residual, tol, info, NLEParams);
#endif

	for (int i = 0; i < 3; i++) returnedparamscaled[i] = x[i];
	localmin = (info == 1) ? 0 : (info - 1);
#else // undefined
	exit(1);
#endif

	//
	// Once we know the contact force, we can run the FK again in FREE_TIP mode using the contact force as the tip force
	//
	int ignore;  // we are returning localmin that is given by the original NL equations solution
	adTypeVector<adType> temp_out;
	for (int i = 0; i < 3; i++) NLEParams.TipForce[i] = IVALUE_SCALE_F * returnedparamscaled[i];
	temp_out = CRM_ForwardKinematics<adType>(in_x, NLEParams, ignore);
	out_y.head(15) = temp_out;
	for (int i = 0; i < 3; i++) out_y(15 + i) = IVALUE_SCALE_F * returnedparamscaled[i];  // tip force for output comes from the solution above

	return out_y;  
}


template <typename adType>
MatrixXd CRM_FKJacobian_Numerical(const adTypeVector<adType> in_x, adTypeVector<adType>& out_y, CRMForwardKinematicsData<adType> in_Params, int& localmin) {

	if (in_Params.ContactMode == ContactModeType::FREE_TIP) {
		return ( CRM_FKJacobian_FreeSpace(in_x, out_y, in_Params, FKFreeJacobianType::ACTUATION_ONLY, localmin) );
	} 
	// else ContactModeType::FIXED_TIP
	CRMForwardKinematicsData<adType>Params = in_Params;					// create a local copy
	Params.FinalValueOnly = true;
	out_y = CRM_ForwardKinematics<adType>(in_x, Params, localmin);	// calculate forward kinematics at x_0 for contact model
	for (int i = 0; i < 3; i++) Params.TipForce[i] = out_y(3 + 9 + 3 + i);  // we need to use the calculated contact force as the tip force for free-space model
	int templocalmin;  // dummy variable
	adTypeVector<adType> tempvector;  // dummy variable
	MatrixXd JBVP = CRM_FKJacobian_FreeSpace(in_x, tempvector, Params, FKFreeJacobianType::ACTUATION_AND_TIP_FORCE, templocalmin);
	MatrixXd JBVPp = JBVP.topRows(3);
	MatrixXd JBVPp_z = JBVPp.leftCols(JBVPp.cols() - 3);
	MatrixXd JBVPp_ftip = JBVPp.rightCols(3);
	MatrixXd J = -JBVPp_ftip.completeOrthogonalDecomposition().pseudoInverse() * JBVPp_z;  // Is this the best option in Eigen ???
	return J;
}


template <typename adType>
MatrixXd CRM_FKJacobian_BruteForce(const adTypeVector<adType> in_x, adTypeVector<adType>& out_y, CRMForwardKinematicsData<adType> in_Params, int& localmin) {

	CRMForwardKinematicsData<adType>Params = in_Params;					// create a local copy
	Params.FinalValueOnly = true;
	adTypeVector<adType> x_0 = in_x;									// create a local copy
	adTypeVector<adType> f_0 = CRM_ForwardKinematics<adType>(x_0, Params, localmin);	// calculate forward kinematics at x_0
	auto f_size = f_0.size();
	auto var_size = x_0.size();
	MatrixXd Df_(f_size, var_size);
	adTypeVector<adType> x_h(var_size), f_h(f_size);
	int lmin;

	// we will use the u0 calculated for FK at x_0 as the initial guess for all of the subsequent FK calculations
	for (int i = 0; i < 3; i++) {
		Params.u0_initialguess[i] = dVal(f_0(3 + 9 + i));
	}
	// in_Params.ContactMode == ContactModeType::FIXED_TIP, we will aslo use the ftip calculated for FK at x_0 as the initial guess for all of the subsequent FK calculations
	if (Params.ContactMode == ContactModeType::FIXED_TIP) {
		for (int i = 0; i < 3; i++) {
			Params.ftip_initialguess[i] = dVal(f_0(3 + 9 + 3 + i));
		}
	}
	f_0 = CRM_ForwardKinematics<adType>(x_0, Params, localmin);	// calculate forward kinematics at x_0 with the new initial conditions for numerical stability

	double h = NUM_JACOBIAN_CURRENT_STEPSIZE;
	double hinv = 1 / NUM_JACOBIAN_CURRENT_STEPSIZE;
	for (int j = 0; j < var_size; ++j) {
		if (j == var_size - 1) {  // last variable is the insertion length
			h = NUM_JACOBIAN_INSERTIONLENGTH_STEPSIZE;
			hinv = 1 / NUM_JACOBIAN_INSERTIONLENGTH_STEPSIZE;
		}
		x_h = x_0;  x_h(j) = x_0(j) + h;
		f_h = CRM_ForwardKinematics<adType>(x_h, Params, lmin);
		Df_.col(j) = (f_h - f_0) * hinv;
		// for (int k = 0; k < f_size; ++k) {
		//	Df_(k, j) = (f_h(k) - f_0(k)) / h;
		//}
	}
	out_y = f_0;

	// while we are at it, let's calculate the spatial angular velocity from \dot{R}
	int JRow;
	if (Params.ContactMode == ContactModeType::FREE_TIP) {
		JRow = 3 + 3 + 3;				// pdot + ws + u_0dot
	}
	else { // FIXED_TIP
		JRow = 3 + 3 + 3 + 3;			//   ... + tip force dot
	}
	MatrixXd T(JRow, f_size);
	MatrixXd J(JRow, var_size);
	T.setZero();
	T.block(0, 0, 3, 3).setIdentity();
	T.block(6, 12, 3, 3).setIdentity();
	T.block(3, 9, 1, 3) << dVal(f_0(6)), dVal(f_0(7)), dVal(f_0(8));
	T.block(4, 3, 1, 3) << dVal(f_0(9)), dVal(f_0(10)), dVal(f_0(11));
	T.block(5, 6, 1, 3) << dVal(f_0(3)), dVal(f_0(4)), dVal(f_0(5));
	if (Params.ContactMode == ContactModeType::FIXED_TIP) 	T.block(9, 15, 3, 3).setIdentity();
	J = T * Df_;

	return J;

}


template <typename adType>
MatrixXd CRM_FKJacobian_FreeSpace(const adTypeVector<adType> in_x, adTypeVector<adType>& out_y, CRMForwardKinematicsData<adType> in_Params, FKFreeJacobianType mode, int& localmin) {
	
	adTypeVector<adType> x_0 = in_x;									// create a local copy
	CRMForwardKinematicsData<adType>Params = in_Params;					// create a local copy
	Params.FinalValueOnly = true;										// we don't need marker locations
	Params.ContactMode = ContactModeType::FREE_TIP;						// we will use Free Space model

	adTypeVector<adType> f_0 = CRM_ForwardKinematics<adType>(x_0, Params, localmin);	// calculate forward kinematics at x_0

	double h = NUM_JACOBIAN_CURRENT_STEPSIZE;
	double hinv = 1/ NUM_JACOBIAN_CURRENT_STEPSIZE;
	auto f_size = f_0.size();
	auto var_size = x_0.size();
	auto col_size = (mode == FKFreeJacobianType::ACTUATION_ONLY) ? var_size : (var_size+3);   // we need extra 3 columns for f_tip
	MatrixXd Df_(f_size, col_size);
	adTypeVector<adType> x_h(var_size), f_h(f_size);
	int lmin;

	// we will use the u0 calculated for FK at x_0 as the initial guess for all of the subsequent FK calculations
	for (int i = 0; i < 3; i++) {
		Params.u0_initialguess[i] = dVal(f_0(12 + i));
	}
	f_0 = CRM_ForwardKinematics<adType>(x_0, Params, localmin);	// calculate forward kinematics at x_0 with the new initial conditions for numerical stability

	for (int j = 0; j < var_size; ++j) {
		if (j == var_size - 1) {  // last variable is the insertion length
			h = NUM_JACOBIAN_INSERTIONLENGTH_STEPSIZE;
			hinv = 1 / NUM_JACOBIAN_INSERTIONLENGTH_STEPSIZE;
		}
		x_h = x_0;  x_h(j) = x_0(j) + h;
		f_h = CRM_ForwardKinematics<adType>(x_h, Params, lmin);
		Df_.col(j) = (f_h - f_0) * hinv; 
		// for (int k = 0; k < f_size; ++k) {
		//	Df_(k, j) = (f_h(k) - f_0(k)) / h;
		//}
	}
	adType dummy;  // temporary storage
	h = NUM_JACOBIAN_TIPFORCE_STEPSIZE;
	hinv = 1 / NUM_JACOBIAN_TIPFORCE_STEPSIZE;
	for (int j = var_size, k=0; j < col_size; j++, k++) {
		dummy = Params.TipForce[k];
		Params.TipForce[k] += h;
		f_h = CRM_ForwardKinematics<adType>(x_0, Params, lmin);
		Df_.col(j) = (f_h - f_0) * hinv;
		Params.TipForce[k] = dummy;
	}
	out_y = f_0;

	// while we are at it, let's calculate the spatial angular velocity from \dot{R}
	int JRow= 3 + 3 + 3;				// pdot + ws + u_0dot
	MatrixXd T(JRow, f_size);
	MatrixXd J(JRow, col_size);
	T.setZero();
	T.block(0, 0, 3, 3).setIdentity();
	T.block(6, 12, 3, 3).setIdentity();
	T.block(3, 9, 1, 3) << dVal(f_0(6)), dVal(f_0(7)),	dVal(f_0(8));
	T.block(4, 3, 1, 3) << dVal(f_0(9)), dVal(f_0(10)), dVal(f_0(11));
	T.block(5, 6, 1, 3) << dVal(f_0(3)), dVal(f_0(4)),	dVal(f_0(5));
	J = T * Df_;

	return J;

}


template <typename adType>
MatrixXd CRM_FKJacobian_Analytical(const adTypeVector<adType> in_x, const adTypeVector<adType>& in_FKouty, CRMForwardKinematicsData<adType> in_Params) {


	CRMShootingMethodParams<adType> BVPParams;
	auto& CathParams = *(in_Params.CathParams);
	auto& CathConfig = *(in_Params.CathConfig);
	adType ActuationCurrents[NUM_ACT_SET][3];
	adType InsertedLength;
	for (int i = 0; i < NUM_ACT_SET; i++)	for (int j = 0; j < 3; j++)	ActuationCurrents[i][j] = in_x(i * 3 + j);
	InsertedLength = in_x(NUM_ACT_SET * 3);

	CRMConstructShootingMethodParamSet(CathParams, CathConfig, InsertedLength, ActuationCurrents, in_Params.ContactMode, in_Params.TipConstraintPoint, in_Params.TipForce, in_Params.IntegrationStepSize, BVPParams);

	adType u0_calc[3];
	adType ftip_calc[3];

	for (int i = 0; i < 3; i++) u0_calc[i] = in_FKouty(i + 3 + 9);
	if (in_Params.ContactMode == ContactModeType::FREE_TIP) {
		for (int i = 0; i < 3; i++) ftip_calc[i] = in_Params.TipForce[i];
	}
	else {
		for (int i = 0; i < 3; i++) ftip_calc[i] = in_FKouty(i + 3 + 9 + 3);
	}

	adType xf[NUM_STATES];
	adType residual[3];

	auto [JBVP_p_z, JBVP_ws_z, JBVP_p_ft, JBVP_ws_ft, Jft_z] = CRMSolverIVPJacobian(BVPParams, u0_calc, ftip_calc, false, xf, residual);

	if (in_Params.ContactMode == ContactModeType::FREE_TIP) {
		MatrixXd JBVP_pws_zft(6, JBVP_p_z.cols()+3);
		JBVP_pws_zft.template topLeftCorner<3, 3 * NUM_ACT_SET + 1>() = JBVP_p_z;
		JBVP_pws_zft.template bottomLeftCorner<3, 3 * NUM_ACT_SET + 1 >() = JBVP_ws_z;
		JBVP_pws_zft.template topRightCorner<3,3>() = JBVP_p_ft;
		JBVP_pws_zft.template bottomRightCorner<3,3>() = JBVP_ws_ft;
		return JBVP_pws_zft;
	}
	else { // (in_Params.ContactMode == ContactModeType::FIXED_TIP)
		return Jft_z;
	}

}


template<typename adType>
VectorXd CRMSolverIVPWrapper(VectorXd in, CRMShootingMethodParams<adType> in_Params) {

	CRMShootingMethodParams<adType> Params = in_Params;
	adType u0_calc[3];
	adType ftip_calc[3];
	adType actvect[3];

	for (int i = 0; i < 3; i++) Params.p0[i] = in(i);
	for (int i = 0; i < 9; i++) Params.R0[i] = in(i + 3);
	for (int i = 0; i < 3; i++) u0_calc[i] = in(i + 3 + 9);
	for (int i = 0; i < 3; i++) ftip_calc[i] = in(i + 3 + 9 + 3);
	for (int i = 0; i < NUM_ACT_SET; i++) {
		for (int j = 0; j < 3; j++) actvect[j] = in(i * 3 + j + 3 + 9 + 3 + 3);
		mMult_AB<3, 3, 1>(Params.CoilAlignmentTurnAreaMatrix[i], actvect, Params.MagMoment[i]);
	}
	Params.Li = in(in.rows()-1);


	adType xf[NUM_STATES];
	__EVT<adType, 15> XF(xf);
	adType residual[3];
	double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]; // dummy

	CRMSolverIVP(Params, u0_calc, ftip_calc, false, xf, residual, out_p_atLocMarkers);

	VectorXd result = XF;
	return result;
}

template <typename adType>
MatrixXd CRM_FK_Jacobian_From_IVP_Numerical(const adTypeVector<adType> in_x, const adTypeVector<adType>& in_FKouty, CRMForwardKinematicsData<adType> in_Params) {


	CRMShootingMethodParams<adType> BVPParams;
	auto& CathParams = *(in_Params.CathParams);
	auto& CathConfig = *(in_Params.CathConfig);
	adType ActuationCurrents[NUM_ACT_SET][3];
	adType InsertedLength;
	for (int i = 0; i < NUM_ACT_SET; i++)	for (int j = 0; j < 3; j++)	ActuationCurrents[i][j] = in_x(i * 3 + j);
	InsertedLength = in_x(NUM_ACT_SET * 3);

	CRMConstructShootingMethodParamSet(CathParams, CathConfig, InsertedLength, ActuationCurrents, in_Params.ContactMode, in_Params.TipConstraintPoint, in_Params.TipForce, in_Params.IntegrationStepSize, BVPParams);

	Matrix<adType, 3 + 9 + 3 + 3 + 3 * NUM_ACT_SET + 1, 1> x0, xh;
	Matrix<adType, 3 + 9 + 3, 1> y, y0;
	Matrix<adType, 3 + 9 + 3, 3 + 9 + 3 + 3 * NUM_ACT_SET + 1 + 3> J_reversezc;

	if (in_Params.ContactMode == ContactModeType::FREE_TIP) {
		x0.template head<3 + 9 + 3>() = in_FKouty;
		for (int i = 0; i < 3; i++) x0(3 + 9 + 3 + i) = in_Params.TipForce[i];
	}
	else {
		x0.template head<3 + 9 + 3 + 3>() = in_FKouty;

	}
	x0.template segment<3 * NUM_ACT_SET + 1>(3 + 9 + 3 + 3) = in_x;

	y0 = CRMSolverIVPWrapper(x0, BVPParams);
	double h;
	for (int i = 0; i < x0.rows(); i++) {
		xh = x0;
		if (i < 3) h = NUM_JACOBIAN_BASEPOSITION_STEPSIZE;
		else if (i < 3 + 9) h = NUM_JACOBIAN_BASEROTATION_STEPSIZE;
		else if (i < 3 + 9 + 3) h = NUM_JACOBIAN_BASECURVATURE_STEPSIZE;
		else if (i < 3 + 9 + 3 + 3) h = NUM_JACOBIAN_TIPFORCE_STEPSIZE;
		else if (i < 3 + 9 + 3 + 3 + 3 * NUM_ACT_SET) h = NUM_JACOBIAN_CURRENT_STEPSIZE;
		else h = NUM_JACOBIAN_INSERTIONLENGTH_STEPSIZE;
		xh(i) += h;
		y = CRMSolverIVPWrapper(xh, BVPParams);
		J_reversezc.col(i) = (y - y0) / h;
	}

	Matrix<adType, 3 + 9 + 3, 3 + 9 + 3 + 3 * NUM_ACT_SET + 1 + 3> JIVP_x_phi;
	JIVP_x_phi.template leftCols<3 + 9 + 3>() = J_reversezc.template leftCols<3 + 9 + 3>();
	// reorder ft
	JIVP_x_phi.template rightCols<3>() = J_reversezc.template middleCols<3>(3 + 9 + 3);
	// we will need to flip the order of actuators in zc, since IVP calculations use proximal to distal ordering
	for (int i = 0; i < NUM_ACT_SET; i++) JIVP_x_phi.template middleCols<3>(3 + 9 + 3 + 3 * i) = J_reversezc.template middleCols<3>(3 + 9 + 3 + 3 + 3 * ((NUM_ACT_SET - 1) - i));
	JIVP_x_phi.template middleCols<1>(3 + 9 + 3 + 3 * NUM_ACT_SET) = J_reversezc.template rightCols<1>();

	// let's calculate the FK Jacobians
	// for simplicity, define an alias
	constexpr unsigned int Cs = CURRENT_ACT_VECTOR_DIM;

	// while we are at it, let's calculate the spatial angular velocity from \dot{R}
	Matrix<adType, 9, 1> _R = in_FKouty.template segment<9>(3);
	Matrix<adType, 3, 9, RowMajor> T;
	T.setZero();
	T.block(0, 6, 1, 3) << dVal(_R(3)), dVal(_R(4)), dVal(_R(5));
	T.block(1, 0, 1, 3) << dVal(_R(6)), dVal(_R(7)), dVal(_R(8));
	T.block(2, 3, 1, 3) << dVal(_R(0)), dVal(_R(1)), dVal(_R(2));

	Matrix<adType, 3, 3, RowMajor> JIVP_p_u0 = JIVP_x_phi.template block<3, 3>(0, 3 + 9);
	Matrix<adType, 3, Cs + 1, RowMajor> JIVP_p_z = JIVP_x_phi.template block<3, Cs + 1>(0, 3 + 9 + 3);
	Matrix<adType, 3, 3, RowMajor> JIVP_p_ft = JIVP_x_phi.template block<3, 3>(0, 3 + 9 + 3 + Cs + 1);
	Matrix<adType, 3, 3, RowMajor> JIVP_u_u0 = JIVP_x_phi.template block<3, 3>(3 + 9, 3 + 9);
	Matrix<adType, 3, Cs + 1, RowMajor> JIVP_u_z = JIVP_x_phi.template block<3, Cs + 1>(3 + 9, 3 + 9 + 3);
	Matrix<adType, 3, 3, RowMajor> JIVP_u_ft = JIVP_x_phi.template block<3, 3>(3 + 9, 3 + 9 + 3 + Cs + 1);
	Matrix<adType, 3, 3, RowMajor> JIVP_ws_u0 = T * JIVP_x_phi.template block<9, 3>(3, 3 + 9);// _R_u0;
	Matrix<adType, 3, 3, RowMajor> JIVP_ws_ft = T * JIVP_x_phi.template block<9, 3>(3, 3 + 9 + 3 + Cs + 1);// _R_ft;
	Matrix<adType, 3, Cs + 1, RowMajor> JIVP_ws_z = T * JIVP_x_phi.template block<9, Cs+1>(3, 3 + 9 +3);// _R_z;

	MatrixXd JIVP_u_u0_pinv = JIVP_u_u0.completeOrthogonalDecomposition().pseudoInverse();
	MatrixXd JBVP_p_z = JIVP_p_z - JIVP_p_u0 * JIVP_u_u0_pinv * JIVP_u_z;
	MatrixXd JBVP_ws_z = JIVP_ws_z - JIVP_ws_u0 * JIVP_u_u0_pinv * JIVP_u_z;

	MatrixXd JBVP_p_ft = JIVP_p_ft - JIVP_p_u0 * JIVP_u_u0_pinv * JIVP_u_ft;
	MatrixXd Jft_z = -JBVP_p_ft.completeOrthogonalDecomposition().pseudoInverse() * JBVP_p_z;  // Is this the best option in Eigen ???


	if (in_Params.ContactMode == ContactModeType::FREE_TIP) {
		MatrixXd JBVP_pws_z(6, JBVP_p_z.cols());
		JBVP_pws_z.topRows(3) = JBVP_p_z;
		JBVP_pws_z.bottomRows(3) = JBVP_ws_z;
		return JBVP_pws_z;
	}
	else { // (in_Params.ContactMode == ContactModeType::FIXED_TIP)
		return Jft_z;
	}

}



#endif //#ifndef DO_NOT_USE_EIGEN


template<typename adType>
void CRM_Contact_Equation(adType in_x[], adType out_y[], CRMContactEquationParams<adType> Params) {

	for (int i = 0; i < 3; i++) Params.TipForce[i] = IVALUE_SCALE_F * in_x[i];			// we will use the in_x as the tip foce
	Params.ContactMode = ContactModeType::FREE_TIP;										//    and calculate FK using FREE_TIP
	adType fkoutput[3 + 9 + 3];		// tip position + R + u_0							//
	CRM_ForwardKinematics<adType>(Params.Actuation, fkoutput, Params);					//
	for (int i = 0; i < 3; i++) out_y[i] = RESIDUAL_SCALE_P * (fkoutput[i] - Params.TipConstraintPoint[i]);	// and return the residual

}

template<typename adType>
void CRM_Contact_Equation_AnalyticalJac(adType in_x[], adType out_y[], adType out_fjac[], CRMContactEquationParams<adType> Params) {

	// we need to solve the free space equiliibrium equation first to send to jacobian calculation 
	for (int i = 0; i < 3; i++) Params.TipForce[i] = IVALUE_SCALE_F * in_x[i];			// we will use the in_x as the tip foce
	Params.ContactMode = ContactModeType::FREE_TIP;										//    and calculate FK using FREE_TIP
	adType fkoutput[3 + 9 + 3];		// tip position + R + u_0							//
	CRM_ForwardKinematics<adType>(Params.Actuation, fkoutput, Params);					//
	for (int i = 0; i < 3; i++) out_y[i] = RESIDUAL_SCALE_P * (fkoutput[i] - Params.TipConstraintPoint[i]);	// and return the residual

	Matrix<adType, 3 * NUM_ACT_SET + 1, 1> u;
	Matrix<adType, 15, 1> y;
	for (int i = 0; i < 3 * NUM_ACT_SET + 1; i++) u(i) = Params.Actuation[i];
	for (int i = 0; i < 3 + 9 + 3; i++) y(i) = fkoutput[i];
	Matrix<adType, 6, 3 * NUM_ACT_SET + 1 + 3 > Jfull = CRM_FKJacobian_Analytical<adType>(u, y, Params);  // [p;ws] vs [zc,zl,ft]
	Matrix<adType, 3, 3> sJp = RESIDUAL_SCALE_P * Jfull.template topRightCorner<3, 3>() * IVALUE_SCALE_F;
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) out_fjac[i + j * 3] = sJp(i, j);  // Note that minpack uses Fortran style column-major ordering, not C sytle row-major ordering

}

