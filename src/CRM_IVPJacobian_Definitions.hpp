#pragma once
#include <cmath>

template <typename adType>
std::tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd> CRMSolverIVPJacobian(
		CRMShootingMethodParams<adType> in_Params,
		adType in_u0[3], adType in_ftip[3],
		bool in_FinalValueOnly,
		adType out_x_N[NUM_STATES], adType out_MomentResidual[3]) {

	adType x_0[NUM_STATES];
	for (int i = 0; i < 3; i++) x_0[i] = in_Params.p0[i];
	for (int i = 0; i < 9; i++) x_0[i + 3] = in_Params.R0[i];
	for (int i = 0; i < 3; i++) x_0[i + 3 + 9] = in_u0[i];

	CRMIVPCoreParams<adType> CoreParams;
	adType u_0[3], ftip[3];
	bool FinalValueOnly = false;

	AugmentedStateVector<adType, IVPJacobiansFull> x_N;

	CRMSolverIVP_Prep(x_0, in_Params.IntegrationStepSize,
		in_Params.Li, in_Params.dlambdainv,
		in_Params.SegEndLambdas, in_Params.LocMarkerLambdas,
		in_Params.K, in_Params.Kinv, in_Params.ustar,
		in_Params.CoilAlignmentTurnAreaMatrix,
		in_Params.MagMoment, in_Params.fcumlambda,
		in_Params.B0,
		FinalValueOnly,
		CoreParams);

	// copy to local variable
	for (int i = 0; i < 3; i++) ftip[i] = in_ftip[i];

	// We need to pass u_0 as input argument as the values in CoreParams will be overriden with the values provided in the input arguments - functionality needed for solving Boundary Value Problems (BVP)
	for (int i = 0; i < 3; i++) u_0[i] = in_u0[i];

	CRMSolverIVP_CoreWithJacobian(CoreParams, u_0, ftip, x_N, out_MomentResidual);

	// NO NEED FOR CRMSolverIVP_Return
	mCopy_AB<3>(x_N._p, out_x_N + 0);
	mCopy_AB<9>(x_N._R, out_x_N + 3);
	mCopy_AB<3>(x_N._u, out_x_N + 3 + 9);

	// for simplicity, define an alias
	constexpr unsigned int Cs = CURRENT_ACT_VECTOR_DIM;

	//using EMT33 = Eigen::Map<Matrix<double, 3, 3, RowMajor>, Eigen::Unaligned, Eigen::Stride<Cs, 1> >;

	__EMT<adType, 3, 3> JIVP_p_u0(x_N._p_u0);
	__EMT<adType, 3, 3> JIVP_ws_u0(x_N._ws_u0);
	__EMT<adType, 3, 3> JIVP_u_u0(x_N._u_u0);

	__EMT<adType, 3, 3> JIVP_p_ft(x_N._p_ft);
	__EMT<adType, 3, 3> JIVP_ws_ft(x_N._ws_ft);
	__EMT<adType, 3, 3> JIVP_u_ft(x_N._u_ft);

	Matrix<adType, 3, Cs + 1, RowMajor> JIVP_p_z;
	Matrix<adType, 3, Cs + 1, RowMajor> JIVP_ws_z;
	Matrix<adType, 3, Cs + 1, RowMajor> JIVP_u_z;
	// we will need to flip the order of actuators in zc, since IVP Core calculations use proximal to distal ordering
	//_p_z
	__EMT<adType, 3, Cs> JIVP_p_zc_reverse(x_N._p_zc);
	__EVT<adType, 3> JIVP_p_zl(x_N._p_zl);
	for (int i = 0; i < NUM_ACT_SET; i++) JIVP_p_z.template middleCols<3>(3 * i) = JIVP_p_zc_reverse.template middleCols<3>(((NUM_ACT_SET - 1) - i) * 3);
	JIVP_p_z.template rightCols<1>() = JIVP_p_zl;
	//_R_z
	__EMT<adType, 3, Cs> JIVP_ws_zc_reverse(x_N._ws_zc);
	__EVT<adType, 3> JIVP_ws_zl(x_N._ws_zl);
	for (int i = 0; i < NUM_ACT_SET; i++) JIVP_ws_z.template middleCols<3>(3 * i) = JIVP_ws_zc_reverse.template middleCols<3>(((NUM_ACT_SET - 1) - i) * 3);
	JIVP_ws_z.template rightCols<1>() = JIVP_ws_zl;
	//_u_z
	__EMT<adType, 3, Cs> JIVP_u_zc_reverse(x_N._u_zc);
	__EVT<adType, 3> JIVP_u_zl(x_N._u_zl);
	for (int i = 0; i < NUM_ACT_SET; i++) JIVP_u_z.template middleCols<3>(3 * i) = JIVP_u_zc_reverse.template middleCols<3>(((NUM_ACT_SET - 1) - i) * 3);
	JIVP_u_z.template rightCols<1>() = JIVP_u_zl;
		
	// for debugging
	//std::cout << "---- x_N ---- \n" << x_N << std::endl;

	MatrixXd JIVP_u_u0_pinv = JIVP_u_u0.completeOrthogonalDecomposition().pseudoInverse();
	MatrixXd JBVP_p_z = JIVP_p_z - JIVP_p_u0 * JIVP_u_u0_pinv * JIVP_u_z;
	MatrixXd JBVP_ws_z = JIVP_ws_z - JIVP_ws_u0 * JIVP_u_u0_pinv * JIVP_u_z;

	MatrixXd JBVP_p_ft = JIVP_p_ft - JIVP_p_u0 * JIVP_u_u0_pinv * JIVP_u_ft;
	MatrixXd JBVP_ws_ft = JIVP_ws_ft - JIVP_ws_u0 * JIVP_u_u0_pinv * JIVP_u_ft;
	MatrixXd Jft_z = -JBVP_p_ft.completeOrthogonalDecomposition().pseudoInverse() * JBVP_p_z;  // Is this the best option in Eigen ???

	return { JBVP_p_z, JBVP_ws_z, JBVP_p_ft, JBVP_ws_ft, Jft_z };
}



template <typename adType, template<typename> typename IVPJacobians>
void CRMSolverIVP_CoreWithJacobian(CRMIVPCoreParams<adType> in_params,
	adType in_u[3], adType in_ftip[3],
	AugmentedStateVector<adType, IVPJacobians>& out_x_N, adType out_MomentResidual[3]) {

	// convenience definitions
	const double Identity3x3[9] = { 1,0,0,0,1,0,0,0,1 };
	const double Zero3[3] = { 0,0,0 };

	double l_zero[3] = { 0.0,0.0,0.0 };
	adType h;						// integration stepsize
	adType RigidSegmentLength;		// Length of the rigid segment - intermediate variable
	adType deltau[3];
	//adType muhat[9];
	//adType RscTB0[3], Tb[3], deltau1[3], K1deltau1[3], K2invResidual[3]; // intermediate variables
	int   fsegno; 					// flexible segment no
	int   actno, fsegi, fsegip1;	// actuator no, flexible segment before, flexible segment after
	bool  LastSegmentIsRigid = true;	// Flag indicating if the last segment processed is rigid (true) or not (false)

	// we need to copy in_ftip to local variable
	adType ftip[3];
	mCopy_AB<3>(in_ftip, ftip);
	// we will copy anything we will access more than once (or write to) to local variables
	int   StartSegmentIndex = in_params.StartSegmentIndex;
	int	  NextLocMarker = in_params.NextLocMarker;
	int	  InitialLocMarker = NextLocMarker;
	// for others, we will create aliases
	auto& SegBounds = in_params.SegBounds;
	auto& SegSteps = in_params.SegSteps;
	auto& InsertedLength = in_params.InsertedLength;
	auto& dlambdainv = in_params.dlambdainv;
	auto& K = in_params.K;
	auto& Kinv = in_params.Kinv;
	auto& ustar = in_params.ustar;
	auto& fcumlambda = in_params.fcumlambda;
	bool FinalValueOnly = true;   // we do not need to calculate the LocMarker positions for Jacobian calculation
	auto& LocMarkers = in_params.LocMarkers;
	auto& B0 = in_params.B0;
	auto& MagMoment = in_params.MagMoment;
	auto& CoilAlignmentTurnAreaMatrix = in_params.CoilAlignmentTurnAreaMatrix;

	AugmentedStateVector<adType, IVPJacobians> xi;  		//  Initial value of the state for the next segment to be integrated
	auto& xf = out_x_N;					//  Final value of the state for the last segment integrated
	auto& Residual = out_MomentResidual;	// Residual at the catheter tip -- will be returned
	double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3];	// we need to allocate the buffer to pass to ABM4 even though it will not be used

	// initial conditions for the regular states
	// we need to copy xi from in_params to the local variable and update it with u[0..2] specified in in_u
	mCopy_AB<3>(in_params.xi + 0, xi._p);
	mCopy_AB<9>(in_params.xi + 3, xi._R);
	mCopy_AB<3>(in_u, xi._u);

	// initial conditions for the augmented states
	//for (int i = 0; i < 3; i++)  xi._p_p0[i * 3 + i] = 1.0;			// _p_p0 = I33
	//for (int i = 0; i < 3 * 3; i++) xi._ws_w0[i] = xi._R[i];		// _ws_w0 = R0
	for (int i = 0; i < 3; i++)  xi._u_u0[i * 3 + i] = 1.0;			// _u_u0 = I33
	if constexpr (std::is_same_v<expr_type<IVPJacobians<adType>>, expr_type<IVPJacobiansFull<adType>>>) {
		for (int i = 0; i < 3; i++) xi._p_zl[i] = xi._R[i * 3 + 2];		// _p_zl = R0_.3
	}

	// IMPORTANT NOTE: most proximal segment is assumed to be always flexible
	//    and the flexible and rigid segments are assumed to be alternating
	//    most distal segment can be flexible or rigid

	// If the starting segment is a rigid segment, then we will need to move initial conditions to the start of the next flexible segment and change the StartingSegment to that segment
	if (StartSegmentIndex % 2 == 1) {
		RigidSegmentLength = SegBounds[StartSegmentIndex + 1] - SegBounds[StartSegmentIndex];	// how far we need to move along the length of the rigid segment to reach the next flexible segment
		for (int i = 0; i < 3; i++) { 										// u[0..2] and R[0..8] remain the same
			xi._p[i] = xi._p[i] + RigidSegmentLength * xi._R[i * 3 + 2];	// p[0..2] will translate along the z direction of the R matrix (3rd column)
		}
		StartSegmentIndex++;  // move the start segment to the next segment

		// initial conditions for the augmented states
		//for (int i = 0; i < 3; i++) {								// _p_w0 = [ -d * R0_2, d* R0_1, 0]  for rigid initial segment 
		//	xi._p_w0[i * 3 + 0] = - xi._R[i * 3 + 1] * RigidSegmentLength;
		//	xi._p_w0[i * 3 + 1] =   xi._R[i * 3 + 0] * RigidSegmentLength;
		//	xi._p_w0[i * 3 + 2] = 0.0;
		//}
		if constexpr (std::is_same_v<expr_type<IVPJacobians<adType>>, expr_type<IVPJacobiansFull<adType>>>) {
			for (int i = 0; i < 9; i++)  xi._ws_zl[i] = 0.0;			// _ws_zl = 0 for rigid initial segment 
			for (int i = 0; i < 3; i++)  xi._u_zl[i] = 0.0;				// _u_zl = 0 for rigid initial segment 
		}
	}
	else {
		// initial conditions for the augmented states
		if constexpr (std::is_same_v<expr_type<IVPJacobians<adType>>, expr_type<IVPJacobiansFull<adType>>>) {
			mMult_AB<3, 3, 1>(xi._R, xi._u, xi._ws_zl);					// _ws_zl = R * u
		}
		StateDerivativeVector<adType> xdot_temp;	// intermediate variable
		CRMIntegrand<adType>(SegBounds[0] /* 0.0 */, static_cast<StateVector<adType>>(xi), InsertedLength, dlambdainv, K[StartSegmentIndex], Kinv[StartSegmentIndex], l_zero, ustar[StartSegmentIndex], fcumlambda, ftip, xdot_temp);
		if constexpr (std::is_same_v<expr_type<IVPJacobians<adType>>, expr_type<IVPJacobiansFull<adType>>>) {
			for (int i = 0; i < 3; i++)  xi._u_zl[i] = xdot_temp._u[i];	// _u_zl = Fu(0)
		}
	}

	// In case there is nothing to integrate or the rigid segment is the only segment (it is both the start and the end segment),
	//   in which case, the main integration loop will not execute, 
	//   we need to have valid return values for xf and 
	//   Residual is initialized to the 0 vector and residual is initialized to the 0 vector
	xf = xi;
	for (int i = 0; i < 3; i++) Residual[i] = 0.0;


	// Integrate each of the remaining segments
	for (int i = StartSegmentIndex; i < NUM_SEGMENTS; i++) {
		if (i % 2 == 0) {  // Flexible Segment
			LastSegmentIsRigid = false;
			// Prepare the CRMIntegrand Parameters
			fsegno = i >> 1; // i/2, flexible segment no

			// Calculate the actual stepsize, based on the number of steps
			h = (SegBounds[i + 1] - SegBounds[i]) / (SegSteps[fsegno] * 1.0);

			// Integrate
			ABM4(xi, SegBounds[i], SegSteps[fsegno], h,
				InsertedLength, dlambdainv, K[fsegno], Kinv[fsegno], l_zero, ustar[fsegno], fcumlambda, ftip,
				FinalValueOnly, LocMarkers, NextLocMarker,
				xf, p_atLocMarkers);

			// residual = K (u1 - u1star)
			mSub_AB<3, 1>(xf._u, ustar[fsegno], deltau);
			mMult_AB<3, 3, 1>(K[fsegno], deltau, Residual);

		}
		else {  // Need to do actuation/rigid segment calculations to transfer Initial Conditions to next flexible segment
			LastSegmentIsRigid = true;
			RigidSegmentLength = (SegBounds[i + 1] - SegBounds[i]);
			actno = (i - 1) >> 1;		// actuator no
			fsegi = actno;				// index of flexible segment before (immediately proximal to) the rigid link
			fsegip1 = actno + 1;		// index of flexible segment after (immediately distal to) the rigid link
			if (i == (NUM_SEGMENTS - 1)) {	// if we are at the last segment, ustar[fsegip1] and Kinv[fsegip1] are assigned to zero and identity matrix, respectively
				CRMSolverIVP_PropagateBCThroughRigidLink(xf, Residual, RigidSegmentLength, actno, MagMoment[actno], CoilAlignmentTurnAreaMatrix[actno], B0, ustar[fsegi], K[fsegi], Zero3, Identity3x3, xi, Residual);
			}
			else
				CRMSolverIVP_PropagateBCThroughRigidLink(xf, Residual, RigidSegmentLength, actno, MagMoment[actno], CoilAlignmentTurnAreaMatrix[actno], B0, ustar[fsegi], K[fsegi], ustar[fsegip1], Kinv[fsegip1], xi, Residual);

		}
		xi = xf; // The calculated final values will be the initial value of the next iteration
	}

	// no need to copy the final values of the state x_f and residual to the output
	// since we have created an alias
	//out_x_N = xf;
	//mCopy_AB<3>(Residual, out_MomentResidual);

}


// each column k of  Pnew_k = P_k + (w_k^ R)_3 * d, where (w_k^ * R)_3 is the third column of the (w_k^ * R) matrix ( 3x3 matrices packed in row major order )
template<unsigned int dim, typename adType>
inline void PpluswhatR3timesD(const adType P[], const adType w[], const adType R[], const adType d, adType Pnew[]) {
	for (int k = 0; k < dim; k++) {  // loop over columns of P
		Pnew[0 * dim + k] = P[0 * dim + k] + (-w[2 * dim + k] * R[1 * 3 + 2] + w[1 * dim + k] * R[2 * 3 + 2]) * d;
		Pnew[1 * dim + k] = P[1 * dim + k] + ( w[2 * dim + k] * R[0 * 3 + 2] - w[0 * dim + k] * R[2 * 3 + 2]) * d;
		Pnew[2 * dim + k] = P[2 * dim + k] + (-w[1 * dim + k] * R[0 * 3 + 2] + w[0 * dim + k] * R[1 * 3 + 2]) * d;
	}
}

// we will overload this function
template<typename adType, template<typename> typename IVPJacobians>
void CRMSolverIVP_PropagateBCThroughRigidLink(AugmentedStateVector<adType, IVPJacobians>& xi_ip1, adType Residual_ip1[3],
	const adType RigidSegmentLength, const unsigned int ActNo, const adType MagMoment[3], const double CoilAlignmentTurnAreaMatrix[9],
	const double B0[3], const double ustar_i[3], const double K_i[9], const double ustar_ip1[3], const double Kinv_ip1[9],
	const AugmentedStateVector<adType, IVPJacobians>& xf_i, const adType Residual_i[3]) {

	// for debugging
	//std::cout << "~~~~~~ xf_i =    ~~~~~~ \n" << xf_i << std::endl;

	// Let's first propagate the regular state boundary conditions
	CRMSolverIVP_PropagateBCThroughRigidLink(xi_ip1, Residual_ip1, RigidSegmentLength, ActNo, MagMoment, CoilAlignmentTurnAreaMatrix, B0, ustar_i, K_i, ustar_ip1, Kinv_ip1, static_cast<StateVector<adType>>(xf_i), Residual_i);

	adType muhat[9]; // intermediate variables
	wHat<adType>(MagMoment, muhat);

	// for simplicity, define an alias
	constexpr unsigned int Cs = CURRENT_ACT_VECTOR_DIM;

	// p
	// dpdPhi_ip1 = dpdPhi_i + d * (ws_Phi_i^ * R)_3
	//PpluswhatR3timesD<3>(xf_i._p_p0, xf_i._ws_p0, xf_i._R, RigidSegmentLength, xi_ip1._p_p0);
	//PpluswhatR3timesD<3>(xf_i._p_w0, xf_i._ws_w0, xf_i._R, RigidSegmentLength, xi_ip1._p_w0);
	PpluswhatR3timesD<3>(xf_i._p_u0, xf_i._ws_u0, xf_i._R, RigidSegmentLength, xi_ip1._p_u0);
	if constexpr (std::is_same_v<expr_type<IVPJacobians<adType>>, expr_type<IVPJacobiansFull<adType>>>) {
		PpluswhatR3timesD<Cs>(xf_i._p_zc, xf_i._ws_zc, xf_i._R, RigidSegmentLength, xi_ip1._p_zc);
		PpluswhatR3timesD<1>(xf_i._p_zl, xf_i._ws_zl, xf_i._R, RigidSegmentLength, xi_ip1._p_zl);
		PpluswhatR3timesD<3>(xf_i._p_ft, xf_i._ws_ft, xf_i._R, RigidSegmentLength, xi_ip1._p_ft);
	}

	// R
	// dwsdPhi_ip1 = dwsdPhi_i
	//mCopy_AB<3 * 3>(xf_i._ws_p0, xi_ip1._ws_p0);
	//mCopy_AB<3 * 3>(xf_i._ws_w0, xi_ip1._ws_w0);
	mCopy_AB<3 * 3>(xf_i._ws_u0, xi_ip1._ws_u0);
	if constexpr (std::is_same_v<expr_type<IVPJacobians<adType>>, expr_type<IVPJacobiansFull<adType>>>) {
		mCopy_AB<3 * Cs>(xf_i._ws_zc, xi_ip1._ws_zc);
		mCopy_AB<3 * 1>(xf_i._ws_zl, xi_ip1._ws_zl);
		mCopy_AB<3 * 3>(xf_i._ws_ft, xi_ip1._ws_ft);
	}

	// u
	// dudPhi_ip1 = Kinv_ip1 * K_i * dudPhi_i - Kinv_ip1 dTaudPhi
	// dTaudPhi = dTaudws * dwsPhi  for all Phi except zc_link  (zc_link : actuation corresponding to current link)
	// dTaudzc_link = dTaudzc_link + dTaudws * dwsdzc_link  for zc_link
	// - Kinv_ip1 dTaudws * dwsPhi = Kinv_ip1 * mu^ * R^T * ws_Phi^ * B0


	// Kinv_ip1K_i = Kinv_ip1 * K_i
	adType Kinv_ip1K_i[9];
	mMult_AB<3, 3, 3>(Kinv_ip1, K_i, Kinv_ip1K_i);

	// common terms
	// Kinv_ip1* K_i* dudPhi_i
	//mMult_AB<3, 3, 3>(Kinv_ip1K_i, xf_i._u_p0, xi_ip1._u_p0);
	//mMult_AB<3, 3, 3>(Kinv_ip1K_i, xf_i._u_w0, xi_ip1._u_w0);
	mMult_AB<3, 3, 3>(Kinv_ip1K_i, xf_i._u_u0, xi_ip1._u_u0);
	if constexpr (std::is_same_v<expr_type<IVPJacobians<adType>>, expr_type<IVPJacobiansFull<adType>>>) {
		mMult_AB<3, 3, Cs>(Kinv_ip1K_i, xf_i._u_zc, xi_ip1._u_zc);
		mMult_AB<3, 3, 1>(Kinv_ip1K_i, xf_i._u_zl, xi_ip1._u_zl);
		mMult_AB<3, 3, 3>(Kinv_ip1K_i, xf_i._u_ft, xi_ip1._u_ft);
	}

	// Kinv_ip1muhatRT = Kinv_ip1 * mu^ * R^T
	adType muhatRT[9], Kinv_ip1muhatRT[9];
	mMult_ABT<3, 3, 3>(muhat, xf_i._R, muhatRT);
	mMult_AB<3, 3, 3>(Kinv_ip1, muhatRT, Kinv_ip1muhatRT);
	// define lambda function to calculate xi_ip1._u_Phi += Kinv_ip1muhatRT * ws_Phi^ * B0, for a given Phi
	auto addmKinvdTaudPhi = [&B0, &Kinv_ip1muhatRT](adType* _ws_Phi, adType* _u_Phi, unsigned int stride) {
		for (unsigned int j = 0; j < stride; j++) {
			adType what[9], whatB0[3], Kinv_ip1muhatRTwhatB0[3];
			wHat(_ws_Phi + j, what, stride);
			mMult_AB<3, 3, 1>(what, B0, whatB0);
			mMult_AB<3, 3, 1>(Kinv_ip1muhatRT, whatB0, Kinv_ip1muhatRTwhatB0);
			for (int i = 0; i < 3; i++) _u_Phi[i * stride + j] += Kinv_ip1muhatRTwhatB0[i];
		}
	};

	// xi_ip1._u_Phi += Kinv_ip1muhatRT * ws_Phi^ * B0, for every Phi
	//addmKinvdTaudPhi(xf_i._ws_p0, xi_ip1._u_p0, 3);
	//addmKinvdTaudPhi(xf_i._ws_w0, xi_ip1._u_w0, 3);
	addmKinvdTaudPhi(xf_i._ws_u0, xi_ip1._u_u0, 3);
	if constexpr (std::is_same_v<expr_type<IVPJacobians<adType>>, expr_type<IVPJacobiansFull<adType>>>) {
		addmKinvdTaudPhi(xf_i._ws_zc, xi_ip1._u_zc, Cs);
		addmKinvdTaudPhi(xf_i._ws_zl, xi_ip1._u_zl, 1);
		addmKinvdTaudPhi(xf_i._ws_ft, xi_ip1._u_ft, 3);

		// dTaudPhi term for the specific actuator that is on the current link
		// - Kinv_ip1 * dTaudzc_link
		//
		// mKinv_ip1dTaudzc_link = Kinv_ip1 * mdTaudzc_link = Kinv_ip1 * -(- (R^T B0)^ (CATA)) , where CATA = CoilAlignmentMatrix * CoilTurnAreaMatrix
		adType RTB0[3];
		mMult_ATB<3, 3, 1>(xf_i._R, B0, RTB0);
		adType RTB0hat[9];
		wHat(RTB0, RTB0hat);
		adType mdTaudzc_link[3 * 3];
		mMult_AB<3, 3, 3>(RTB0hat, CoilAlignmentTurnAreaMatrix, mdTaudzc_link);
		adType mKinv_ip1dTaudzc_link[3 * 3];
		mMult_AB<3, 3, 3>(Kinv_ip1, mdTaudzc_link, mKinv_ip1dTaudzc_link);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				xi_ip1._u_zc[i * Cs + (ActNo * 3 + j)] += mKinv_ip1dTaudzc_link[i * 3 + j];
	}

	// for debugging
	//std::cout << "~~~~~~  xi_ip1 =   ~~~~~~\n" << xi_ip1 << std::endl;
}



template <typename adType, template<typename> typename IVPJacobians>
void CRMIntegrand (	const adType s, const AugmentedStateVector<adType, IVPJacobians>& in_x, const adType in_Li, const double in_dlambdainv,
					const double in_K[9], const double in_Kinv[9], const double in_l[3], const double in_ustar[3], const double in_fcumlambda[NUM_FCUM_LAMBDA+1][3],
					const adType in_ftip[3],
					AugmentedStateDerivativeVector<adType, IVPJacobians>& out_xdot) {

	using namespace Eigen;

	// for simplicity, define an alias
	constexpr unsigned int Cs = CURRENT_ACT_VECTOR_DIM;

	// first evaluate partial derivative of x w.r.t. s using CRMIntegrand
	CRMIntegrand(s, static_cast<StateVector<adType>>(in_x), in_Li, in_dlambdainv, in_K, in_Kinv, in_l, in_ustar, in_fcumlambda, in_ftip, out_xdot);

	// for simplicity, create aliases
	auto& Length = in_Li;
	auto& deltalambdainv = in_dlambdainv;
	auto& K = in_K;
	auto& Kinv = in_Kinv;
	//__EMT<adType, 3, 3> Kinv(in_Kinv);
	auto& ustar = in_ustar;
	auto& l = in_l;
	//  assume ustardot = 0
	auto& x = in_x;
	auto& p = in_x._p;
	auto& R = in_x._R;
	auto& u = in_x._u;

	adType fcum[3];
	// calculate interpolated value of fcum
	adType lambda = Length - s;
	adType ix = lambda * deltalambdainv;
	double ird_f = floor(dVal(ix)); // index for round down  -- doubleing point
	if (ird_f < 0) ird_f = 0;
	int ird = (int)ird_f;	//    integer index
	double iru_f = ceil(dVal(ix)); 	// index for round up  -- doubleing point
	if (iru_f > NUM_FCUM_LAMBDA) iru_f = NUM_FCUM_LAMBDA;
	int iru = (int)iru_f;	//    integer index
	adType ixmird = ix - ird;    // weight for interpolation
	adType irumix = iru - ix;	// weight for interpolation
	for (int i = 0; i < 3; i++) {
		fcum[i] = in_fcumlambda[iru][i] * ixmird + in_fcumlambda[ird][i] * irumix;
	}
	// add the tip force to fcum
	for (int i = 0; i < 3; i++) {
		fcum[i] += in_ftip[i];
	}


	//
	// now evaluate partial derivative of ( \partial x / \partial \phi ) x w.r.t. s 
	//

	// intermediate terms needed
	adType uhat[3 * 3];
	wHat(u, uhat);

	//
	// _p terms
	//

	// define lambda function for calculating derivative of _p_Phi
	//   _p_phi = \hat{ws_phi} * R_3
	auto dot_p_Phi = [&R](adType* x_ws_Phi, adType* xdot_p_Phi, unsigned int stride) {
		for (unsigned int j = 0; j < stride; j++) {
			xdot_p_Phi[0 * stride + j] = -x_ws_Phi[2 * stride + j] * R[1 * 3 + 2] + x_ws_Phi[1 * stride + j] * R[2 * 3 + 2];
			xdot_p_Phi[1 * stride + j] = x_ws_Phi[2 * stride + j] * R[0 * 3 + 2] - x_ws_Phi[0 * stride + j] * R[2 * 3 + 2];
			xdot_p_Phi[2 * stride + j] = -x_ws_Phi[1 * stride + j] * R[0 * 3 + 2] + x_ws_Phi[0 * stride + j] * R[1 * 3 + 2];
		}
	};

	//_p_p0
	//dot_p_Phi(x._ws_p0, out_xdot._p_p0,3);

	//_p_w0
	//dot_p_Phi(x._ws_w0, out_xdot._p_w0, 3);

	//_p_u0
	dot_p_Phi(x._ws_u0, out_xdot._p_u0, 3);

	if constexpr (std::is_same_v<expr_type<IVPJacobians<adType>>, expr_type<IVPJacobiansFull<adType>>>) {
		//_p_zc
		dot_p_Phi(x._ws_zc, out_xdot._p_zc, Cs);

		//_p_zl
		dot_p_Phi(x._ws_zl, out_xdot._p_zl, 1);

		//_p_ft
		dot_p_Phi(x._ws_ft, out_xdot._p_ft, 3);
	}

	//
	// _ws terms
	//

	//   _ws_Phi = [ dFR_du * _u_Phi * R^T ]V;   where V is the 'vee' operator
	//           =  R _u_Phi

	//_ws_p0
	//mMult_AB<3, 3, 3>(R, x._u_p0, out_xdot._ws_p0);

	//_ws_w0
	//mMult_AB<3, 3, 3>(R, x._u_w0, out_xdot._ws_w0);

	//_ws_u0
	mMult_AB<3, 3, 3>(R, x._u_u0, out_xdot._ws_u0);

	if constexpr (std::is_same_v<expr_type<IVPJacobians<adType>>, expr_type<IVPJacobiansFull<adType>>>) {
		//_ws_zc
		mMult_AB<3, 3, Cs>(R, x._u_zc, out_xdot._ws_zc);

		//_ws_zl
		mMult_AB<3, 3, 1>(R, x._u_zl, out_xdot._ws_zl);

		//_ws_ft
		mMult_AB<3, 3, 3>(R, x._u_ft, out_xdot._ws_ft);
	}

	//
	// _u terms
	//

	// common terms
	adType Kinve3hat[3 * 3] = { Kinv[0 * 3 + 1], -Kinv[0 * 3 + 0], 0, /*;*/ Kinv[1 * 3 + 1], -Kinv[1 * 3 + 0], 0, /*;*/ Kinv[2 * 3 + 1], -Kinv[2 * 3 + 0], 0 };
	adType Kinve3hatRT[3 * 3], KinvRT[3 * 3];
	mMult_ABT<3, 3, 3>(Kinve3hat, R, Kinve3hatRT);
	mMult_ABT<3, 3, 3>(Kinv, R, KinvRT);

	// define lambda function to calculate dFu_dws * _ws_Phi term of dFu_dPhi
	//   
	auto dFu_dws_dws_dPhi = [&fcum, &l, &Kinve3hatRT, &KinvRT](adType* _ws_Phi, adType* xdot_u_Phi, unsigned int stride) {
		adType _ws_Phihat[9], _ws_Phihatfcum[3], _ws_Phihatl[3], temp[3];
		for (unsigned int j = 0; j < stride; j++) {
			wHat(_ws_Phi + j, _ws_Phihat, stride);
			mMult_AB<3, 3, 1>(_ws_Phihat, fcum, _ws_Phihatfcum);
			mMult_AB<3, 3, 1>(_ws_Phihat, l, _ws_Phihatl);
			mMult_AB<3, 3, 1>(Kinve3hatRT, _ws_Phihatfcum, temp);
			mMultAdd_AB<3, 3, 1>(KinvRT, _ws_Phihatl, temp);
			for (int i = 0; i < 3; i++)	xdot_u_Phi[i * stride + j] = temp[i];
		}
	};

	//Matrix<adType, 3, 3, RowMajor> dFu_du = ;
	// FOR THIS CALCULATION, WE ARE ASSUMING THAT: K and Kinv matrices are diagonal
	auto& K0 = K[0 * 3 + 0];
	auto& K1 = K[1 * 3 + 1];
	auto& K2 = K[2 * 3 + 2];
	auto& K0inv = Kinv[0 * 3 + 0];
	auto& K1inv = Kinv[1 * 3 + 1];
	auto& K2inv = Kinv[2 * 3 + 2];
	auto& u0 = u[0];
	auto& u1 = u[1];
	auto& u2 = u[2];
	auto& u0st = ustar[0];
	auto& u1st = ustar[1];
	auto& u2st = ustar[2];
	adType dFu_du[3 * 3] = { 0.0,                                          K0inv * K1 * u2 - K0inv * K2 * (u2 - u2st),   -K0inv * K2 * u1 + K0inv * K1 * (u1 - u1st), /*;*/
		                     -K1inv * K0 * u2 + K1inv * K2 * (u2 - u2st),  0.0,                                          K1inv * K2 * u0 - K1inv * K0 * (u0 - u0st), /*;*/
		                     K2inv * K0 * u1 - K2inv * K1 * (u1 - u1st),   -K2inv * K1 * u0 + K2inv * K0 * (u0 - u0st),  0.0                                          };

	
	//_u_p0
	//out_xdot._u_p0 = dFu_dws * _ws_p0 + dFu_du * _u_p0; 
	//dFu_dws_dws_dPhi(x._ws_p0, out_xdot._u_p0, 3);
	//mMultAdd_AB<3, 3, 3>(dFu_du, x._u_p0, out_xdot._u_p0);

	//_u_w0
	//out_xdot._u_w0 = dFu_dws * _ws_w0 + dFu_du * _u_w0; 
	//dFu_dws_dws_dPhi(x._ws_w0, out_xdot._u_w0, 3);
	//mMultAdd_AB<3, 3, 3>(dFu_du, x._u_w0, out_xdot._u_w0);

	//_u_u0
	//out_xdot._u_u0 = dFu_dws * _ws_u0 + dFu_du * _u_u0; 
	dFu_dws_dws_dPhi(x._ws_u0, out_xdot._u_u0, 3);
	mMultAdd_AB<3, 3, 3>(dFu_du, x._u_u0, out_xdot._u_u0);

	if constexpr (std::is_same_v<expr_type<IVPJacobians<adType>>, expr_type<IVPJacobiansFull<adType>>>) {
		//_u_zc
		//out_xdot._u_zc = dFu_dws * _ws_zc + dFu_du * _u_zc; 
		dFu_dws_dws_dPhi(x._ws_zc, out_xdot._u_zc, Cs);
		mMultAdd_AB<3, 3, Cs>(dFu_du, x._u_zc, out_xdot._u_zc);

		//_u_zl
		//out_xdot._u_zl = dFu_dws * _ws_zl + dFu_du * _u_zl; 
		dFu_dws_dws_dPhi(x._ws_zl, out_xdot._u_zl, 1);
		mMultAdd_AB<3, 3, 1>(dFu_du, x._u_zl, out_xdot._u_zl);

		//  dFu_dft = - Kinv * e3hat * R'
		//   (- e3hat * R') = [ r12 r22 r32; -r11 -r21 -r31; 0 0 0];
		//Matrix<adType, 3, 3, RowMajor> me3hatRT{ { R[1], R[4], R[7] }, { -R[0], -R[3], -R[6] }, { 0.0, 0.0, 0.0 } };
		//Matrix<adType, 3, 3, RowMajor> dFu_dft = Kinv * me3hatRT;
		adType me3hatRT[3 * 3] = { R[1], R[4], R[7], /*;*/ -R[0], -R[3], -R[6], /*;*/ 0.0, 0.0, 0.0 };
		adType dFu_dft[3 * 3];
		mMult_AB<3, 3, 3>(Kinv, me3hatRT, dFu_dft);

		//_u_ft
		//out_xdot._u_ft = dFu_dws * _ws_ft + dFu_du * _u_ft + dFu_dft; 
		dFu_dws_dws_dPhi(x._ws_ft, out_xdot._u_ft, 3);
		mMultAdd_AB<3, 3, 3>(dFu_du, x._u_ft, out_xdot._u_ft);
		mAdd_AB<3, 3>(out_xdot._u_ft, dFu_dft);
	}

}


