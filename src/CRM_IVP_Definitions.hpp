#pragma once
#include <cmath>


template <typename adType>
void CRMSolverIVP (	double in_x_0[NUM_STATES], double in_IntegrationStepSize,
					double in_Li, double in_dlambdainv,
					double in_SegEndLambdas[NUM_SEGMENTS], double	in_LocMarkerLambdas[NUM_LOCALIZATION_MARKERS], 
					double in_K[NUM_FLEX_SEG][9], double in_Kinv[NUM_FLEX_SEG][9], double in_ustar[NUM_FLEX_SEG][3], 
					double in_CoilAlignmentTurnAreaMatrix[NUM_ACT_SET][9],
					adType in_MagMoment[NUM_ACT_SET][3], double in_fcumlambda[NUM_FCUM_LAMBDA+1][3], adType in_ftip[3],
					double in_B0[3],
					bool in_FinalValueOnly,
					double out_x_N[NUM_STATES], double out_MomentResidual[3],
					double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]
					) {


	CRMIVPCoreParams<adType> CoreParams;
	adType x_0[NUM_STATES];
	mCopy_AB<NUM_STATES>(in_x_0, x_0);
	StateVector<adType> x_N;
	adType u_0[3], ftip[3];
	adType MomentResidual[3];
	double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3];

	// copy to local variable
	adType MagMoment[NUM_ACT_SET][3];
	mCopy_ABm<NUM_ACT_SET,3>(in_MagMoment, MagMoment);
	adType Li = in_Li;

	CRMSolverIVP_Prep(x_0, in_IntegrationStepSize,
		Li, in_dlambdainv,
		in_SegEndLambdas, in_LocMarkerLambdas,
		in_K, in_Kinv, in_ustar,
		in_CoilAlignmentTurnAreaMatrix,
		MagMoment, in_fcumlambda,
		in_B0,
		in_FinalValueOnly,
		CoreParams);
	// copy to local variable
	for (int i=0; i<3; i++) ftip[i]=in_ftip[i];

	// We need to pass u_0 as input argument as the values in CoreParams will be overriden with the values provided in the input arguments - functionality needed for solving Boundary Value Problems (BVP)
	for (int i = 0; i < 3; i++) {
		u_0[i] = in_x_0[3 + 9 + i];
	}
	CRMSolverIVP_Core ( CoreParams, u_0, ftip, x_N, MomentResidual, p_atLocMarkers );

	CRMSolverIVP_Return( x_N, MomentResidual, p_atLocMarkers, out_x_N, out_MomentResidual, out_p_atLocMarkers);

}


template <typename adType>
void CRMSolverIVP(	CRMShootingMethodParams<adType> in_Params,
					adType in_u0[3], adType in_ftip[3],
					bool in_FinalValueOnly,
					adType out_x_N[NUM_STATES], adType out_MomentResidual[3], 
					double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]		) {

	adType x_0[NUM_STATES];
	for (int i = 0; i < 3; i++) x_0[i] = in_Params.p0[i];
	for (int i = 0; i < 9; i++) x_0[i + 3] = in_Params.R0[i];
	for (int i = 0; i < 3; i++) x_0[i + 3 + 9] = in_u0[i];

	CRMIVPCoreParams<adType> CoreParams;
	adType u_0[3], ftip[3];

	StateVector<adType> x_N;

	CRMSolverIVP_Prep(x_0, in_Params.IntegrationStepSize,
		in_Params.Li, in_Params.dlambdainv,
		in_Params.SegEndLambdas, in_Params.LocMarkerLambdas,
		in_Params.K, in_Params.Kinv, in_Params.ustar,
		in_Params.CoilAlignmentTurnAreaMatrix,
		in_Params.MagMoment, in_Params.fcumlambda,
		in_Params.B0,
		in_FinalValueOnly,
		CoreParams);

	// copy to local variable
	for (int i = 0; i < 3; i++) ftip[i] = in_ftip[i];

	// We need to pass u_0 as input argument as the values in CoreParams will be overriden with the values provided in the input arguments - functionality needed for solving Boundary Value Problems (BVP)
	for (int i = 0; i < 3; i++) u_0[i] = in_u0[i];

	CRMSolverIVP_Core(CoreParams, u_0, ftip, x_N, out_MomentResidual, out_p_atLocMarkers);

	// NO NEED FOR CRMSolverIVP_Return
	mCopy_AB<3>(x_N._p, out_x_N + 0);
	mCopy_AB<9>(x_N._R, out_x_N + 3);
	mCopy_AB<3>(x_N._u, out_x_N + 3 + 9);

}


template <typename adType>
void CRMSolverIVP_Prep ( adType in_x_0[NUM_STATES], double in_IntegrationStepSize,
						 adType in_Li, double in_dlambdainv,
						 double in_SegEndLambdas[NUM_SEGMENTS], double in_LocMarkerLambdas[NUM_LOCALIZATION_MARKERS],
						 double in_K[NUM_FLEX_SEG][9], double in_Kinv[NUM_FLEX_SEG][9], double in_ustar[NUM_FLEX_SEG][3],
						 double in_CoilAlignmentTurnAreaMatrix[NUM_ACT_SET][9],
						 adType in_MagMoment[NUM_ACT_SET][3], double in_fcumlambda[NUM_FCUM_LAMBDA+1][3],
						 double in_B0[3],
						 bool in_FinalValueOnly,
						 CRMIVPCoreParams<adType> &out_CoreParams	) {
	// Process the incoming parameters (including changing from distal-proximal order to proximal-distal order)
	//   and package them to be passed to CRMSolverIVP_Core
	// Everything is copied to local variables (inside out_CoreParams), therefore, this function
	//   can be used to transfer data from global memory to device memory for subsequent computations
	// The packaged data can be used to call CRMSolverIVP_Core multiple times by changing
	//   only u[0..2] components of xi --- other parameters should not change

	// local copies for variables accessed out-of-order
	double SegEndLambdas[NUM_SEGMENTS];
	mCopy_AB<NUM_SEGMENTS>(in_SegEndLambdas,SegEndLambdas);

	// create aliases for variables in CoreParams
	auto & xi = out_CoreParams.xi;

	auto & SegBounds = out_CoreParams.SegBounds;
	auto & SegSteps = out_CoreParams.SegSteps;

	auto & InsertedLength = out_CoreParams.InsertedLength;
	auto & dlambdainv = out_CoreParams.dlambdainv;
	auto & K = out_CoreParams.K;
	auto & Kinv = out_CoreParams.Kinv;
	auto & ustar = out_CoreParams.ustar;
	auto & CoilAlignmentTurnAreaMatrix = out_CoreParams.CoilAlignmentTurnAreaMatrix;
	auto & fcumlambda = out_CoreParams.fcumlambda;
//	auto & ftip = out_CoreParams.ftip;
	auto & FinalValueOnly = out_CoreParams.FinalValueOnly;
	auto & LocMarkers = out_CoreParams.LocMarkers;

	auto & B0 = out_CoreParams.B0;
	auto & MagMoment = out_CoreParams.MagMoment;

	auto & StartSegmentIndex = out_CoreParams.StartSegmentIndex;
	auto & NextLocMarker = out_CoreParams.NextLocMarker;
	auto & p_atLocMarkers = out_CoreParams.p_atLocMarkers;

	// local variables
	double IntegrationStepSize=in_IntegrationStepSize;
	double DeltaSInv=1.0/IntegrationStepSize;
	adType tempadType;

	// process parameters as needed and copy into CoreParams
	mCopy_AB<NUM_STATES>(in_x_0,xi);		// Initial value of the state for the next segment to be integrated
											//  States are packed p[0..2],R[0..8],u[0..2]  (R: 3x3 matrix stored in row major order R11 R12 R13 R21 R22 R23 R31 R32 R33)
	mCopy_AB<3>(in_B0,B0);					//  B0 field vector of the MRI scanner (in spatial coordinates)

	dlambdainv=in_dlambdainv;				// reciprocal of dlambda (lambda stepsize used in discretizing fcumlambda)
	FinalValueOnly=in_FinalValueOnly;		// Flag used to indicate if only final value (xf) is returned (true) or if Marker Locations are returned as well (false)
	InsertedLength=in_Li;					// Inserted Length (length of the catheter from the entry point to the tip)
	// If the catheter is inserted more than the length of the catheter, clamp it to catheter length
	if (InsertedLength > SegEndLambdas[NUM_SEGMENTS - 1]) InsertedLength = SegEndLambdas[NUM_SEGMENTS - 1];

	// WE ARE GOING TO REORDER SEGMENT AND ACTUATOR UNITS SO THAT THEY ARE ORDERED FROM THE INSERTION POINT TO THE TIP
	//   I.E., SWITCH TO PROXIMAL TO DISTAL ORDERING

	StartSegmentIndex=NUM_SEGMENTS-1;		// Index of the segment where the integration to solve IVP will start -- the segment located at the entry point; note that segment indices start at 0
	SegBounds[NUM_SEGMENTS]=InsertedLength;	// End s value of last segment is s=InsertedLength
											// Entry point has a value of s=0 -- we may not simulate full length of the most proximal segment in the chamber
	for (int i=0; i<NUM_SEGMENTS; i++) {
		tempadType =InsertedLength-SegEndLambdas[i];
		if (tempadType >0.0) {
			SegBounds[(NUM_SEGMENTS-1)-i]= tempadType;
			StartSegmentIndex--;				// Integration will start at the previous segment
		}
		else {
			SegBounds[(NUM_SEGMENTS-1)-i]=0.0;	// This segment boundary is still inside the sheath
		}
	}

	for (int i=0; i<NUM_FLEX_SEG; i++) { 		// flexible catheter segment
		// Calculate the number of integration steps based on the given IntegrationStepSize
		SegSteps[i]=int(ceil( dVal(SegBounds[2*i+1]-SegBounds[2*i])*DeltaSInv ));
	}

	NextLocMarker=NUM_LOCALIZATION_MARKERS;
	if (!FinalValueOnly) {
		// While changing the order and converting from lambda to s,
		//   also find the index of the first localization marker after the entry point
		//   and assign (extrapolated) locations to markers which are still inside the sheath
		for (int i=0; i<NUM_LOCALIZATION_MARKERS; i++) {
			tempadType =InsertedLength-in_LocMarkerLambdas[i];
			LocMarkers[(NUM_LOCALIZATION_MARKERS-1)-i]= dVal(tempadType);
			if (tempadType >0.0) NextLocMarker--;
			else {   // and assign (extrapolated) locations to markers which are still inside the sheath
				LocMarkerUpdate(p_atLocMarkers[(NUM_LOCALIZATION_MARKERS - 1) - i], xi, tempadType);
			}
		}
	}

	// be careful - order is reversed in K, Kinv, ustar, CoilAlignmentTurnAreaMatrix, and MagMoment
	for (int i=0; i<NUM_FLEX_SEG; i++) {
		for (int j=0; j<9; j++) {
			K[(NUM_FLEX_SEG-1)-i][j]		=	in_K[i][j];
			Kinv[(NUM_FLEX_SEG-1)-i][j]		=	in_Kinv[i][j];
		}
		for (int j=0; j<3; j++) {
			ustar[(NUM_FLEX_SEG-1)-i][j]	=	in_ustar[i][j];
		}
	}
	for (int i=0; i<NUM_ACT_SET; i++) {
		for (int j=0; j<3; j++) {
			MagMoment[(NUM_ACT_SET-1)-i][j]	=	in_MagMoment[i][j];
		}
		mCopy_AB<9>(in_CoilAlignmentTurnAreaMatrix[i], CoilAlignmentTurnAreaMatrix[(NUM_ACT_SET - 1) - i]);
	}

	mCopy_ABm<(NUM_FCUM_LAMBDA+1),3>(in_fcumlambda,fcumlambda);
	//mCopy_AB<3>(in_ftip, ftip);

}


template <typename adType>
void CRMSolverIVP_Core ( const CRMIVPCoreParams<adType>& in_params,
						 adType in_u[3], adType in_ftip[3],
						 StateVector<adType>& out_x_N, adType out_MomentResidual[3],
						 double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3] ){

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
	auto& FinalValueOnly = in_params.FinalValueOnly;
	auto& LocMarkers = in_params.LocMarkers;
	auto& B0 = in_params.B0;
	auto& MagMoment = in_params.MagMoment;
	auto& CoilAlignmentTurnAreaMatrix = in_params.CoilAlignmentTurnAreaMatrix;

	StateVector<adType> xi;  		//  Initial value of the state for the next segment to be integrated
	auto& xf = out_x_N;					//  Final value of the state for the last segment integrated
	auto& Residual = out_MomentResidual;	// Residual at the catheter tip -- will be returned
	double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3];	// Positions at markers (ordered proximal to distal)

	// we will also copy the already filled entries of p_atLocMarkers from params
	for (int i = 0; i < InitialLocMarker; i++) {
		for (int j = 0; j < 3; j++) {
			p_atLocMarkers[i][j] = in_params.p_atLocMarkers[i][j];
		}
	}

	// initial conditions for the regular states
	// we need to copy xi from in_params to the local variable and update it with u[0..2] specified in in_u
	mCopy_AB<3>(in_params.xi + 0, xi._p);
	mCopy_AB<9>(in_params.xi + 3, xi._R);
	mCopy_AB<3>(in_u, xi._u);


	// IMPORTANT NOTE: most proximal segment is assumed to be always flexible
	//    and the flexible and rigid segments are assumed to be alternating
	//    most distal segment can be flexible or rigid

	// If the starting segment is a rigid segment, then we will need to move initial conditions to the start of the next flexible segment and change the StartingSegment to that segment
	if ( StartSegmentIndex % 2 == 1 ) {
		RigidSegmentLength	=	SegBounds[StartSegmentIndex+1]-SegBounds[StartSegmentIndex];	// how far we need to move along the length of the rigid segment to reach the next flexible segment
		for (int i = 0; i < 3; i++) { 										// u[0..2] and R[0..8] remain the same
			xi._p[i] = xi._p[i] + RigidSegmentLength * xi._R[i * 3 + 2];	// p[0..2] will translate along the z direction of the R matrix (3rd column)
		}
		StartSegmentIndex++;  // move the start segment to the next segment
		if (!FinalValueOnly) {
			// are there any localization markers?  If so, calculate their positions
			CalculateLocMarkers(NextLocMarker, xf, LocMarkers, SegBounds[StartSegmentIndex], p_atLocMarkers);  //	xf has already been updated
		}
	}

	// In case there is nothing to integrate or the rigid segment is the only segment (it is both the start and the end segment),
	//   in which case, the main integration loop will not execute, 
	//   we need to have valid return values for xf and 
	//   Residual is initialized to the 0 vector and residual is initialized to the 0 vector
	xf = xi;
	for (int i = 0; i < 3; i++) Residual[i] = 0.0;


	// Integrate each of the remaining segments
	for (int i=StartSegmentIndex; i<NUM_SEGMENTS; i++){
		if ( i%2 == 0 ) {  // Flexible Segment
			LastSegmentIsRigid=false;
			// Prepare the CRMIntegrand Parameters
			fsegno=i>>1; // i/2, flexible segment no

			// Calculate the actual stepsize, based on the number of steps
			h=(SegBounds[i+1]-SegBounds[i])/(SegSteps[fsegno]*1.0);

			// Integrate
			ABM4( 	xi, SegBounds[i], SegSteps[fsegno], h,
					InsertedLength, dlambdainv, K[fsegno], Kinv[fsegno], l_zero, ustar[fsegno], fcumlambda, ftip,
					FinalValueOnly, LocMarkers, NextLocMarker,
					xf, p_atLocMarkers	);

			// residual = K (u1 - u1star)
			mSub_AB<3, 1>(xf._u, ustar[fsegno], deltau);
			mMult_AB<3, 3, 1>(K[fsegno], deltau, Residual);

		}
		else {  // Need to do actuation/rigid segment calculations to transfer Initial Conditions to next flexible segment
			LastSegmentIsRigid=true;
			RigidSegmentLength=(SegBounds[i+1]-SegBounds[i]);
			actno = (i - 1) >> 1;		// actuator no
			fsegi = actno;				// index of flexible segment before (immediately proximal to) the rigid link
			fsegip1 = actno + 1;		// index of flexible segment after (immediately distal to) the rigid link
			if (i == (NUM_SEGMENTS - 1)) {	// if we are at the last segment, ustar[fsegip1] and Kinv[fsegip1] are assigned to zero and identity matrix, respectively
				CRMSolverIVP_PropagateBCThroughRigidLink(xf, Residual, RigidSegmentLength, actno, MagMoment[actno], CoilAlignmentTurnAreaMatrix[actno], B0, ustar[fsegi], K[fsegi], Zero3, Identity3x3, xi, Residual);
			}
			else
				CRMSolverIVP_PropagateBCThroughRigidLink(xf, Residual, RigidSegmentLength, actno, MagMoment[actno], CoilAlignmentTurnAreaMatrix[actno], B0, ustar[fsegi], K[fsegi], ustar[fsegip1], Kinv[fsegip1], xi, Residual);

			if (!FinalValueOnly) {
				// are there any localization markers?  If so, calculate their positions
				CalculateLocMarkers(NextLocMarker, xf, LocMarkers, SegBounds[i + 1], p_atLocMarkers);
			}

		}
		xi = xf; // The calculated final values will be the initial value of the next iteration
	}

	// Copy marker locations to the output
	//    note that the order is being reversed
	if (!FinalValueOnly) {
		for (int i=0; i<NUM_LOCALIZATION_MARKERS; i++) {
			for (int j=0; j<3; j++) {
				out_p_atLocMarkers[i][j]=p_atLocMarkers[(NUM_LOCALIZATION_MARKERS-1)-i][j];
			}
		}
	}

	// no need to copy the final values of the state x_f and residual to the output
	// since we have created an alias
	//out_x_N = xf;
	//mCopy_AB<3>(Residual, out_MomentResidual);

}



template<typename adType>
void CRMSolverIVP_PropagateBCThroughRigidLink(StateVector<adType>& xi_ip1, adType Residual_ip1[3],
	const adType RigidSegmentLength, const unsigned int ActNo, const adType MagMoment[3], const double CoilAlignmentTurnAreaMatrix[9],
	const double B0[3], const double ustar_i[3], const double K_i[9], const double ustar_ip1[3], const double Kinv_ip1[9],
	const StateVector<adType>& xf_i, const adType Residual_i[3]) {

	adType muhat[9], RscTB0[3], Tb[3], K2invResidual[3]; // intermediate variables

	// p
	for (int j = 0; j < 3; j++) {
		xi_ip1._p[j] = xf_i._p[j] + xf_i._R[j * 3 + 2] * RigidSegmentLength;
	}
	// R
	mCopy_AB<9>(xf_i._R, xi_ip1._R);
	// u
	// Tb=\mu_c \cross R_sc^T B0,s
	wHat<adType>(MagMoment, muhat);
	mMult_ATB<3, 3, 1>(xf_i._R, B0, RscTB0);
	mMult_AB<3, 3, 1>(muhat, RscTB0, Tb);
	// u2=u2star + ( K2inv K1 (u1 - u1star ) - K2inv Tb )
	mSub_AB<3, 1>(Residual_i, Tb, Residual_ip1);  // Residual_ip1 = K1 (u1 - u1star) - Tb = Residual_i - Tb
	mMult_AB<3, 3, 1>(Kinv_ip1, Residual_ip1, K2invResidual);
	mAdd_AB<3, 1>(ustar_ip1, K2invResidual, xi_ip1._u);

}


template <typename adType>
void CRMSolverIVP_Return (  const StateVector<adType>& in_x_N, adType in_MomentResidual[3],
							double in_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3],
							double out_x_N[NUM_STATES], double out_MomentResidual[3],
							double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]) {

	// Copy local variables to the output variables
	for (int i = 0; i < 3; i++) out_x_N[i] = dVal(in_x_N._p[i]);
	for (int i = 0; i < 9; i++) out_x_N[i + 3] = dVal(in_x_N._R[i]);
	for (int i = 0; i < 3; i++) out_x_N[i + 3 + 9] = dVal(in_x_N._u[i]);
	for (int i = 0; i < 3; i++) out_MomentResidual[i] = dVal(in_MomentResidual[i]);
	mCopy_ABm<NUM_LOCALIZATION_MARKERS,3>(in_p_atLocMarkers,out_p_atLocMarkers);

}


template<typename adType>
void CalculateLocMarkers(	int& NextLocMarker, const StateVector<adType>& xnext,
							const double LocMarkers[NUM_LOCALIZATION_MARKERS], adType SegBounds_ip1,
							double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]) {

	bool loopcondition = ((NextLocMarker < NUM_LOCALIZATION_MARKERS) && (LocMarkers[NextLocMarker] <= SegBounds_ip1));
	adType tempadType;				// intermediate variable
	while (loopcondition) {
		tempadType = LocMarkers[NextLocMarker] - SegBounds_ip1;				// this would be a negative value
		LocMarkerUpdate(p_atLocMarkers[NextLocMarker], xnext, tempadType);		//  xi has already been updated, it is the end point position
		NextLocMarker++;
		loopcondition = ((NextLocMarker < NUM_LOCALIZATION_MARKERS) && (LocMarkers[NextLocMarker] <= SegBounds_ip1));
	}
}


template <typename adType>
void LocMarkerUpdate(double p[3], const StateVector<adType>& xi, adType t) {
	p[0] = xi._p[0] + t * xi._R[2];
	p[1] = xi._p[1] + t * xi._R[5];
	p[2] = xi._p[2] + t * xi._R[8];
}


template <typename adType>
void LocMarkerUpdate(double p[3], adType xi[NUM_STATES], adType t) {
	p[0] = xi[0] + t * xi[5];
	p[1] = xi[1] + t * xi[8];
	p[2] = xi[2] + t * xi[11];
}


template <typename adType>
void CRMIntegrand (	adType s, const StateVector<adType>& in_x, const adType in_Li, const double in_dlambdainv,
					const double in_K[9], const double in_Kinv[9], const double in_l[3], const double in_ustar[3], const double in_fcumlambda[NUM_FCUM_LAMBDA+1][3],
					const adType in_nL[3],
					StateDerivativeVector<adType>& out_xdot) {

	adType Length;
	double deltalambdainv;
	adType fcum[3], ustardot[3];  //  We are assuming Kdot=0.0 (K=const)

	// copy inputs and parameters to local variables
	Length = in_Li;
	deltalambdainv = in_dlambdainv;
	auto& K = in_K;
	auto& Kinv = in_Kinv;
	auto& ustar = in_ustar;
	auto& l = in_l;
	for (int i = 0; i < 3; i++) {
		ustardot[i] = 0.0; //in_ustardot[i]; // we assume ustardot=0.0 since our rest shape model is piecewise constant curvature
	}

	// for simplicity, create aliases
	auto& u = in_x._u;
	auto& R = in_x._R;
	auto& udot = out_xdot._u;
#ifndef ANALYTICAL_SE3_STEP
	auto& p = in_x._p;				// we will not need this for analytical calculation
	auto& pdot = out_xdot._p;		// we will not need this for analytical calculation
	auto& Rdot = out_xdot._R;		// we will not need this for analytical calculation
#endif	


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

    adType nL_spatial[3];
    mMult_AB<3,3,1>(R, in_nL, nL_spatial);
    // add the tip force to fcum
    for (int i = 0; i < 3; i++) {
        fcum[i]		+=	nL_spatial[i];
    }

	// calculate u_hat
	adType u_hat[9];
	wHat(u, u_hat);

	// udot = ustardot - Kinv*((um*K+Kdot)*(u-ustar_s) + e3m*R'*intf + R'*l); % udot
	//
	//   e3hat*R' = [ -r12 -r22 -r32; r11 r21 r31; 0 0 0];
	adType e3hatRT[9];
	e3hatRT[0] = -R[1];   e3hatRT[1] = -R[4];   e3hatRT[2] = -R[7];
	e3hatRT[3] = R[0];    e3hatRT[4] = R[3];    e3hatRT[5] = R[6];
	e3hatRT[6] = 0.0;     e3hatRT[7] = 0.0;     e3hatRT[8] = 0.0;
	adType e3hatRTfcum[3];
	mMult_AB<3, 3, 1>(e3hatRT, fcum, e3hatRTfcum);			//  e3m*R'*intf
	adType RTl[3];
	mMult_ATB<3, 3, 1>(R, l, RTl);								// R'*l
	adType umustar[3];
	mSub_AB<3, 1>(u, ustar, umustar);							// (u-ustar_s)
	adType Kumustar[3], uhatKumustar[3];
	mMult_AB<3, 3, 1>(K, umustar, Kumustar);
	mMult_AB<3, 3, 1>(u_hat, Kumustar, uhatKumustar); 		 	//(um*K+Kdot)*(u-ustar_s)  assuming Kdot=0
	adType sumterm[3];
	mAdd_ABC<3, 1>(uhatKumustar, e3hatRTfcum, RTl, sumterm);	// ((um*K+Kdot)*(u-ustar_s) + e3m*R'*intf + R'*l)
	adType KinvSum[3];
	mMult_AB<3, 3, 1>(Kinv, sumterm, KinvSum);					// Kinv*((um*K+Kdot)*(u-ustar_s) + e3m*R'*intf + R'*l)
	mSub_AB<3, 1>(ustardot, KinvSum, udot);					// udot = ustardot - Kinv*((um*K+Kdot)*(u-ustar_s) + e3m*R'*intf + R'*l);

#ifndef ANALYTICAL_SE3_STEP
	// we will not need these for analytical calculation
	// Rdot = R*u_hat
	mMult_AB<3, 3, 3>(R, u_hat, Rdot);
	// pdot = R*e3,
	for (int i = 0; i < 3; i++) {
		pdot[i] = R[i * 3 + 2];
	}
#endif

	//xdot(1:3) = R*e3;             % pdot
	//xdot(4:12) = reshape(R*um,9,1); % Rdot   --- Note that matlab code reshapes in column major order while we are saving in row major order
	//xdot(13:15) = ustardot - Kinv*((um*K+Kdot)*(u-ustar_s) + e3m*R'*intf + R'*l); % udot
	//  note: the sample code has matlab indexing starting from 1 to 15

}


//
//
//	NUMERICAL INTEGRATION SUPPORT FUNCTIONS
//
//


template <typename adType>
void Project_State_to_Manifold(StateVector<adType>& State) {
#ifdef ANALYTICAL_SE3_STEP

	// we don't need to do anything for analytical SE3 step

#else

	// project R to SO(3)
#error("Functionality not implemented!...\n");

#endif
}



//
//
//	NUMERICAL INTEGRATION FUNCTIONS
//
//

//function [x_1toN] = ABM4(x_0, t_0, N, h, Integrand, initmethod)
template <typename adType, typename StVecType>
void ABM4 (	const StVecType& in_x_0, const adType t_0, const int N, const adType h,
			const adType Li, const double dlambdainv, const double in_K[9], const double in_Kinv[9], const double in_l[3], 
			const double in_ustar[3], const double in_fcumlambda[NUM_FCUM_LAMBDA+1][3], const adType in_nL[3],
			const bool FinalValueOnly, const double in_LocMarkers[NUM_LOCALIZATION_MARKERS], int &inout_NextLocMarkerIdx,
			StVecType& out_x_N, double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]	) {

	using _SVT = StVecType;
	using _DVT = StDerivativeVectType<StVecType, adType>;

	_SVT x_nm3;
	_SVT x_nm2;
	_SVT x_nm1;
	_SVT x_n(in_x_0);
	_SVT x_np1;
	adType t_n;
	_DVT xdot_nm3;
	_DVT xdot_nm2;
	_DVT xdot_nm1;
	_DVT xdot_n;
	bool loopcondition;

	int NextLocMarkerIdx = inout_NextLocMarkerIdx;
	double LocMarkers[NUM_LOCALIZATION_MARKERS];
	double delta_n_overh,delta_np1_overh;

	// initialize the iteration items
	t_n = t_0;
	for (int i=0; i<NUM_LOCALIZATION_MARKERS; i++) {
		LocMarkers[i]=in_LocMarkers[i];
	}

	for (int idx=0; idx<N; idx++) {

		if (idx<3) {  // RK2 initialization steps
			RK2_step(x_n, t_n, h, Li, dlambdainv, in_K, in_Kinv, in_l, in_ustar, in_fcumlambda, in_nL, x_np1, xdot_n);
		}
		else { 		 // ABM4 steps
			ABM4_step(x_n, t_n, h, xdot_nm1, xdot_nm2, xdot_nm3, x_nm1, x_nm2, x_nm3, Li, dlambdainv, in_K, in_Kinv, in_l, in_ustar, in_fcumlambda, in_nL, x_np1, xdot_n);
		}

		//Project_State_to_Manifold(x_np1);
		// increment "time"
		t_n=t_n+h;

		if (!FinalValueOnly) {
			// are there any localization markers?  If so, calculate their positions
			// remember, t_n has already been incremented
			loopcondition = (NextLocMarkerIdx<NUM_LOCALIZATION_MARKERS ) && (LocMarkers[NextLocMarkerIdx]<=t_n);
			while ( loopcondition ) {
				delta_np1_overh=(dVal(t_n)-LocMarkers[NextLocMarkerIdx])/dVal(h);
				delta_n_overh=1.0-delta_np1_overh;
				out_p_atLocMarkers[NextLocMarkerIdx][0]=delta_np1_overh*dVal(x_n._p[0])+delta_n_overh*dVal(x_np1._p[0]);
				out_p_atLocMarkers[NextLocMarkerIdx][1]=delta_np1_overh*dVal(x_n._p[1])+delta_n_overh*dVal(x_np1._p[1]);
				out_p_atLocMarkers[NextLocMarkerIdx][2]=delta_np1_overh*dVal(x_n._p[2])+delta_n_overh*dVal(x_np1._p[2]);
				NextLocMarkerIdx++;
				loopcondition = (NextLocMarkerIdx<NUM_LOCALIZATION_MARKERS ) && (LocMarkers[NextLocMarkerIdx]<=t_n);
			}
		}

		// update the iteration items
		x_nm3 = x_nm2;
		x_nm2 = x_nm1;
		x_nm1 = x_n;
		x_n = x_np1;
		xdot_nm3 = xdot_nm2;
		xdot_nm2 = xdot_nm1;
		xdot_nm1 = xdot_n;

		//for (int i = 0; i < NUM_STATES; i++) {
		//	x_nm3[i] = x_nm2[i];
		//	x_nm2[i] = x_nm1[i];
		//	x_nm1[i] = x_n[i];
		//	x_n[i] = x_np1[i];
		//}
		//for (int i = 0; i < NUM_INTEGRATION_STATES; i++) {
		//	xdot_nm3[i]=xdot_nm2[i];
		//	xdot_nm2[i]=xdot_nm1[i];
		//	xdot_nm1[i]=xdot_n[i];
		//}

	}

	// Copy final value
	out_x_N=x_n;
	inout_NextLocMarkerIdx=NextLocMarkerIdx;

}


// [x_np1, xdot_n, xdot_nm1, xdot_nm2] = ABM4_step(x_n, t_n, xdot_nm1, xdot_nm2, xdot_nm3, h, Integrand)
template <typename adType, typename StVecType>
void ABM4_step(	const StVecType& in_x_n, adType t_n, adType h,
				const StDerivativeVectType<StVecType, adType>& in_xdot_nm1, const StDerivativeVectType<StVecType, adType>& in_xdot_nm2, const StDerivativeVectType<StVecType, adType>& in_xdot_nm3,
				const StVecType& in_x_nm1, const StVecType& in_x_nm2, const StVecType& in_x_nm3,
				const adType Li, const double dlambdainv, const double in_K[9], const double in_Kinv[9], const double in_l[3], 
				const double in_ustar[3], const double in_fcumlambda[NUM_FCUM_LAMBDA + 1][3], const adType in_nL[3],
				StVecType& out_x_np1, StDerivativeVectType<StVecType, adType>& out_xdot_n) {

	using _SVT = StVecType;
	using _DVT = StDerivativeVectType<StVecType, adType>;

	const double P_COEFF_N=55.0/24.0, P_COEFF_Nm1=-59.0/24.0, P_COEFF_Nm2=37.0/24.0, P_COEFF_Nm3=-9.0/24.0;  // AB4 Predictor Coefficients
	const double C_COEFF_Np1=9.0/24.0, C_COEFF_N=19.0/24.0, C_COEFF_Nm1=-5.0/24.0, C_COEFF_Nm2=1.0/24.0;     // AM4 Corrector Coefficients
	_SVT x_n(in_x_n);       		// from input
	_SVT x_nm1(in_x_nm1);       	// from input
	_SVT x_nm2(in_x_nm2);       	// from input
	_SVT x_nm3(in_x_nm3);       	// from input
	_SVT x_np1_hat;    				// intermediate
	_DVT xdot_np1_hat;				// intermediate
	auto& xdot_n = out_xdot_n;
	_DVT xdot_nm1(in_xdot_nm1);  	// from input
	_DVT xdot_nm2(in_xdot_nm2);  	// from input
	_DVT xdot_nm3(in_xdot_nm3);  	// from input


	//ABM4_STEP_STEP1:
	CRMIntegrand(t_n, x_n, Li, dlambdainv, in_K, in_Kinv, in_l, in_ustar, in_fcumlambda, in_nL, xdot_n);
	x_np1_hat    = x_n + h * ( P_COEFF_N * xdot_n + P_COEFF_Nm1 * xdot_nm1 + P_COEFF_Nm2 * xdot_nm2 + P_COEFF_Nm3 * xdot_nm3 );
#ifdef ANALYTICAL_SE3_STEP
	//      calculate R_np1_hat and p_np1_hat analytically, without numerical integration
	adType u_n_pred[3];
	for (int i = 0; i < 3; i++) u_n_pred[i] = (P_COEFF_N * x_n._u[i] + P_COEFF_Nm1 * x_nm1._u[i] + P_COEFF_Nm2 * x_nm2._u[i] + P_COEFF_Nm3 * x_nm3._u[i]);
	SE3_Analytical_Step(x_n._R, x_n._p, u_n_pred, h, x_np1_hat._R /*R_np1_hat*/, x_np1_hat._p /*p_np1_hat*/);
#endif
	//ABM4_STEP_STEP2:
	CRMIntegrand(t_n+h, x_np1_hat, Li, dlambdainv, in_K, in_Kinv, in_l, in_ustar, in_fcumlambda, in_nL, xdot_np1_hat);
	out_x_np1    = x_n + h * ( C_COEFF_Np1 * xdot_np1_hat + C_COEFF_N * xdot_n + C_COEFF_Nm1 * xdot_nm1 + C_COEFF_Nm2 * xdot_nm2 );
#ifdef ANALYTICAL_SE3_STEP
	//      calculate R_np1 and p_np1 analytically, without numerical integration
	adType u_n_corr[3];
	for (int i = 0; i < 3; i++) u_n_corr[i] = (C_COEFF_Np1 * x_np1_hat._u[i] + C_COEFF_N * x_n._u[i] + C_COEFF_Nm1 * x_nm1._u[i] + C_COEFF_Nm2 * x_nm2._u[i]);
	SE3_Analytical_Step(x_n._R, x_n._p, u_n_corr, h, out_x_np1._R /*R_np1*/, out_x_np1._p /*p_np1*/);
#endif

}


//[x_np1, xdot_n] = RK2_step(x_n, t_n, h, Integrand)
template <typename adType, typename StVecType>
void RK2_step(	const StVecType& in_x_n, const adType t_n, const adType h,
				const adType Li, const double dlambdainv, const double in_K[9], const double in_Kinv[9], const double in_l[3], 
				const double in_ustar[3], const double in_fcumlambda[NUM_FCUM_LAMBDA+1][3], const adType in_nL[3],
				StVecType& out_x_np1, StDerivativeVectType<StVecType, adType>& out_xdot_n ) {

	using _SVT = StVecType;
	using _DVT = StDerivativeVectType<StVecType, adType>;

	_SVT x_n(in_x_n);       // from input
	_DVT k1;				// intermediate
	_DVT k2oh;				// intermediate
	_SVT x_n_p_k1o2;		// intermediate
	auto& xdot_n = out_xdot_n;// for output

	//RK2_STEP_STEP1:
	CRMIntegrand(t_n, x_n, Li, dlambdainv, in_K, in_Kinv, in_l, in_ustar, in_fcumlambda, in_nL, xdot_n);
	k1			= h * xdot_n;
	x_n_p_k1o2 	= x_n + k1 * 0.5;
#ifdef ANALYTICAL_SE3_STEP
	// we will calculate R_n_p_k1o2 and p_n_p_k1o2 analytically, without numerical integration
	SE3_Analytical_Step(x_n._R, x_n._p, x_n._u, h * 0.5, x_n_p_k1o2._R, x_n_p_k1o2._p);
#endif

	//RK2_STEP_STEP2:
	CRMIntegrand(t_n + h * 0.5, x_n_p_k1o2, Li, dlambdainv, in_K, in_Kinv, in_l, in_ustar, in_fcumlambda, in_nL, k2oh);
	out_x_np1	= x_n + h * k2oh;
#ifdef ANALYTICAL_SE3_STEP
	// we will calculate R_np1 and p_np1 analytically, without numerical integration
	SE3_Analytical_Step(x_n._R, x_n._p, x_n_p_k1o2._u/*u_np1half*/, h, out_x_np1._R, out_x_np1._p);
#endif

}


// definitions needed for twist exponential calculation
#define EPS 1.0e-12   // the threshold for assuming ||u|| to be approximately 0, so that we should use pure translation equation


// calculate R_np1 and p_np1 analytically using twist exponential, without numerical integration
//   g_np1 = g_n * expm ( \hat{\xi}^b *h ),  where \xi^b= [ 0 0 1 u_n^T ]^T
//   g = [R p; 0 0 0 1];
template <typename adType>
void SE3_Analytical_Step(adType in_R_n[9], adType in_p_n[3], adType in_u_n[3], adType h, adType out_R_np1[9], adType out_p_np1[3]) {
	adType R_n[9], p_n[3], u_n[3];

	mCopy_AB<3*3>(in_R_n, R_n);
	mCopy_AB<3>(in_p_n, p_n);
	mCopy_AB<3>(in_u_n, u_n);

	// calculate R_np1 and p_np1 analytically, without numerical integration
	adType Rdelta[9], pdelta[3];
	adType umagsq, umagsqresp, umag, umagresp, unorm[3], delsumag, uu3dels[3], ImRuxv[3], ImRuxvpuuTvds[3], Rnpd[3];
	umagsq = vNormSq<3>(u_n);
	if (umagsq < EPS) {
		mCopy_AB<3*3>(R_n, out_R_np1);
		out_p_np1[0] = p_n[0] + R_n[2] * h;
		out_p_np1[1] = p_n[1] + R_n[5] * h;
		out_p_np1[2] = p_n[2] + R_n[8] * h;
	}
	else {
		umagsqresp = 1.0 / umagsq;
		umag = sqrt(umagsq);
		umagresp = 1.0 / umag;
		mMult_sA<3, 1>(umagresp, u_n, unorm);
		delsumag = h * umag;
		RodriguesExpanded(unorm, delsumag, Rdelta);
		mMult_sA<3, 1>(u_n[2] * h, u_n, uu3dels);
		ImRuxv[0] = Rdelta[0 * 3 + 1] * u_n[0] - Rdelta[0 * 3 + 0] * u_n[1] + u_n[1];
		ImRuxv[1] = Rdelta[1 * 3 + 1] * u_n[0] - Rdelta[1 * 3 + 0] * u_n[1] - u_n[0];
		ImRuxv[2] = Rdelta[2 * 3 + 1] * u_n[0] - Rdelta[2 * 3 + 0] * u_n[1];
		mAdd_AB<3, 1>(ImRuxv, uu3dels, ImRuxvpuuTvds);
		mMult_sA<3, 1>(umagsqresp, ImRuxvpuuTvds, pdelta);
		mMult_AB<3, 3, 3>(R_n, Rdelta, out_R_np1);
		mMult_AB<3, 3, 1>(R_n, pdelta, Rnpd);
		mAdd_AB<3, 1>(p_n, Rnpd, out_p_np1);
	}
}

#undef EPS

//
//
//  Robotic Kinematics Related Functions
//
//

template <typename adType>
void wHat(const adType in_w[3], adType out_what[9]) {  // what is stored as a 1-dim array in row major order

	out_what[0] = 0.0;
	out_what[1] = -in_w[2];
	out_what[2] = in_w[1];
	out_what[3] = in_w[2];
	out_what[4] = 0.0;
	out_what[5] = -in_w[0];
	out_what[6] = -in_w[1];
	out_what[7] = in_w[0];
	out_what[8] = 0.0;

}


template <typename adType>
void wHat(const adType in_w[3], adType out_what[9], unsigned int stride) {  // what is stored as a 1-dim array in row major order

	out_what[0] = 0.0;
	out_what[1] = -in_w[2*stride];
	out_what[2] = in_w[1*stride];
	out_what[3] = in_w[2*stride];
	out_what[4] = 0.0;
	out_what[5] = -in_w[0*stride];
	out_what[6] = -in_w[1*stride];
	out_what[7] = in_w[0*stride];
	out_what[8] = 0.0;

}


template <typename adType>
void vee_from_so3(const adType in_what[9], adType out_w[3]) {

	out_w[0] = 0.5 * (in_what[2 * 3 + 1] - in_what[1 * 3 + 2]);
	out_w[1] = 0.5 * (in_what[0 * 3 + 2] - in_what[2 * 3 + 0]);
	out_w[3] = 0.5 * (in_what[1 * 3 + 0] - in_what[0 * 3 + 1]);

}



template <typename adType>
void RodriguesFormula(adType in_w[3], adType in_theta, adType out_R[9]) {

	adType I3x3[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

	adType w[3], what[9], whatsq[9], temp1[9], temp2[9], temp3[9];
	adType c = cos(in_theta);
	adType s = sin(in_theta);
	adType v = 1.0 - c;

	for (int i = 0; i < 3; i++) {
		w[i] = in_w[i];
	}

	wHat(w, what);
	mMult_sA<3, 3>(s, what, temp1);
	mMult_AB<3, 3, 3>(what, what, whatsq);
	mMult_sA<3, 3>(v, whatsq, temp2);
	mAdd_AB<3, 3>(I3x3, temp1, temp3);
	mAdd_AB<3, 3>(temp2, temp3, out_R);

}


template <typename adType>
void RodriguesExpanded(adType in_w[3], adType in_theta, adType out_R[9]) {

	adType w1 = in_w[0];
	adType w2 = in_w[1];
	adType w3 = in_w[2];
	adType c = cos(in_theta);
	adType s = sin(in_theta);
	adType v = 1.0 - c;
	adType w1s = w1 * s;
	adType w2s = w2 * s;
	adType w3s = w3 * s;
	adType w1w2v = w1 * w2 * v;
	adType w1w3v = w1 * w3 * v;
	adType w2w3v = w2 * w3 * v;
	out_R[0 * 3 + 0] = w1 * w1 * v + c;
	out_R[0 * 3 + 1] = w1w2v - w3s;
	out_R[0 * 3 + 2] = w1w3v + w2s;
	out_R[1 * 3 + 0] = w1w2v + w3s;
	out_R[1 * 3 + 1] = w2 * w2 * v + c;
	out_R[1 * 3 + 2] = w2w3v - w1s;
	out_R[2 * 3 + 0] = w1w3v - w2s;
	out_R[2 * 3 + 1] = w2w3v + w1s;
	out_R[2 * 3 + 2] = w3 * w3 * v + c;

}
