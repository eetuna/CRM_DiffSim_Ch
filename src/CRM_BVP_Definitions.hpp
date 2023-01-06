#pragma once
#include <cmath>
#include "numerical/minpack_Declarations.hpp"

#define M_PI 3.14159265358979323846

template <typename adType>
void CRMShootingMethodBVP(	CRMShootingMethodParams<adType> in_Params, 
							double in_u0_initialguess[3], double in_ftip_initialguess[3],
							adType out_u0[3], adType out_ftip[3], int& out_localmin) {

	ContactModeType ContactMode = in_Params.ContactMode;
	int NLEq_Dim;  // Dimension of the Nonlinear Equation to Solve
	if (ContactMode == ContactModeType::FREE_TIP) {
		NLEq_Dim = 3;
	}
	else { // FIXED_TIP
		NLEq_Dim = 6;
	}
	// Call CRMSolverIVP_Prep, to pre-process parameters
	adType x_0[NUM_STATES];
	for (int i = 0; i < 3; i++) x_0[i] = in_Params.p0[i];
	for (int i = 0; i < 9; i++) x_0[i + 3] = in_Params.R0[i];
	for (int i = 0; i < 3; i++) x_0[i + 3 + 9] = in_u0_initialguess[i];
	bool FinalValueOnly = true;
	NLEqnParams<adType> NLEParams;

	CRMSolverIVP_Prep(x_0, in_Params.IntegrationStepSize,
		in_Params.Li, in_Params.dlambdainv,
		in_Params.SegEndLambdas, in_Params.LocMarkerLambdas,
		in_Params.K, in_Params.Kinv, in_Params.ustar,
		in_Params.CoilAlignmentTurnAreaMatrix,
		in_Params.MagMoment, in_Params.fcumlambda,
		in_Params.B0,
		FinalValueOnly,
		NLEParams);
	NLEParams.ContactMode = ContactMode;
	mCopy_AB<3>(in_Params.TipConstraintPoint, NLEParams.TipConstraintPoint);

	// Scale parameters and call the nonlinear equation solver
	const double USCALE_INV = 1.0 / IVALUE_SCALE_U;
	const double FSCALE_INV = 1.0 / IVALUE_SCALE_F;
	double* initialguessscaled = new double [NLEq_Dim];
	adType* returnedparamscaled = new adType [NLEq_Dim];
	if (ContactMode == ContactModeType::FREE_TIP) {
		for (int i = 0; i < 3; i++) {
			initialguessscaled[i] = USCALE_INV * in_u0_initialguess[i];
			NLEParams.TipForce[i] = in_Params.TipForce[i];  // if the catheter is not in contact, the tip force specified within in_Params needs to be used; this would not be scaled as it is not changed by the solver
		}
	}
	else { // FIXED_TIP
		for (int i = 0; i < 3; i++) {
			initialguessscaled[i] = USCALE_INV * in_u0_initialguess[i];
			initialguessscaled[i+3] = FSCALE_INV * in_ftip_initialguess[i];
		}

	}

	int localmin = 0;

#if defined( FK_TRUSTREGION )
	adType* x = new adType[NLEq_Dim]; // we will create a new variable here and not use initial guess scaled since truss-region-dogleg algorithm uses the same variable for both input and output
	adType* residual = new adType[NLEq_Dim];
	int info;
	double tol = TRUSTREGION_TOLERANCE;
	for (int i = 0; i < NLEq_Dim; i++) x[i] = initialguessscaled[i];

#if defined  (FK_TRUSTREGION_ANALYTICALJAC)
	TrustRegionDogleg_GivenJacobian<adType, NLEqnParams<adType>>(CRM_NLEquation<adType>, CRM_NLEquation_AnalyticalJac<adType>, 3, x, residual, tol, info, NLEParams);
#else 	
	TrustRegionDogleg(CRM_NLEquation<adType>, NLEq_Dim, x, residual, tol, info, NLEParams);
#endif

	localmin = (info == 1) ? 0 : (info - 1);
	for (int i = 0; i < NLEq_Dim; i++) returnedparamscaled[i] = x[i];
	delete[] residual;
	delete[] x;
#else // undefined
	exit(1);
#endif

	/*
	// AUTODIFF TEST CODE
	VectorXadType u(NLEq_Dim), y(NLEq_Dim);
	for (int i=0;i<NLEq_Dim;i++) u(i)=x[i];
	MatrixXd J = jacobian(NLEquationAD, wrt(u), at(u, NLEParams), y);
	std::cout << "y = \n" << y << std::endl;    // print the evaluated output vector F
	std::cout << "J = \n" << J << std::endl;    // print the evaluated Jacobian matrix dF/dx
	*/

	if (ContactMode == ContactModeType::FREE_TIP) {
		for (int i = 0; i < 3; i++) {
			out_u0[i] = IVALUE_SCALE_U * returnedparamscaled[i];
			out_ftip[i] = in_Params.TipForce[i];  // if it is free-tip, return the tip force specified within in_Params
		}
	}
	else { // FIXED_TIP
		for (int i = 0; i < 3; i++) {
			out_u0[i] = IVALUE_SCALE_U * returnedparamscaled[i];
			out_ftip[i] = IVALUE_SCALE_F * returnedparamscaled[i+3];
		}
	}
	out_localmin = localmin;
	delete[] initialguessscaled;
	delete[] returnedparamscaled;
}


template <typename adType>
void CRMShootingMethodBVP(
	ContactModeType in_ContactMode,					// in_ContactMode == FREE_TIP if the catheter is not in contact with a surface, FIXED_TIP if catheter tip is constrained to be at TipContraintPoint
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
	int& out_localmin							// out_localmin!=0 if algorithms is stuck at a local minimum, or cannot make further progress
) {

	// Copy input parameters to local variables --- needed for dataflow optimization
	CRMCatheterModelParams CathParams;
	CatheterConfiguration CathConfig;
	CRMShootingMethodBVP_Prep(in_B0, in_g, in_p0, in_R0,
		in_SegLengths, in_LocMarkers, in_InnerRadius, in_OuterRadius, in_YoungsModulus, in_ShearModulus,
		in_ustar, in_CoilAlignmentAngles, in_CoilTurnAreaMat, in_rho, in_ActMass,
		CathParams, CathConfig);

	adType InsertedLength = in_InsertedLength;
	adType ActuationCurrents[NUM_ACT_SET][3];
	mCopy_AB<NUM_ACT_SET * 3>(&(in_ActuationCurrents[0][0]), &(ActuationCurrents[0][0]));
	double IntegrationStepSize = in_IntegrationStepSize;
	// Calculate Shooting Method Parameter Set from model and configuration parameters
	CRMShootingMethodParams<adType> ShootingParams;
	CRMConstructShootingMethodParamSet(CathParams, CathConfig,
		InsertedLength, ActuationCurrents, in_ContactMode, in_TipConstraintPoint, in_TipForce,
		IntegrationStepSize,
		ShootingParams);

	adType u0_calc[3], ftip_calc[3];

	CRMShootingMethodBVP(ShootingParams, in_u0_initialguess, in_ftip_initialguess, u0_calc, ftip_calc, out_localmin);

	for (int i = 0; i < 3; i++) out_u0[i] = u0_calc[i];
	for (int i = 0; i < 3; i++) out_ftip[i] = ftip_calc[i];

}


template <typename adType>
void CRM_NLEquation(adType in_x[], adType out_y[], NLEqnParams<adType> Params) {

	StateVector<adType> x_N;
	double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3];
	adType MomentResidual[3];

	adType u_0[3], ftip[3];
	// don't forget to scale parameters before passing to the CRMSolverIVP
	if (Params.ContactMode == ContactModeType::FREE_TIP) {
		for (int i = 0; i < 3; i++) {
			u_0[i] = IVALUE_SCALE_U * in_x[i];
			ftip[i] = Params.TipForce[i];  // for free-tip, this parameter is not given by the nonlinear equation solver, and hence, does not need to be scaled
		}
	}
	else { // FIXED_TIP
		for (int i = 0; i < 3; i++) {
			u_0[i] = IVALUE_SCALE_U * in_x[i];
			ftip[i] = IVALUE_SCALE_F * in_x[i + 3];
		}
	}
	// We will only call the IVP_Core, since preprocessing is already done
	CRMSolverIVP_Core(Params, u_0, ftip, x_N, MomentResidual, p_atLocMarkers);

	// don't forget to scale parameters before returning to the nonlinear equation solver
	if (Params.ContactMode == ContactModeType::FREE_TIP) {
		for (int i = 0; i < 3; i++) {
			out_y[i] = RESIDUAL_SCALE_M * MomentResidual[i];
		}
	}
	else { // FIXED_TIP
		for (int i = 0; i < 3; i++) {
			out_y[i] = RESIDUAL_SCALE_M * MomentResidual[i];
			out_y[i + 3] = RESIDUAL_SCALE_P * (x_N._p[i] - Params.TipConstraintPoint[i]);
		}
	}

}


template<typename adType>
void CRM_NLEquation_AnalyticalJac(adType in_x[], adType out_y[], adType out_fjac[], NLEqnParams<adType> Params) {

	AugmentedStateVector<adType, IVPJacobiansMini> x_N;
	adType MomentResidual[3], kJuu0[9];

	adType u_0[3], ftip[3];
	// don't forget to scale parameters before passing to the CRMSolverIVP
	if (Params.ContactMode == ContactModeType::FREE_TIP) {
		for (int i = 0; i < 3; i++) {
			u_0[i] = IVALUE_SCALE_U * in_x[i];
			ftip[i] = Params.TipForce[i];  // for free-tip, this parameter is not given by the nonlinear equation solver, and hence, does not need to be scaled
		}
	}
	else { // FIXED_TIP
		for (int i = 0; i < 3; i++) {
			u_0[i] = IVALUE_SCALE_U * in_x[i];
			ftip[i] = IVALUE_SCALE_F * in_x[i + 3];
		}
	}
	// We will only call the IVP_Core, since preprocessing is already done
	CRMSolverIVP_CoreWithJacobian(Params, u_0, ftip, x_N, MomentResidual);

	// don't forget to scale parameters before returning to the nonlinear equation solver
	if (Params.ContactMode == ContactModeType::FREE_TIP) {
		for (int i = 0; i < 3; i++) {
			out_y[i] = RESIDUAL_SCALE_M * MomentResidual[i];
		}
		if (NUM_SEGMENTS % 2 == 1) { // if last segment is flexible
			mMult_AB<3, 3, 3>(Params.K[NUM_FLEX_SEG - 1], x_N._u_u0, kJuu0);
		}
		else {  // if last segment is rigid
			mCopy_AB<3 * 3>(x_N._u_u0, kJuu0);
		}
		for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) out_fjac[i + j * 3] = RESIDUAL_SCALE_M * kJuu0[i * 3 + j];  // Note that minpack uses Fortran style column-major ordering, not C sytle row-major ordering

	}
	else { // FIXED_TIP
		for (int i = 0; i < 3; i++) {
			out_y[i] = RESIDUAL_SCALE_M * MomentResidual[i];
			out_y[i + 3] = RESIDUAL_SCALE_P * (x_N._p[i] - Params.TipConstraintPoint[i]);
		}
		std::cerr << "Functionality Not Implemented!...\n";
		exit(1);
	}


}



void CRMShootingMethodBVP_Prep (
		double in_B0[3], double in_g[3], double in_p0[3], double in_R0[9],
		double in_SegLengths[NUM_SEGMENTS], double in_LocMarkers[NUM_LOCALIZATION_MARKERS],
		double in_InnerRadius[NUM_FLEX_SEG],	double in_OuterRadius[NUM_FLEX_SEG],
		double in_YoungsModulus[NUM_FLEX_SEG], double in_ShearModulus[NUM_FLEX_SEG],
		double in_ustar[NUM_FLEX_SEG][3],
		double in_CoilAlignmentAngles[NUM_ACT_SET][2], double in_CoilTurnAreaMat[NUM_ACT_SET][9],
		double in_rho[NUM_SEGMENTS],	double in_ActMass[NUM_ACT_SET],
		CRMCatheterModelParams &CathParams, CatheterConfiguration &CathConfig ) {

	mCopy_AB<NUM_SEGMENTS>(in_SegLengths,CathParams.SegLengths);
	mCopy_AB<NUM_LOCALIZATION_MARKERS>(in_LocMarkers,CathParams.LocMarkers);
	mCopy_AB<NUM_FLEX_SEG>(in_InnerRadius,CathParams.InnerRadius);
	mCopy_AB<NUM_FLEX_SEG>(in_OuterRadius,CathParams.OuterRadius);
	mCopy_AB<NUM_FLEX_SEG>(in_YoungsModulus,CathParams.YoungsModulus);
	mCopy_AB<NUM_FLEX_SEG>(in_ShearModulus,CathParams.ShearModulus);
	mCopy_ABm<NUM_FLEX_SEG,3>(in_ustar,CathParams.ustar);
	mCopy_ABm<NUM_ACT_SET,2>(in_CoilAlignmentAngles,CathParams.CoilAlignmentAngles);
	mCopy_ABm<NUM_ACT_SET,9>(in_CoilTurnAreaMat,CathParams.CoilTurnAreaMat);
	mCopy_AB<NUM_SEGMENTS>(in_rho,CathParams.rho);
	mCopy_AB<NUM_ACT_SET>(in_ActMass,CathParams.ActMass);
	mCopy_AB<3>(in_B0,CathConfig.B0);
	mCopy_AB<3>(in_g,CathConfig.g);
	mCopy_AB<3>(in_p0,CathConfig.p0);
	mCopy_AB<9>(in_R0,CathConfig.R0);

}

template <typename adType>
void CRMConstructShootingMethodParamSet(	CRMCatheterModelParams CathParams, CatheterConfiguration CathConfig,
											adType InsertionLength, adType ActuationCurrents[NUM_ACT_SET][3],
											ContactModeType ContactMode,
											double TipConstraintPoint[3], adType TipForce[3],
											double IntegrationStepSize,
											CRMShootingMethodParams<adType> &ShootingParams) {
	double length = 0.0;
	for (int i = 0; i < NUM_SEGMENTS; i++) {
		length += CathParams.SegLengths[i];
		ShootingParams.SegEndLambdas[i] = length;
	}
	ShootingParams.Li = MIN(InsertionLength, length);
	double dlambda = length / NUM_FCUM_LAMBDA;
	ShootingParams.dlambdainv = 1.0 / dlambda;
	mCopy_AB<NUM_LOCALIZATION_MARKERS>(CathParams.LocMarkers, ShootingParams.LocMarkerLambdas);
	double oR, iR, mI, pmI, E, G, K[9], Kinv[9];
	for (int i = 0; i < NUM_FLEX_SEG; i++) {
		oR = CathParams.OuterRadius[i];
		iR = CathParams.InnerRadius[i];
		E = CathParams.YoungsModulus[i];
		G = CathParams.ShearModulus[i];
		mI = 0.25 * M_PI * (POW4(oR) - POW4(iR));	// Area moment of inertia along x - axis
		pmI = 0.5 * M_PI * (POW4(oR) - POW4(iR));	// Polar moment of inertia of area
		K[0] = E * mI;		K[1] = 0.0;			K[2] = 0.0;
		K[3] = 0.0;			K[4] = E * mI;		K[5] = 0.0;
		K[6] = 0.0;			K[7] = 0.0;			K[8] = G * pmI;
		Kinv[0] = 1.0 / (E * mI);		Kinv[1] = 0.0;				Kinv[2] = 0.0;
		Kinv[3] = 0.0;				Kinv[4] = 1.0 / (E * mI);		Kinv[5] = 0.0;
		Kinv[6] = 0.0;				Kinv[7] = 0.0;				Kinv[8] = 1.0 / (G * pmI);
		mCopy_AB<9>(K, ShootingParams.K[i]);
		mCopy_AB<9>(Kinv, ShootingParams.Kinv[i]);
	}
	mCopy_ABm<NUM_FLEX_SEG, 3>(CathParams.ustar, ShootingParams.ustar);
	double CoilAlignMat[9], c0, s0, c1, s1;
	adType tempf[3];
	for (int i = 0; i < NUM_ACT_SET; i++) {
		c0 = cos(CathParams.CoilAlignmentAngles[i][0]);
		s0 = sin(CathParams.CoilAlignmentAngles[i][0]);
		c1 = cos(CathParams.CoilAlignmentAngles[i][1]);
		s1 = sin(CathParams.CoilAlignmentAngles[i][1]);
		mMult_AB<3, 3, 1>(CathParams.CoilTurnAreaMat[i], ActuationCurrents[i], tempf);
		CoilAlignMat[0] = c0;	CoilAlignMat[1] = -s1;	CoilAlignMat[2] = 0.0;
		CoilAlignMat[3] = s0;	CoilAlignMat[4] = c1;	CoilAlignMat[5] = 0.0;
		CoilAlignMat[6] = 0.0;	CoilAlignMat[7] = 0.0;	CoilAlignMat[8] = 1.0;
		mMult_AB<3, 3, 1>(CoilAlignMat, tempf, ShootingParams.MagMoment[i]);
		mMult_AB<3, 3, 3>(CoilAlignMat, CathParams.CoilTurnAreaMat[i], ShootingParams.CoilAlignmentTurnAreaMatrix[i]);
	}

	ShootingParams.fcumlambda[0][0] = 0.0;
	ShootingParams.fcumlambda[0][1] = 0.0;
	ShootingParams.fcumlambda[0][2] = 0.0;
	double cumpos = 0.0, mass = 0.0, prevactmass = 0.0, lastsegmentend, thissegmentend;
	int csegment = 0, act_off, flex_off;
	bool RigidTip = false;
	if (NUM_SEGMENTS % 2 == 0) RigidTip = true;  // if there are even number of segments, then the tip (most distal) segment is rigid, since the most proximal segment is assumed to be always flexible
	// calculation loop
	for (int i = 1; i < NUM_FCUM_LAMBDA + 1; i++) {
		mass = 0.0;
		cumpos = i * dlambda;
		while ((cumpos > ShootingParams.SegEndLambdas[csegment]) && (csegment < NUM_SEGMENTS)) csegment++;		// let's find the current segment number,
		if (csegment == 0)	lastsegmentend = 0.0;								// the end point of the previous
		else lastsegmentend = ShootingParams.SegEndLambdas[csegment - 1];		// segment, and
		thissegmentend = ShootingParams.SegEndLambdas[csegment];				// the end point of the current segment
		act_off = 0; flex_off = 0;						// these offsets will
		if (RigidTip) act_off = 1;	else flex_off = 1;	// help with index calculations
		if ((csegment % 2) == (NUM_SEGMENTS % 2)) {		// we are at a rigid segment (the most proximal segment is assumed to be always flexible)
			mass += CathParams.ActMass[(csegment + act_off) / 2] * (cumpos - lastsegmentend) / (thissegmentend - lastsegmentend);  // we will add the partial mass of the actuator
		}
		mass += (cumpos - lastsegmentend) * CathParams.rho[csegment];			// we will add the partial mass of the flexible substrate
		for (int j = 0; j < (csegment + act_off) / 2; j++)		mass += CathParams.ActMass[j];			// add masses the distal actuators
		for (int j = 0; j < csegment; j++)		mass += CathParams.rho[j] * CathParams.SegLengths[j];	// add masses of the flexible substrate from distal segments
		mMult_sA<3, 1>(mass, CathConfig.g, ShootingParams.fcumlambda[i]);
	}

	mCopy_AB<3>(CathConfig.B0, ShootingParams.B0);
	ShootingParams.IntegrationStepSize = IntegrationStepSize;
	mCopy_AB<9>(CathConfig.R0, ShootingParams.R0);
	mCopy_AB<3>(CathConfig.p0, ShootingParams.p0);

	ShootingParams.ContactMode = ContactMode;
	mCopy_AB<3>(TipConstraintPoint, ShootingParams.TipConstraintPoint);
	mCopy_AB<3>(TipForce, ShootingParams.TipForce);
}
