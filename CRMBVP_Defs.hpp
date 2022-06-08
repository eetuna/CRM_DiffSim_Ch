#pragma once
#include <cmath>
#include "minpack.hpp"
#include "CRMDYN.hpp"

#define M_PI 3.14159265358979323846


template <typename adType>
void CRMShootingMethodBVP(	CRMShootingMethodParams<adType> in_Params,
                              double in_u0_initialguess[3], double in_ftip_initialguess[3],
                              adType out_u0[3], adType out_ftip[3], int& out_localmin) {

    /*Now we just assume no contact */
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
	for (int i = 0; i < NUM_STATES; i++) {
		if (i < 3) {
			x_0[i] = in_u0_initialguess[i];
		}
		else if (i < 12) {
			x_0[i] = in_Params.R0[i - 3];
		}
		else if (i < 15){
			x_0[i] = in_Params.p0[i - 12];
		}
        else if (i < 18){
            x_0[i] = in_Params.v0[i-15];
        }
        else {
            x_0[i] = in_Params.w0[i - 18];
        }
	}
	bool FinalValueOnly = true;
	NLEqnParams<adType> NLEParams;

	CRMSolverIVP_Prep(x_0, in_Params.IntegrationStepSize,
		in_Params.Li, in_Params.dlambdainv,
		in_Params.SegEndLambdas, in_Params.LocMarkerLambdas,
		in_Params.K, in_Params.Kinv, in_Params.ustar,
		in_Params.MagMoment, in_Params.fcumlambda,
		in_Params.B0, in_Params.g,
		in_Params.v_L_pre, in_Params.w_L_pre, in_Params.actMass, in_Params.actInertia,
		FinalValueOnly,
		NLEParams);
	NLEParams.ContactMode = ContactMode;
	mCopy_AB<3>(in_Params.TipConstraintPoint, NLEParams.TipConstraintPoint);

	// Scale parameters and call the nonlinear equation solver
	const double uscaleinv = 1.0 / IVALUE_SCALE_U;
	const double fscaleinv = 1.0 / IVALUE_SCALE_F;

    const double vscaleinv = 1.0 / IVALUE_SCALE_V;
    const double wscaleinv = 1.0 / IVALUE_SCALE_W;

	double* initialguessscaled = new double [NLEq_Dim];
	adType* returnedparamscaled = new adType [NLEq_Dim];
	if (ContactMode == ContactModeType::FREE_TIP) {
		for (int i = 0; i < 3; i++) {
			initialguessscaled[i] = uscaleinv * in_u0_initialguess[i];
//            initialguessscaled[i+3] = vscaleinv * in_v0_initialguess[i];
//            initialguessscaled[i+6] = wscaleinv * in_w0_initialguess[i];
            NLEParams.TipForce[i] = in_Params.TipForce[i];  // if the catheter is not in contact, the tip force specified within in_Params needs to be used; this would not be scaled as it is not changed by the solver
		}
	}
	else { // FIXED_TIP
		for (int i = 0; i < 3; i++) {
			initialguessscaled[i] = uscaleinv * in_u0_initialguess[i];
			initialguessscaled[i+3] = fscaleinv * in_ftip_initialguess[i];
		}

	}



//    for (int i = 0; i < NUM_ACT_SET; ++i) {
//        std::cout << "ActInertia: " << NLEParams.actInertia[i][0] << " " << NLEParams.actInertia[i][4] << " " << NLEParams.actInertia[i][8] << std::endl;
//        std::cout << "ActMass: " << NLEParams.actMass[i]  << std::endl;
//    }


    int localmin = 0, errorcode = 0;

	adType* x = new adType[NLEq_Dim]; // we will create a new variable here and not use initial guess scaled since truss-region-dogleg algorithm uses the same variable for both input and output
	adType* residual = new adType[3];

	int info;
	int lwa = (NLEq_Dim * (3 * NLEq_Dim + 13)) / 2; // what is this?
	double tol = 0.00001;
	adType* wa = new adType[lwa];
	for (int i = 0; i < NLEq_Dim; i++) x[i] = initialguessscaled[i];


#if defined( TRUSTREGION )
	TrustRegionDogleg(NLEq_Dim, x, residual, tol, info, wa, lwa, NLEParams);
#else // undefined
	exit(1);
#endif
	localmin = (info == 1) ? 0 : (info - 1);
	for (int i = 0; i < NLEq_Dim; i++) returnedparamscaled[i] = x[i];

//    for (int i = 0; i < RESIDUALDIM; ++i) {
//        std::cout << "residual " << residual[i]<< std::endl;
//    }

//    for (int i = 0; i < NLEq_Dim; ++i) {
//        std::cout << "returnedparamscaled " << returnedparamscaled[i]<< std::endl;
//    }


    delete[] residual;
    delete[] x;
    delete[] wa;

	if (ContactMode == ContactModeType::FREE_TIP) {
		for (int i = 0; i < 3; i++) {
			out_u0[i] = IVALUE_SCALE_U * returnedparamscaled[i];
//            out_v0[i] = IVALUE_SCALE_V * returnedparamscaled[i+3];
//            out_w0[i] = IVALUE_SCALE_W * returnedparamscaled[i+6];

            out_ftip[i] = in_Params.TipForce[i];  // if it is free-tip, return the tip force specified within in_Params
		}
	}
	else { // FIXED_TIP
		for (int i = 0; i < EQNDIMENSION; i++) {
			out_u0[i] = IVALUE_SCALE_U * returnedparamscaled[i];
			out_ftip[i] = IVALUE_SCALE_F * returnedparamscaled[i+3];
		}
	}

	out_localmin = localmin;
	delete[] initialguessscaled;
	delete[] returnedparamscaled;
}


//
//template <typename adType>
//void CRMShootingMethodBVP(
//	ContactModeType in_ContactMode,					// in_ContactMode == FREE_TIP if the catheter is not in contact with a surface, FIXED_TIP if catheter tip is constrained to be at TipContraintPoint
//	double 	in_u0_initialguess[3],					// Initial guess for the local curvature vector at the entry point
//    double  in_v0[3],
//    double  in_w0[3],
//	double 	in_ftip_initialguess[3],				// Initial guess for the tip force (\lambda=0), used when in_ContactMode == FIXED_TIP
//	double	in_InsertedLength, 						// Inserted Length (length of the catheter from the entry point to the tip)
//	double 	in_ActuationCurrents[NUM_ACT_SET][3], 	// Actuation current for each of the actuation coils
//													//    actuator units are numbered/ordered from the tip of the catheter towards the base (distal to proximal)
//	double	in_TipConstraintPoint[3],				// The spatial coordinates of the point where the catheter tip is constrained to be (used if in_ContactMode == FIXED_TIP)
//	double	in_TipForce[3],							// External point force (in spatial coordinates) applied at the tip of the catheter (\lambda = 0) - this value will not be used if in_ContactMode == FIXED_TIP
//	double	in_IntegrationStepSize,					// Stepsize used in numerical integration along the length of the catheter
//	// Catheter Configuration Parameters
//	double 	in_B0[3],								// B0 field vector of the MRI scanner (in spatial coordinates)
//	double 	in_g[3],								// Gravity vector (in spatial coordinates)
//	double 	in_p0[3],								// Catheter entry port position (in spatial coordinates)
//	double 	in_R0[9],								// Catheter orientation at the entry port (relative to the spatial frame); 3x3 matrix stored in row major order
//	// Catheter Model Parameters
//	//	For all parameters below, segments and actuator units are numbered/ordered from the tip of the catheter towards the base (distal to proximal)
//	double 	in_SegLengths[NUM_SEGMENTS],			// Array of segment lengths; NUM_SEGMENTS long array
//	double	in_LocMarkers[NUM_LOCALIZATION_MARKERS],// Array of localization marker locations (in lambda coordinates); NUM_LOCALIZATION_MARKERS long array
//	double	in_InnerRadius[NUM_FLEX_SEG],			// Inner radii of the flexible catheter segments
//	double	in_OuterRadius[NUM_FLEX_SEG],			// Outer radii of the flexible catheter segments
//	double	in_YoungsModulus[NUM_FLEX_SEG],			// Youngs Moduli of the flexible catheter segments
//	double	in_ShearModulus[NUM_FLEX_SEG],			// Shear Moduli of the flexible catheter segments
//	double 	in_ustar[NUM_FLEX_SEG][3],				// Local curvature in unloaded configuration for each of the flexible segments; (NUM_FLEX_SEG)*3 long array, (NUM_FLEX_SEG) 3x1 vectors
//													//		Actuator segments are assumed to be straight
//	double 	in_CoilAlignmentAngles[NUM_ACT_SET][2],	// Coil Alignment Angles; NUM_ACT_SET*2 long array, for each actuator set, the angle for the first coil is relative to x axis, and the angle for the second coil is relative to y axis
//	double 	in_CoilTurnAreaMat[NUM_ACT_SET][9], 	// Coil Turn Area matrices; NUM_ACT_SET*9 long array, NUM_ACT_SET 3x3 matrices stored in row major order
//	double 	in_rho[NUM_SEGMENTS],					// Length density (mass per unit length) of the flexible catheter substrate (tubing); (NUM_SEGMENTS) long array
//	double 	in_ActMass[NUM_ACT_SET],				// Actuator segment masses, does not include the flexible substrate; (NUM_ACT_SET) long array
//    double  in_ActInertia[NUM_ACT_SET][9],
//	double  in_v_L_pre[3],
//	double  in_w_L_pre[3],
//	// Outputs
//	double 	out_u0[3], 								// Local curvature vector at the entry point calculated through the solution of BVP
//    double 	out_ftip[3],							// Tip force calculated through the solution of BVP (\lambda=0) --- used when in_ContactMode == FIXED_TIP
//	int& out_localmin							// out_localmin!=0 if algorithms is stuck at a local minimum, or cannot make further progress
//) {
//
//	// Copy input parameters to local variables --- needed for dataflow optimization
//    CRMDynamicsModelParams CathParams;
//	CatheterConfiguration CathConfig;
//	CRMShootingMethodBVP_Prep(in_B0, in_g, in_p0, in_R0, in_v0, in_w0,
//		in_SegLengths, in_LocMarkers, in_InnerRadius, in_OuterRadius, in_YoungsModulus, in_ShearModulus,
//		in_ustar, in_CoilAlignmentAngles, in_CoilTurnAreaMat, in_rho, in_ActMass, in_ActInertia,
//		CathParams, CathConfig);
//
//	adType InsertedLength = in_InsertedLength;
//	adType ActuationCurrents[NUM_ACT_SET][3];
//	mCopy_AB<NUM_ACT_SET * 3>(&(in_ActuationCurrents[0][0]), &(ActuationCurrents[0][0]));
//	double IntegrationStepSize = in_IntegrationStepSize;
//	// Calculate Shooting Method Parameter Set from model and configuration parameters
//	CRMShootingMethodParams<adType> ShootingParams;
//	CRMConstructShootingMethodParamSet(CathParams, CathConfig,
//		InsertedLength, ActuationCurrents, in_ContactMode, in_TipConstraintPoint, in_TipForce,
//		IntegrationStepSize, in_v_L_pre, in_w_L_pre,
//		ShootingParams);
//
//	adType u0_calc[3], ftip_calc[3];
//
//	CRMShootingMethodBVP(ShootingParams, in_u0_initialguess, in_ftip_initialguess, u0_calc, ftip_calc, out_localmin);
//
//	for (int i = 0; i < 3; i++) out_u0[i] = dVal(u0_calc[i]);
//	for (int i = 0; i < 3; i++) out_ftip[i] = dVal(ftip_calc[i]);
//
//}


template <typename adType>
void NLEquation(adType in_x[], adType out_y[], NLEqnParams<adType> Params) {

	adType x_N[NUM_STATES];
	double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3];
	adType WrenchResidual[RESIDUALDIM];
    UHistory u_history[NUM_FLEX_SEG];

	adType u_0[3], v_0[3], w_0[3], ftip[3];
	// don't forget to scale parameters before passing to the CRMSolverIVP
	if (Params.ContactMode == ContactModeType::FREE_TIP) {
		for (int i = 0; i < 3; i++) {
			u_0[i] = IVALUE_SCALE_U * in_x[i];
//            v_0[i] = IVALUE_SCALE_V * in_x[i+3];
//            w_0[i] = IVALUE_SCALE_W * in_x[i+6];
			ftip[i] = Params.TipForce[i];  // for free-tip, this parameter is not given by the nonlinear equation solver, and hence, does not need to be scaled
		}
	}
	else { // FIXED_TIP
		for (int i = 0; i < 3; i++) {
			u_0[i] = IVALUE_SCALE_U * in_x[i];
			ftip[i] = IVALUE_SCALE_F * in_x[i + 9];
		}
	}
	// We will only call the IVP_Core, since preprocessing is already done
    CRMSolverIVP_Core(Params, u_0, ftip, x_N, WrenchResidual, p_atLocMarkers, u_history );

	// don't forget to scale parameters before returning to the nonlinear equation solver
	if (Params.ContactMode == ContactModeType::FREE_TIP) {
		for (int i = 0; i < 3; i++) {
            out_y[i] = RESIDUAL_SCALE_F * (WrenchResidual[i] * WrenchResidual[i])  + RESIDUAL_SCALE_M * (WrenchResidual[i+3] * WrenchResidual[i+3] ) ;
//            out_y[i+3] = RESIDUAL_SCALE_M * WrenchResidual[i+3];
		}
	}
	else { // FIXED_TIP
		for (int i = 0; i < RESIDUALDIM; i++) {
			out_y[i] = RESIDUAL_SCALE_M * WrenchResidual[i];
			out_y[i + RESIDUALDIM] = RESIDUAL_SCALE_P * (x_N[i + 12] - Params.TipConstraintPoint[i]);
		}
	}


}

void CRMShootingMethodBVP_Prep (
		double in_B0[3], double in_g[3], double in_p0[3], double in_R0[9], double in_v0[3], double in_w0[3],
		double in_SegLengths[NUM_SEGMENTS], double in_LocMarkers[NUM_LOCALIZATION_MARKERS],
		double in_InnerRadius[NUM_FLEX_SEG],	double in_OuterRadius[NUM_FLEX_SEG],
		double in_YoungsModulus[NUM_FLEX_SEG], double in_ShearModulus[NUM_FLEX_SEG],
		double in_ustar[NUM_FLEX_SEG][3],
		double in_CoilAlignmentAngles[NUM_ACT_SET][2], double in_CoilTurnAreaMat[NUM_ACT_SET][9],
		double in_rho[NUM_SEGMENTS],	double in_ActMass[NUM_ACT_SET], double in_ActInertia[NUM_ACT_SET][9],
        CRMDynamicsModelParams &CathParams, CatheterConfiguration &CathConfig ) {

	mCopy_AB<NUM_SEGMENTS>(in_SegLengths,CathParams.SegLengths);
	mCopy_AB<NUM_LOCALIZATION_MARKERS>(in_LocMarkers,CathParams.LocMarkers);
	mCopy_AB<NUM_FLEX_SEG>(in_InnerRadius,CathParams.InnerRadius);
	mCopy_AB<NUM_FLEX_SEG>(in_OuterRadius,CathParams.OuterRadius);
	mCopy_AB<NUM_FLEX_SEG>(in_YoungsModulus,CathParams.YoungsModulus);
	mCopy_AB<NUM_FLEX_SEG>(in_ShearModulus,CathParams.ShearModulus);
	mCopy_AB<NUM_FLEX_SEG,3>(in_ustar,CathParams.ustar);
	mCopy_AB<NUM_ACT_SET,2>(in_CoilAlignmentAngles,CathParams.CoilAlignmentAngles);
	mCopy_AB<NUM_ACT_SET,9>(in_CoilTurnAreaMat,CathParams.CoilTurnAreaMat);
	mCopy_AB<NUM_SEGMENTS>(in_rho,CathParams.rho);
	mCopy_AB<NUM_ACT_SET>(in_ActMass,CathParams.ActMass);
    mCopy_AB<NUM_ACT_SET, 9>(in_ActInertia,CathParams.ActInertia);

	mCopy_AB<3>(in_B0,CathConfig.B0);
	mCopy_AB<3>(in_g,CathConfig.g);
	mCopy_AB<3>(in_p0,CathConfig.p0);
	mCopy_AB<9>(in_R0,CathConfig.R0);
    mCopy_AB<3>(in_v0,CathConfig.v0);
    mCopy_AB<3>(in_w0,CathConfig.w0);
}

template <typename adType>
void CRMConstructShootingMethodParamSet(	CRMDynamicsModelParams CathParams, CatheterConfiguration CathConfig,
											adType InsertionLength, adType ActuationCurrents[NUM_ACT_SET][3],
											ContactModeType ContactMode,
											double TipConstraintPoint[3], double TipForce[3],
											double IntegrationStepSize,
                                            double in_v_L_pre[3], double in_w_L_pre[3],
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
	mCopy_AB<NUM_FLEX_SEG, 3>(CathParams.ustar, ShootingParams.ustar);
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

        ShootingParams.actMass[i] = CathParams.ActMass[i];
        for (int j = 0; j < 9; ++j) {
            ShootingParams.actInertia[i][j] = CathParams.ActInertia[i][j];
        }
	}

    mCopy_AB<3>(in_v_L_pre, ShootingParams.v_L_pre);
    mCopy_AB<3>(in_w_L_pre, ShootingParams.w_L_pre);

    mCopy_AB<3>(CathConfig.v0, ShootingParams.v0);
    mCopy_AB<3>(CathConfig.w0, ShootingParams.w0);


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
    mCopy_AB<3>(CathConfig.g, ShootingParams.g);
	ShootingParams.IntegrationStepSize = IntegrationStepSize;
	mCopy_AB<9>(CathConfig.R0, ShootingParams.R0);
	mCopy_AB<3>(CathConfig.p0, ShootingParams.p0);

	ShootingParams.ContactMode = ContactMode;
	mCopy_AB<3>(TipConstraintPoint, ShootingParams.TipConstraintPoint);
	mCopy_AB<3>(TipForce, ShootingParams.TipForce);
}
