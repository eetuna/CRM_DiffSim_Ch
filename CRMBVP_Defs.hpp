#pragma once
#include <cmath>
#include "Numerical_methods/minpack.hpp"

#define M_PI 3.14159265358979323846


template <typename adType>
void CRMShootingMethodBVP(	CRMShootingMethodParams<adType> in_Params, 
							const double in_u0_initialguess[NUM_FLEX_SEG][3],
							adType out_u0[NUM_FLEX_SEG][3], adType out_ftip[3], int& out_localmin) {

	ContactModeType ContactMode = in_Params.ContactMode;
	int NLEq_Dim;  // Dimension of the Nonlinear Equation to Solve
    NLEq_Dim = NUM_RESIDUAL;


	// Call CRMSolverIVP_Prep, to pre-process parameters
	adType x_0[NUM_STATES];
	for (int i = 0; i < NUM_STATES; i++) {
		if (i < 3) {
			x_0[i] = in_u0_initialguess[0][i];
		}
		else if (i < 12) {
			x_0[i] = in_Params.R0[i - 3];
		}
		else if (i < 15) {
			x_0[i] = in_Params.p0[i - 12];
		}
	}
	bool FinalValueOnly = true;
	NLEqnParams<adType> NLEParams;
    double in_mL[NUM_ACT_SET][3] ={0.0,0.0,0.0};
    double in_nL[NUM_ACT_SET][3] ={0.0,0.0,0.0};

    CRMSolverIVP_Prep(x_0, in_Params.IntegrationStepSize,
		in_Params.Li, in_Params.dlambdainv,
		in_Params.SegEndLambdas, in_Params.LocMarkerLambdas,
		in_Params.K, in_Params.Kinv, in_Params.ustar,
		in_Params.MagMoment, in_Params.fcumlambda,
		in_Params.B0, in_Params.g, in_Params.actMass, in_Params.actInertia, in_Params.damping, in_Params.DELTA_T,
                      in_Params.v_L_pre, in_Params.w_L_pre, in_Params.p_pre, in_Params.R_pre,
                      in_mL, in_nL, false, FinalValueOnly,NLEParams);

    NLEParams.ContactMode = ContactMode;
	mCopy_AB<3>(in_Params.TipConstraintPoint, NLEParams.TipConstraintPoint);

	// Scale parameters and call the nonlinear equation solver
	const double uscaleinv = 1.0 / IVALUE_SCALE_U;

//    const double fscaleinv = 1.0 / IVALUE_SCALE_F;
	auto* initialguessscaled = new double [NLEq_Dim];
	auto* returnedparamscaled = new adType [NLEq_Dim];

    if (ContactMode == ContactModeType::FREE_TIP) {
        for (int i = 0; i < NUM_FLEX_SEG; ++i) {
            for (int j = 0; j < 3; ++j) {
                initialguessscaled[j+3*i] = uscaleinv * in_u0_initialguess[i][j];
            }
        }
		for (int i = 0; i < 3; i++) {
            NLEParams.TipForce[i] = in_Params.TipForce[i];  // if the catheter is not in contact, the tip force specified within in_Params needs to be used; this would not be scaled as it is not changed by the solver
		}
	}

    int localmin = 0;// errorcode = 0;

	auto* x = new adType[NLEq_Dim]; // we will create a new variable here and not use initial guess scaled since truss-region-dogleg algorithm uses the same variable for both input and output

	auto* residual = new adType[NLEq_Dim];
	int info;
	int lwa = (NLEq_Dim * (3 * NLEq_Dim + 13)) / 2;
	double tol = 0.00001;
	auto* wa = new adType [lwa];
	for (int i = 0; i < NLEq_Dim; i++) x[i] = initialguessscaled[i];
#if defined( TRUSTREGION )
	TrustRegionDogleg(NLEq_Dim, x, residual, tol, info, wa, lwa, NLEParams);
#else // undefined
	exit(1);
#endif
	localmin = (info == 1) ? 0 : (info - 1);
	for (int i = 0; i < NLEq_Dim; i++) returnedparamscaled[i] = x[i];
	delete[] residual;
	delete[] x;
	delete[] wa;

     for (int i = 0; i < NUM_FLEX_SEG; ++i) {
            for (int j = 0; j < 3; ++j) {
                out_u0[i][j] = IVALUE_SCALE_U * returnedparamscaled[j+3*i];
            }
        }
		for (int i = 0; i < 3; i++) {
            out_ftip[i] = in_Params.TipForce[i];  // if it is free-tip, return the tip force specified within in_Params
		}

	out_localmin = localmin;
	delete[] initialguessscaled;
    delete[] returnedparamscaled;

}


template <typename adType>
void NLEquation(adType in_x[], adType out_y[], NLEqnParams<adType> Params) {

//    adType x_N[NUM_STATES];
//    double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3];
//    adType OutResidual[NUM_RESIDUAL];
    auto* x_N = new adType [NUM_STATES];
    auto* p_atLocMarkers = new adType [NUM_LOCALIZATION_MARKERS][3];
    auto* OutResidual = new adType [NUM_RESIDUAL];

    // output for time advance, not used in BVP, just placeholders
    adType u_0[NUM_FLEX_SEG][3], ftip[3];
    // don't forget to scale parameters before passing to the CRMSolverIVP
    if (Params.ContactMode == ContactModeType::FREE_TIP) {
        for (int i = 0; i < NUM_FLEX_SEG; ++i) {
            for (int j = 0; j < 3; ++j) {
                u_0[i][j] = IVALUE_SCALE_U * in_x[j+3*i];
            }
        }
        for (int i = 0; i < 3; i++) {
            ftip[i] = Params.TipForce[i];  // for free-tip, this parameter is not given by the nonlinear equation solver, and hence, does not need to be scaled
        }
    }


    adType pcoil[NUM_ACT_SET][3], Rcoil[NUM_ACT_SET][9];


    // We will only call the IVP_Core, since preprocessing is already done
    CRMSolverIVP_Core(Params, u_0, ftip, x_N, OutResidual, pcoil, Rcoil, p_atLocMarkers );

    // don't forget to scale parameters before returning to the nonlinear equation solver
    if (Params.ContactMode == ContactModeType::FREE_TIP) {
        for (int i = 0; i < NUM_RESIDUAL; i++) {
            out_y[i] = RESIDUAL_SCALE_M * OutResidual[i];// RESIDUAL_SCALE_F * WrenchResidual[i] ;
        }
    }

    delete[] x_N;
    delete[] p_atLocMarkers;
    delete[] OutResidual;

}


void CRMShootingMethodBVP_Prep (
        double in_B0[3], double in_g[3], double in_p0[3], double in_R0[9],
        double in_SegLengths[NUM_SEGMENTS], double in_LocMarkers[NUM_LOCALIZATION_MARKERS],
        double in_InnerRadius[NUM_FLEX_SEG],	double in_OuterRadius[NUM_FLEX_SEG],
        double in_YoungsModulus[NUM_FLEX_SEG], double in_ShearModulus[NUM_FLEX_SEG],
        double in_ustar[NUM_FLEX_SEG][3],
        double in_CoilAlignmentAngles[NUM_ACT_SET][2], double in_CoilTurnAreaMat[NUM_ACT_SET][9],
        double in_rho[NUM_SEGMENTS],	double in_ActMass[NUM_ACT_SET], double in_ActInertia[NUM_ACT_SET][9],
        CRMCatheterModelParams &CathParams, CatheterConfiguration &CathConfig ) {

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
}

template <typename adType>
void CRMConstructShootingMethodParamSet(	CRMCatheterModelParams CathParams, CatheterConfiguration CathConfig,
                                            adType InsertionLength, adType ActuationCurrents[NUM_ACT_SET][3],
                                            ContactModeType ContactMode, double TipConstraintPoint[3], double TipForce[3], double IntegrationStepSize,
                                            double in_v_L_pre[NUM_ACT_SET][3], double in_w_L_pre[NUM_ACT_SET][3], double in_p_pre[NUM_ACT_SET][3],
                                            double in_R_pre[NUM_ACT_SET][9], double in_damping[NUM_ACT_SET][6], double in_DELTA_T,
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

//    mCopy_AB<3>(in_v_L_pre, ShootingParams.v_L_pre);
//    mCopy_AB<3>(in_w_L_pre, ShootingParams.w_L_pre);
//    mCopy_AB<3>(in_p_pre, ShootingParams.p_pre);
//    mCopy_AB<9>(in_R_pre, ShootingParams.R_pre);
//    mCopy_AB<6>(in_damping, ShootingParams.damping);

    for (int i = 0; i < NUM_ACT_SET; ++i) {
        for (int j = 0; j < 3; ++j) {
            ShootingParams.v_L_pre[i][j] = in_v_L_pre[i][j];
            ShootingParams.w_L_pre[i][j] = in_w_L_pre[i][j];
            ShootingParams.p_pre[i][j] = in_p_pre[i][j];
        }
        for (int j = 0; j < 9; ++j) {
            ShootingParams.R_pre[i][j] = in_R_pre[i][j];
        }
        for (int j = 0; j < 6; ++j) {
            ShootingParams.damping[i][j] = in_damping[i][j];
        }
    }

    ShootingParams.DELTA_T = in_DELTA_T;

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
        for (int j = 0; j < (csegment + act_off) / 2; j++)
            mass += CathParams.ActMass[j];			// add masses the distal actuators
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
