#pragma once
#include <cmath>
#include "minpack.hpp"

#define M_PI 3.14159265358979323846


template <typename adType>
void CRMShootingMethodBVP(	CRMShootingMethodParams<adType> in_Params,  adType in_initial_guess[NUM_RESIDUAL], double in_ftip_initialguess[3],
							adType out_calc[NUM_RESIDUAL], adType out_ftip[3], int& out_localmin) {

    double in_u0_initialguess[3], in_n0_initialguess[3], in_nL_initialguess[3], in_uL_initialguess[3];
    for (int i = 0; i < 3; ++i) {
        in_u0_initialguess[i] = in_initial_guess[i];
        in_n0_initialguess[i] = in_initial_guess[i+3];
    }

    if (NUM_RESIDUAL > 6)
    {
        for (int i = 0; i < 3; ++i) {
            in_uL_initialguess[i] = in_initial_guess[i+6];
            in_nL_initialguess[i] = in_initial_guess[i+9];
        }
    }

    ContactModeType ContactMode = in_Params.ContactMode;
	int NLEq_Dim;  // Dimension of the Nonlinear Equation to Solve
	if (ContactMode == ContactModeType::FREE_TIP) {
		NLEq_Dim = NUM_RESIDUAL;
	}
	else { // FIXED_TIP
		NLEq_Dim = 9;
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
		else if (i < 15) {
			x_0[i] = in_Params.p0[i - 12];
		}
        else if (i < 18) {
            x_0[i] = in_Params.v0[i - 15];
        }
        else if (i < 21) {
            x_0[i] = in_Params.w0[i - 18];
        }
        else if (i < 24) {
            x_0[i] = in_n0_initialguess[i - 21];
        }
        else{
            x_0[i] = 0.0; //else, m_0 initialize with 0, initialize in IVPcore
        }
	}
	bool FinalValueOnly = true;
	NLEqnParams<adType> NLEParams;


    CRMSolverIVP_Prep(x_0, in_Params.IntegrationStepSize,
		in_Params.Li, in_Params.dlambdainv,
		in_Params.SegEndLambdas, in_Params.LocMarkerLambdas,
		in_Params.K, in_Params.Kinv, in_Params.ustar,
		in_Params.MagMoment, in_Params.fcumlambda,
		in_Params.B0, in_Params.g, in_Params.u_history, in_Params.v_history, in_Params.w_history,
        in_Params.v_L_pre, in_Params.w_L_pre, in_Params.actMass, in_Params.actInertia,
		in_Params.rho, in_Params.tubingInertia,
        in_Params.DELTA_T, in_Params.damping_tubing_mat, in_Params.damping_coil_diag,
                      in_Params.inv_K_D, FinalValueOnly, NLEParams);

    NLEParams.ContactMode = ContactMode;
	mCopy_AB<3>(in_Params.TipConstraintPoint, NLEParams.TipConstraintPoint);

	// Scale parameters and call the nonlinear equation solver
	const double uscaleinv = 1.0 / IVALUE_SCALE_U;
    const double nscaleinv = 1.0 / IVALUE_SCALE_N;

    const double fscaleinv = 1.0 / IVALUE_SCALE_F;
	auto* initialguessscaled = new double [NLEq_Dim];
	auto* returnedparamscaled = new adType [NLEq_Dim];

    if (ContactMode == ContactModeType::FREE_TIP) {

		for (int i = 0; i < 3; i++) {
			initialguessscaled[i] = uscaleinv * in_u0_initialguess[i];
            initialguessscaled[i+3] = nscaleinv * in_n0_initialguess[i];
            if (NUM_RESIDUAL > 6){
                initialguessscaled[i+6] = uscaleinv * in_uL_initialguess[i];
                initialguessscaled[i+9] = nscaleinv * in_nL_initialguess[i];
            }
            NLEParams.TipForce[i] = in_Params.TipForce[i];  // if the catheter is not in contact, the tip force specified within in_Params needs to be used; this would not be scaled as it is not changed by the solver
		}
	}

    int localmin = 0, errorcode = 0;

	adType* x = new adType[NLEq_Dim]; // we will create a new variable here and not use initial guess scaled since truss-region-dogleg algorithm uses the same variable for both input and output

	adType* residual = new adType[NLEq_Dim];
	int info;
	int lwa = (NLEq_Dim * (3 * NLEq_Dim + 13)) / 2;
	double tol = 0.00001;
	adType* wa = new adType [lwa];
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


	/*
	// AUTODIFF TEST CODE
	using Eigen::MatrixXd;
	VectorXadType u(NLEq_Dim), y(NLEq_Dim);
	for (int i=0;i<NLEq_Dim;i++) u(i)=x[i];
	MatrixXd J = jacobian(NLEquationAD, wrt(u), at(u, NLEParams), y);
	std::cout << "y = \n" << y << std::endl;    // print the evaluated output vector F
	std::cout << "J = \n" << J << std::endl;    // print the evaluated Jacobian matrix dF/dx
	*/

	if (ContactMode == ContactModeType::FREE_TIP) {
		for (int i = 0; i < 3; i++) {
			out_calc[i] = IVALUE_SCALE_U * returnedparamscaled[i];
            out_calc[i+3] = IVALUE_SCALE_N * returnedparamscaled[i+3];
            if (NUM_RESIDUAL > 6){
                out_calc[i+6] = IVALUE_SCALE_U * returnedparamscaled[i+6];
                out_calc[i+9] = IVALUE_SCALE_N * returnedparamscaled[i+9];
            }
            out_ftip[i] = in_Params.TipForce[i];  // if it is free-tip, return the tip force specified within in_Params
		}
	}
	else { // FIXED_TIP
		for (int i = 0; i < 3; i++) {
            out_calc[i] = IVALUE_SCALE_U * returnedparamscaled[i];
			out_ftip[i] = IVALUE_SCALE_F * returnedparamscaled[i+3];
		}
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
//    paramHistory u_history[NUM_FLEX_SEG], v_history[NUM_FLEX_SEG], w_history[NUM_FLEX_SEG];
//    adType v_L[3], w_L[3], pL[3], RL[9];
//    adType ftip[3], in_initials[NUM_RESIDUAL];

    auto* x_N = new adType [NUM_STATES];
    auto* p_atLocMarkers = new adType [NUM_LOCALIZATION_MARKERS][3];
    auto* OutResidual = new adType [NUM_RESIDUAL];
    auto* u_history = new paramHistory [NUM_FLEX_SEG];
    auto* v_history = new paramHistory [NUM_FLEX_SEG];
    auto* w_history = new paramHistory [NUM_FLEX_SEG];

    auto* v_L = new adType [3];
    auto* w_L = new adType [3];
    auto* pL = new adType [3];
    auto* RL = new adType [9];
    auto* ftip = new adType [3];
    auto* in_initials = new adType [NUM_RESIDUAL];


    // don't forget to scale parameters before passing to the CRMSolverIVP
    if (Params.ContactMode == ContactModeType::FREE_TIP) {
        for (int i = 0; i < 3; i++) {
            in_initials[i] = IVALUE_SCALE_U * in_x[i];
            in_initials[i+3] = IVALUE_SCALE_N * in_x[i+3];
            if (NUM_RESIDUAL > 6){
                in_initials[i+6] = IVALUE_SCALE_U * in_x[i+6];
                in_initials[i+9] = IVALUE_SCALE_N * in_x[i+9];
            }

            ftip[i] = Params.TipForce[i];  // for free-tip, this parameter is not given by the nonlinear equation solver, and hence, does not need to be scaled
        }
    }
    else { // FIXED_TIP
        for (int i = 0; i < 3; i++) {
            in_initials[i] = IVALUE_SCALE_U * in_x[i];
            ftip[i] = IVALUE_SCALE_F * in_x[i + 9];
        }
    }

    // We will only call the IVP_Core, since preprocessing is already done
    CRMSolverIVP_Core(Params, in_initials, ftip, x_N, OutResidual, p_atLocMarkers,
                      u_history, v_history, w_history, v_L, w_L , pL, RL);

    delete[] u_history;
    delete[] v_history;
    delete[] w_history;
    delete[] x_N;
    delete[] p_atLocMarkers;
    delete[] v_L;
    delete[] w_L;
    delete[] pL;
    delete[] RL;
    delete[] ftip;
    delete[] in_initials;

    // don't forget to scale parameters before returning to the nonlinear equation solver
    if (Params.ContactMode == ContactModeType::FREE_TIP) {
        for (int i = 0; i < 3; i++) {
            out_y[i] = RESIDUAL_SCALE_M * OutResidual[i+3];// RESIDUAL_SCALE_F * WrenchResidual[i] ;
            out_y[i+3] = RESIDUAL_SCALE_F * OutResidual[i] ;
            if (NUM_RESIDUAL > 6){
                out_y[i+6] = RESIDUAL_SCALE_F * OutResidual[i+6];
                out_y[i+9] = RESIDUAL_SCALE_M * OutResidual[i+9];
            }
        }
    }

    delete[] OutResidual;

}


void CRMShootingMethodBVP_Prep (
        double in_B0[3], double in_g[3], double in_p0[3], double in_R0[9], double in_v0[3], double in_w0[3],
        double in_SegLengths[NUM_SEGMENTS], double in_LocMarkers[NUM_LOCALIZATION_MARKERS],
        double in_InnerRadius[NUM_FLEX_SEG],	double in_OuterRadius[NUM_FLEX_SEG],
        double in_YoungsModulus[NUM_FLEX_SEG], double in_ShearModulus[NUM_FLEX_SEG],
        double in_ustar[NUM_FLEX_SEG][3],
        double in_CoilAlignmentAngles[NUM_ACT_SET][2], double in_CoilTurnAreaMat[NUM_ACT_SET][9],
        double in_rho[NUM_SEGMENTS],	double in_ActMass[NUM_ACT_SET], double in_ActInertia[NUM_ACT_SET][9],
        double in_tubingInertia[NUM_ACT_SET][9],
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
    mCopy_AB<NUM_FLEX_SEG, 9>(in_tubingInertia,CathParams.tubingInertia);

    mCopy_AB<3>(in_B0,CathConfig.B0);
    mCopy_AB<3>(in_g,CathConfig.g);
    mCopy_AB<3>(in_p0,CathConfig.p0);
    mCopy_AB<9>(in_R0,CathConfig.R0);
    mCopy_AB<3>(in_v0,CathConfig.v0);
    mCopy_AB<3>(in_w0,CathConfig.w0);
}

template <typename adType>
void CRMConstructShootingMethodParamSet(	CRMCatheterModelParams CathParams, CatheterConfiguration CathConfig,
                                            adType InsertionLength, adType ActuationCurrents[NUM_ACT_SET][3],
                                            ContactModeType ContactMode,
                                            double TipConstraintPoint[3], double TipForce[3],
                                            double IntegrationStepSize,
                                            double in_v_L_pre[3], double in_w_L_pre[3],
                                            paramHistory u_history[NUM_FLEX_SEG],paramHistory v_history[NUM_FLEX_SEG],paramHistory w_history[NUM_FLEX_SEG],
                                            double in_DELTA_T, double in_damping_tubing[3], double in_damping_coil[6],
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

    mCopy_AB<NUM_SEGMENTS>(CathParams.rho, ShootingParams.rho);
    for (int i = 0; i < NUM_FLEX_SEG; ++i)
    {
        for (int j = 0; j < 9; ++j) {ShootingParams.tubingInertia[i][j] = CathParams.tubingInertia[i][j];}
    }

    //parameter history
    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
        int length_h = u_history[i].length;
        if (length < NUM_HISTORY_LENGTH) {
            ShootingParams.u_history[i].length = length_h;
            for (int j = 0; j < length_h; ++j) {
                ShootingParams.u_history[i].data[j] = u_history[i].data[j];
            }
            ShootingParams.v_history[i].length = length_h;
            for (int j = 0; j < length_h; ++j) {
                ShootingParams.v_history[i].data[j] = v_history[i].data[j];
            }
            ShootingParams.w_history[i].length = length_h;
            for (int j = 0; j < length_h; ++j) {
                ShootingParams.w_history[i].data[j] = w_history[i].data[j];
            }
        }else{std::cout <<" Maximum history length exceeded!";
            exit(1);
        }
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
    mCopy_AB<3>(CathConfig.g, ShootingParams.g);
    ShootingParams.IntegrationStepSize = IntegrationStepSize;
    mCopy_AB<9>(CathConfig.R0, ShootingParams.R0);
    mCopy_AB<3>(CathConfig.p0, ShootingParams.p0);

//    std:: cout << "damping_tubing " << in_damping_tubing[0] << " " << in_damping_tubing[1] << " " << in_damping_tubing[2] << std::endl;

    double damping_tubing[9], inv_K_D[9]; //(K+D)^-1
    damping_tubing[0] = in_damping_tubing[0];		damping_tubing[1] = 0.0;			damping_tubing[2] = 0.0;
    damping_tubing[3] = 0.0;			damping_tubing[4] = in_damping_tubing[1];		damping_tubing[5] = 0.0;
    damping_tubing[6] = 0.0;			damping_tubing[7] = 0.0;			damping_tubing[8] = in_damping_tubing[2];
    mCopy_AB<9>(damping_tubing, ShootingParams.damping_tubing_mat);

    double C_0 = 1 / in_DELTA_T;
    inv_K_D[0] = 1.0 / (K[0] + C_0 * in_damping_tubing[0]);		inv_K_D[1] = 0.0;				inv_K_D[2] = 0.0;
    inv_K_D[3] = 0.0;				inv_K_D[4] = 1.0 / (K[4] + C_0 * in_damping_tubing[1]);		inv_K_D[5] = 0.0;
    inv_K_D[6] = 0.0;				inv_K_D[7] = 0.0;				inv_K_D[8] = 1.0 / (K[8] + C_0 * in_damping_tubing[2]);
    mCopy_AB<9>(inv_K_D, ShootingParams.inv_K_D);
//    std:: cout << "inv_K_D " << inv_K_D[0] << " " << inv_K_D[4] << " " << inv_K_D[8] << std::endl;

    mCopy_AB<6>(in_damping_coil, ShootingParams.damping_coil_diag);

    ShootingParams.DELTA_T = in_DELTA_T;

    ShootingParams.ContactMode = ContactMode;
    mCopy_AB<3>(TipConstraintPoint, ShootingParams.TipConstraintPoint);
    mCopy_AB<3>(TipForce, ShootingParams.TipForce);
}
