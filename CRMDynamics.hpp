#pragma once
#include "CRMBVP_Defs.hpp"
#include "CRMIVP_Defs.hpp"

template <typename adType>
void CRM_Dynamics(const double in_x_t[NUM_STATES], const double control_t[NUM_CONTROL], CRMShootingMethodParams<double> BVPParams,  paramHistory in_u_history[NUM_FLEX_SEG],
                  paramHistory in_v_history[NUM_FLEX_SEG], paramHistory in_w_history[NUM_FLEX_SEG],
                  double in_v_L_pre[3], double in_w_L_pre[3], double in_h0_pre[NUM_FLEX_SEG], CRMDynamicsParams<adType> in_Params,
                  adType out_x_t[NUM_STATES], adType out_calc[NUM_RESIDUAL], adType out_ReportedMarkerPos[NUM_LOCALIZATION_MARKERS][3],
                  paramHistory out_u_history[NUM_FLEX_SEG], paramHistory out_v_history[NUM_FLEX_SEG], paramHistory out_w_history[NUM_FLEX_SEG],
                  double out_v_L[3], double out_w_L[3], double out_h0[NUM_FLEX_SEG], adType residual[NUM_RESIDUAL]){


    adType ActuationCurrents[NUM_ACT_SET][3];
    adType InsertedLength;
    adType ftip_calc[3];
    adType xf[NUM_STATES];
    int localmin; //// not using it right now

    for (int i = 0; i < NUM_ACT_SET; i++)	for (int j = 0; j < 3; j++)	ActuationCurrents[i][j] = control_t[i * 3 + j];
    InsertedLength = control_t[NUM_CONTROL -1 ];

    for (int i = 0; i < NUM_ACT_SET; ++i) {
        std::cout << "ActuationCurrents: " << ActuationCurrents[i][0] << " " << ActuationCurrents[i][1] << " " << ActuationCurrents[i][2] << std::endl;
    }

    for (int i = 0; i < NUM_RESIDUAL; ++i) {
        std::cout << "new initial guess: " << in_Params.initial_guess[i] << std::endl;
    }

    adType p0[3], R0[9], v0[3], w0[3];
    for (int i = 0; i < 3; ++i) {
        p0[i] = in_x_t[i+3+9];
        v0[i] = in_x_t[i+3+9+3];
        w0[i] = in_x_t[i+3+9+3+3];
    }
    for (int i = 0; i < 9; ++i) {
        R0[i] = in_x_t[i+3];
    }

    CRMDYNCoreDataUpdate(*(in_Params.CathParams), p0, R0, v0, w0, InsertedLength, ActuationCurrents,
                         in_Params.ContactMode,in_Params.TipConstraintPoint, in_Params.TipForce, in_Params.IntegrationStepSize,
                         in_u_history, in_v_history, in_w_history, in_v_L_pre, in_w_L_pre, in_h0_pre, BVPParams );

    /**
     * Shooting method calculating u0 and n0
     */
    // Cosserat Rod Model - Solve the Boundary Value Problem to calculate the equilibrium configuration of the catheter
    CRMShootingMethodBVP(BVPParams, in_Params.initial_guess,in_Params.ftip_initialguess, out_calc, ftip_calc, localmin);


    /**
     * Get the final results
     */
    CRMSolverIVP(BVPParams, out_calc, ftip_calc, in_Params.FinalValueOnly, xf, residual, out_ReportedMarkerPos, out_u_history, out_v_history, out_w_history,
                 out_v_L, out_w_L, out_h0);


    for (int i = 0; i < NUM_STATES; ++i) {
        out_x_t[i] = xf[i];
    }

}


template <typename adType>
void CRMDynamicParamPrep( CRMCatheterModelParams CathParams, CatheterConfiguration CathConfig,
                          adType InsertionLength, adType ActuationCurrents[NUM_ACT_SET][3],
                          ContactModeType ContactMode,
                          double TipConstraintPoint[3], double TipForce[3],
                          double IntegrationStepSize, paramHistory u_history[NUM_FLEX_SEG],
                          double u0_initialguess[3], double n0_initialguess[3], double ftip_initialguess[3],
                          double in_v_L_pre[3], double in_w_L_pre[3],
                          CRMShootingMethodParams<adType> &BVPParams, CRMDynamicsParams<adType> &DynamicsParams){


    DynamicsParams.CathParams = &CathParams;
    DynamicsParams.ContactMode = ContactMode;
    DynamicsParams.FinalValueOnly = false;  // we want the localization coil locations, too

	for (int i = 0; i < 3; i++) DynamicsParams.TipConstraintPoint[i] = TipConstraintPoint[i];
	for (int i = 0; i < 3; i++) DynamicsParams.TipForce[i] = TipForce[i];
	for (int i = 0; i < 3; i++) DynamicsParams.u0_initialguess[i] = u0_initialguess[i];
    for (int i = 0; i < 3; i++) DynamicsParams.n0_initialguess[i] = n0_initialguess[i];
    for (int i = 0; i < 3; i++) DynamicsParams.ftip_initialguess[i] = ftip_initialguess[i];
    DynamicsParams.IntegrationStepSize = IntegrationStepSize;

    CRMConstructShootingMethodParamSet(CathParams, CathConfig, InsertionLength, ActuationCurrents, ContactMode,
                                       TipConstraintPoint, TipForce, IntegrationStepSize, u_history, in_v_L_pre, in_w_L_pre, BVPParams);

    DynamicsParams.ShootingParams = &BVPParams;

};

/**
 * Shooting params that stays the same during time advance:
 * ShootingParams.SegEndLambdas, ShootingParams.dlambdainv, ShootingParams.LocMarkerLambdas,
 * ShootingParams.K,  ShootingParams.Kinv, ShootingParams.ustar,
 * ShootingParams.actMass, ShootingParams.actInertia,
 * ShootingParams.rho, ShootingParams.Area, ShootingParams.B0, ShootingParams.g, ShootingParams.fcumlambda[i]
 */
template <typename adType>
void CRMDYNCoreDataUpdate(CRMCatheterModelParams CathParams, double in_p0[3], double in_R0[9], double in_v0[3], double in_w0[3],
                        adType InsertionLength, adType ActuationCurrents[NUM_ACT_SET][3],
                        ContactModeType ContactMode,
                        double TipConstraintPoint[3], double TipForce[3],
                        double IntegrationStepSize, paramHistory in_u_history[NUM_FLEX_SEG],
                        paramHistory in_v_history[NUM_FLEX_SEG], paramHistory in_w_history[NUM_FLEX_SEG],
                        double in_v_L_pre[3], double in_w_L_pre[3], double in_h0_pre[NUM_FLEX_SEG],
                        CRMShootingMethodParams<adType> &ShootingParams){

    double length = ShootingParams.SegEndLambdas[NUM_SEGMENTS -1];

    ShootingParams.Li = MIN(InsertionLength, length);

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
    }
    ShootingParams.IntegrationStepSize = IntegrationStepSize;

    mCopy_AB<9>(in_R0, ShootingParams.R0);
    mCopy_AB<3>(in_p0, ShootingParams.p0);

    mCopy_AB<3>(in_v_L_pre, ShootingParams.v_L_pre);
    mCopy_AB<3>(in_w_L_pre, ShootingParams.w_L_pre);

    mCopy_AB<3>(in_v0, ShootingParams.v0);
    mCopy_AB<3>(in_w0, ShootingParams.w0);

    // history update
    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
        int length_u =  in_u_history[i].length;
        ShootingParams.u_history[i].length = length_u;
        ShootingParams.u_history[i].data = new double[length_u];
        for (int j = 0; j < length_u; ++j) {
            ShootingParams.u_history[i].data[j] = in_u_history[i].data[j];
        }

        ShootingParams.v_history[i].length = length_u;
        ShootingParams.v_history[i].data = new double[length_u];
        for (int j = 0; j < length_u; ++j) {
            ShootingParams.v_history[i].data[j] = in_v_history[i].data[j];
        }

        ShootingParams.w_history[i].length = length_u;
        ShootingParams.w_history[i].data = new double[length_u];
        for (int j = 0; j < length_u; ++j) {
            ShootingParams.w_history[i].data[j] = in_w_history[i].data[j];
        }

        ShootingParams.h0_pre[i] = in_h0_pre[i];
    }

    ShootingParams.ContactMode = ContactMode;
    mCopy_AB<3>(TipConstraintPoint, ShootingParams.TipConstraintPoint);
    mCopy_AB<3>(TipForce, ShootingParams.TipForce);

}
