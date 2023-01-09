#pragma once
#include "CRMBVP_Defs.hpp"
#include "CRMIVP_Defs.hpp"

template <typename adType>
void CRM_Dynamics(const double in_x_t[NUM_STATES], const double control_t[NUM_CONTROL], CRMShootingMethodParams<double> BVPParams,
                  double in_v_L_pre[3], double in_w_L_pre[3],  double pL_pre[3], double RL_pre[9], CRMDynamicsParams<adType> in_Params,
                  adType out_x_t[NUM_STATES], adType u0_calc[3], adType nL_calc[3], adType mL_calc[3], adType out_ReportedMarkerPos[NUM_LOCALIZATION_MARKERS][3],
                  double out_v_L[3], double out_w_L[3],	double out_pL_pre[3], double out_RL_pre[9]){

	adType ActuationCurrents[NUM_ACT_SET][3];
	adType InsertedLength;
	adType ftip_calc[3];
	adType xf[NUM_STATES];
    int localmin; //// not using it right now

	for (int i = 0; i < NUM_ACT_SET; i++)	for (int j = 0; j < 3; j++)	ActuationCurrents[i][j] = control_t[i * 3 + j];
	InsertedLength = control_t[NUM_CONTROL -1 ];

    adType p0[3], R0[9];
    for (int i = 0; i < 3; ++i) {
        p0[i] = in_x_t[i+3+9];
    }
    for (int i = 0; i < 9; ++i) {
        R0[i] = in_x_t[i+3];
    }
    CRMDYNCoreDataUpdate(*(in_Params.CathParams), p0, R0, InsertedLength, ActuationCurrents,
                         in_Params.ContactMode,in_Params.TipConstraintPoint, in_Params.TipForce,
                         in_Params.IntegrationStepSize, in_v_L_pre, in_w_L_pre, pL_pre, RL_pre, BVPParams );

    DynamicsBVP(BVPParams, in_Params.u0_initialguess, in_Params.mL_initialguess,
                in_Params.nL_initialguess, in_Params.ftip_initialguess,
                u0_calc, mL_calc, nL_calc, ftip_calc, localmin);

    double x_coil[NUM_COIL_STATES];
    DYNSolverIVP(BVPParams, u0_calc, mL_calc, nL_calc, ftip_calc,
                 true, xf, x_coil,out_ReportedMarkerPos);

    for (int i = 0; i < 3; ++i) {
        out_v_L[i] = x_coil[i];
        out_w_L[i] = x_coil[i+3];
    }

    for (int i = 0; i < 3; ++i) {
        out_pL_pre[i] = x_coil[i+6];
    }

    for (int i = 0; i < 9; ++i) {
        out_RL_pre[i] = x_coil[i+9];
    }

    for (int i = 0; i < NUM_STATES; ++i) {
        out_x_t[i] = xf[i];
    }



}

/**
 * Shooting params that stays the same during time advance:
 * ShootingParams.SegEndLambdas, ShootingParams.dlambdainv, ShootingParams.LocMarkerLambdas,
 * ShootingParams.K,  ShootingParams.Kinv, ShootingParams.ustar,
 * ShootingParams.actMass, ShootingParams.actInertia,
 * ShootingParams.rho, ShootingParams.Area, ShootingParams.B0, ShootingParams.g, ShootingParams.fcumlambda[i]
 */
template <typename adType>
void CRMDYNCoreDataUpdate(CRMCatheterModelParams CathParams, double in_p0[3], double in_R0[9],
                        adType InsertionLength, adType ActuationCurrents[NUM_ACT_SET][3],
                        ContactModeType ContactMode, double TipConstraintPoint[3], double TipForce[3], double IntegrationStepSize,
                        double in_v_L_pre[3], double in_w_L_pre[3], double pL_pre[3], double RL_pre[9],
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

    mCopy_AB<3>(pL_pre, ShootingParams.p_pre);
    mCopy_AB<9>(RL_pre, ShootingParams.R_pre);

    ShootingParams.ContactMode = ContactMode;
    mCopy_AB<3>(TipConstraintPoint, ShootingParams.TipConstraintPoint);
    mCopy_AB<3>(TipForce, ShootingParams.TipForce);

}
