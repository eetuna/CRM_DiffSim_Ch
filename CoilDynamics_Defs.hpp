#pragma once
#include "CRMBVP_Defs.hpp"
#include "CRMIVP_Defs.hpp"
#include "minpack_DYN.hpp"

#define t_step 0.0005

 /**
  * Coil integration function in coil frame
  * @tparam adType
  * @param in_twist
  * @param n_L
  * @param g
  * @param R
  * @param actMass
  * @param actInertia
  * @param input_tau
  * @param vdot
  * @param wdot
  */
template <typename adType>
void CoilIntegrad(adType in_twist[6], adType in_n[3], adType g[3], adType R[9], adType actMass, adType actInertia[9], adType damping[6],
                  adType in_B0[3], adType in_muhat[9], adType in_mL[3], adType twistdot[6]){

    adType v[3], w[3], n_L[3], m_L[3], B0[3], muhat[9], Tb[3], tau[3];
    adType RTg[3], RscTB0[3], w_v[3], w_hat[9], inertiaw[3], w_inertia_w[3], diff_tau_w[3];
    adType vdot[3], wdot[3];

     mMult_ATB<3,3,1>(R, g, RTg);

     for (int i = 0; i < 3; ++i) {
        v[i] = in_twist[i];
        w[i] = in_twist[i+3];
    }
     for (int i = 0; i < 3; ++i) {
         n_L[i] = in_n[i];
         B0[i] = in_B0[i];
         m_L[i] = in_mL[i];
     }

     for (int i = 0; i < 9; ++i) {
         muhat[i] = in_muhat[i];
     }

     wHat(w,w_hat);
     mMult_AB<3,3,1>(w_hat, v, w_v);

     adType damping_vec[3];
     for (int i = 0; i < 3; ++i) {
         damping_vec[i] = damping[i] * v[i] ;
     }

     for (int i = 0; i < 3; ++i) {
         vdot[i] = RTg[i]  - n_L[i] / actMass - w_v[i] - damping_vec[i];// / actMass; //RTg[i] - n_L[i] / actMass
     }

//     std::cout << "RTg[i] - n_L[i] / actMass: " << RTg[0] - n_L[0] / actMass << " " << RTg[1] - n_L[1] / actMass << " " << RTg[2] - n_L[2] / actMass <<  std::endl;

     adType damping_wec[3], residual_w[3];
     for (int i = 0; i < 3; ++i) {
         damping_wec[i] = damping[i+3] * w[i] ;
     }

     mMult_AB<3,3,1>(actInertia, w, inertiaw);
     mMult_AB<3,3,1>(w_hat, inertiaw, w_inertia_w);

     mMult_ATB<3,3,1>(R,B0,RscTB0);
     mMult_AB<3,3,1>(muhat,RscTB0,Tb);
     mSub_AB<3, 1>(Tb, m_L, tau);

     mSub_AB<3,1>(tau, w_inertia_w, diff_tau_w);
     mSub_AB<3,1>(diff_tau_w, damping_wec, residual_w);

     ///  actInertia is diagonal
     wdot[0] = residual_w[0] / actInertia[0];
     wdot[1] = residual_w[1] / actInertia[4];
     wdot[2] = residual_w[2] / actInertia[8];

     for (int i = 0; i < 3; ++i) {
         twistdot[i] = vdot[i];
         twistdot[i+3] = wdot[i];
     }

}


/**
 * Main Dynamics function of the coil
 * @tpar
 * am adType
 * @param T
 * @param v_L_pre
 * @param w_L_pre
 * @param in_p
 * @param in_R
 * @param in_n
 * @param g
 * @param actMass
 * @param actInertia
 * @param in_tau
 * @param out_coil_state
 */
template <typename adType>
void CoilDynamics( adType in_coil_state[NUM_COIL_STATES], adType in_n[3], adType g[3],
                  adType actMass, adType actInertia[9], adType damping[6], double DELTA_T, adType in_B0[3], adType in_muhat[9], adType in_mL[3], adType out_coil_state[NUM_COIL_STATES]){

//    for (int i = 0; i < 3; ++i) {
//        std::cout << "in_n: " << in_n[i] << std::endl;
//    }
//    for (int i = 0; i < 3; ++i) {
//        std::cout << "in_tau: " << in_tau[i] << std::endl;
//    }

    adType x_nm3[NUM_COIL_STATES];
    adType x_nm2[NUM_COIL_STATES];
    adType x_nm1[NUM_COIL_STATES];
    adType x_n[NUM_COIL_STATES];
    adType x_np1[NUM_COIL_STATES];
    adType xdot_nm3[6];
    adType xdot_nm2[6];
    adType xdot_nm1[6];
    adType xdot_n[6];
    adType nL[3], B0[3], muhat[9], mL[3];

    for (int i = 0; i < 3; ++i) {
        nL[i] = in_n[i];
        B0[i] = in_B0[i];
        mL[i] = in_mL[i];
    }
    for (int i = 0; i < 9; ++i) {
        muhat[i] = in_muhat[i];
    }


    // initialize the iteration items
    for (int i = 0; i < NUM_COIL_STATES; i++) {
            x_n[i] = in_coil_state[i];
    }


    int N = ceil(DELTA_T / t_step);

    for (int idx=0; idx<N; idx++) {

        if (idx<3) {  // RK2 initialization steps
            RK2_coildyn(x_n, nL, g,  actMass, actInertia, damping, B0, muhat, mL,x_np1, xdot_n );
        }
        else { 		 // ABM4 steps
            ABM4_coildyn(	x_n, xdot_nm1, xdot_nm2, xdot_nm3, x_nm1, x_nm2, x_nm3,
                             nL, g,  actMass, actInertia, damping, B0, muhat, mL,x_np1, xdot_n);
            if ( isnan(x_n[0]) ) {
                std::cout << "FLY ME TO THE MOON!! " << std::endl;
//                exit( 3 );
            }
        }

        // update the iteration items
        for (int i = 0; i < NUM_COIL_STATES; i++) {
            x_nm3[i] = x_nm2[i];
            x_nm2[i] = x_nm1[i];
            x_nm1[i] = x_n[i];
            x_n[i] = x_np1[i];
        }
        for (int i = 0; i < 6; i++) {
            xdot_nm3[i]=xdot_nm2[i];
            xdot_nm2[i]=xdot_nm1[i];
            xdot_nm1[i]=xdot_n[i];
        }

    }

    // Copy final value
    for (int i=0; i<NUM_COIL_STATES; i++) {
        out_coil_state[i]=x_n[i];
    }

}



//[x_np1, xdot_n] = RK2_step(x_n, t_n, h, Integrand)
template <typename adType>
void RK2_coildyn(adType in_x_n[NUM_COIL_STATES], adType in_n[3], adType g[3],  adType actMass, adType actInertia[9], adType damping[6], adType in_B0[3], adType in_muhat[9], adType in_mL[3],
                 adType out_x_np1[NUM_COIL_STATES], adType out_xdot_n[6] ) {


    adType nL[3], B0[3], muhat[9], mL[3], twist_n[6];       // from input
    adType k1[6];
    adType k2oh[6];
    adType x_n_p_k1o2[6];
    adType xdot_n[6];

    adType R_n[9], p_n[3];				// variables used in analytical calculation
    for (int i = 0; i < 3; ++i) p_n[i] = in_x_n[i+6];
    for (int i = 0; i < 9; ++i) R_n[i] = in_x_n[i+9];
    for (int i = 0; i < 6; ++i) twist_n[i] = in_x_n[i];

    for (int i = 0; i < 3; ++i) {
        nL[i] = in_n[i];
        B0[i] = in_B0[i];
        mL[i] = in_mL[i];
    }
    for (int i = 0; i < 9; ++i) {
        muhat[i] = in_muhat[i];
    }

//    std::cout << "input v_n: " << twist_n[0] << " " << twist_n[1] << " " << twist_n[2] <<  std::endl;
//    std::cout << "input w_n: " << twist_n[3] << " " << twist_n[4] << " " << twist_n[5] <<  std::endl;

    //RK2_STEP_STEP1:
    CoilIntegrad(twist_n, nL, g, R_n, actMass, actInertia, damping, B0, muhat, mL, xdot_n);
    for (int i=0; i< 6 ; i++) {
        k1[i] 			= t_step * xdot_n[i];
        x_n_p_k1o2[i] 	= twist_n[i] + k1[i] * 0.5;
    }
#ifdef ANALYTICAL_SE3_STEP
    // we will calculate R_np1half and p_np1half analytically, without numerical integration
    adType R_np1half[9], p_np1half[3];

    SE3_TimeSpace(R_n, p_n, t_step*0.5, twist_n, R_np1half, p_np1half) ;


#endif

    //RK2_STEP_STEP2:
    CoilIntegrad(x_n_p_k1o2, nL, g, R_np1half, actMass, actInertia, damping,B0, muhat, mL, k2oh);

    for (int i = 0; i < 6; i++) {
        out_x_np1[i] = twist_n[i] + t_step * k2oh[i];
    }

#ifdef ANALYTICAL_SE3_STEP
    // we will calculate R_np1 and p_np1 analytically, without numerical integration
    adType R_np1[9], p_np1[3];

    SE3_TimeSpace(R_n, p_n, t_step,  x_n_p_k1o2, R_np1, p_np1) ;

    // copy these to the output state
    for (int i = 0; i < 3; i++) out_x_np1[i+6] = p_np1[i];
    for (int i = 0; i < 9; i++) out_x_np1[i+9] = R_np1[i];

#endif
    for (int i = 0; i < 6; i++) {
        out_xdot_n[i] = xdot_n[i];
    }

}


// [x_np1, xdot_n, xdot_nm1, xdot_nm2] = ABM4_step(x_n, t_n, xdot_nm1, xdot_nm2, xdot_nm3, h, Integrand)
template <typename adType>
void ABM4_coildyn(	adType in_x_n[NUM_COIL_STATES],adType in_xdot_nm1[6], adType in_xdot_nm2[6], adType in_xdot_nm3[6],
                   adType in_x_nm1[NUM_COIL_STATES], adType in_x_nm2[NUM_COIL_STATES], adType in_x_nm3[NUM_COIL_STATES],
                   adType in_n[3], adType g[3],  adType actMass, adType actInertia[9], adType damping[6], adType in_B0[3], adType in_muhat[9], adType in_mL[3],
                   adType out_x_np1[NUM_COIL_STATES], adType out_xdot_n[6]) {

    const double P_COEFF_N=55.0/24.0, P_COEFF_Nm1=-59.0/24.0, P_COEFF_Nm2=37.0/24.0, P_COEFF_Nm3=-9.0/24.0;  // AB4 Predictor Coefficients
    const double C_COEFF_Np1=9.0/24.0, C_COEFF_N=19.0/24.0, C_COEFF_Nm1=-5.0/24.0, C_COEFF_Nm2=1.0/24.0;     // AM4 Corrector Coefficients
//    adType x_n[NUM_COIL_STATES];       				// from input
    adType twist_nm1[6];       				// from input
    adType twist_nm2[6];       				// from input
    adType twist_nm3[6];       				// from input
    adType twist_n[6], R_n[9], p_n[3], nL[3], tau[3];				// variables used in analytical calculation
    adType x_np1_hat[6];    			// intermediate twist
    adType xdot_np1_hat[6];	// intermediate
    adType xdot_n[6];  		// for output
    adType xdot_nm1[6];  	// from input
    adType xdot_nm2[6];  	// from input
    adType xdot_nm3[6];  	// from input

    for (int i = 0; i < 6; i++) {
        twist_nm1[i] = in_x_nm1[i];
        twist_nm2[i] = in_x_nm2[i];
        twist_nm3[i] = in_x_nm3[i];
    }
    for (int i = 0; i < 6; i++) {
        xdot_nm1[i]=in_xdot_nm1[i];
        xdot_nm2[i]=in_xdot_nm2[i];
        xdot_nm3[i]=in_xdot_nm3[i];
    }


    for (int i = 0; i < 3; ++i) p_n[i] = in_x_n[i+6];
    for (int i = 0; i < 9; ++i) R_n[i] = in_x_n[i+9];
    for (int i = 0; i < 6; ++i) twist_n[i] = in_x_n[i];

    adType B0[3], muhat[9], mL[3];
    for (int i = 0; i < 3; ++i) {
        nL[i] = in_n[i];
        B0[i] = in_B0[i];
        mL[i] = in_mL[i];
    }
    for (int i = 0; i < 9; ++i) {
        muhat[i] = in_muhat[i];
    }


    adType R_np1_hat[9], p_np1_hat[3];

    //ABM4_STEP_STEP1:
    CoilIntegrad(twist_n, nL, g, R_n, actMass, actInertia, damping,B0, muhat, mL,xdot_n);

    for (int i=0; i<6; i++) {
        x_np1_hat[i]    = twist_n[i] + t_step * ( P_COEFF_N * xdot_n[i] + P_COEFF_Nm1 * xdot_nm1[i] + P_COEFF_Nm2 * xdot_nm2[i] + P_COEFF_Nm3 * xdot_nm3[i] );
    }
#ifdef ANALYTICAL_SE3_STEP
    //      calculate R_np1_hat and p_np1_hat analytically, without numerical integration
    adType twist_n_pred[6];
    for (int i = 0; i < 6; i++) twist_n_pred[i] = (P_COEFF_N * twist_n[i] + P_COEFF_Nm1 * twist_nm1[i] + P_COEFF_Nm2 * twist_nm2[i] + P_COEFF_Nm3 * twist_nm3[i]);
    SE3_TimeSpace(R_n, p_n, t_step, twist_n, R_np1_hat, p_np1_hat) ;

//    for (int i = 0; i < 3; ++i) x_np1_hat[i+6] = p_np1_hat[i];
//    for (int i = 0; i < 9; ++i) x_np1_hat[i+9] = R_np1_hat[i];
#endif
    //ABM4_STEP_STEP2:
    CoilIntegrad(x_np1_hat, nL, g, R_np1_hat, actMass, actInertia, damping,B0, muhat, mL,xdot_np1_hat);

    for (int i=0; i<6; i++) {
        out_x_np1[i]    = twist_n[i] + t_step * ( C_COEFF_Np1 * xdot_np1_hat[i] + C_COEFF_N * xdot_n[i] + C_COEFF_Nm1 * xdot_nm1[i] + C_COEFF_Nm2 * xdot_nm2[i] );
    }
#ifdef ANALYTICAL_SE3_STEP
    //      calculate R_np1 and p_np1 analytically, without numerical integration
    adType twist_n_corr[6], R_x_np1[9], p_x_np1[3];
    for (int i = 0; i < 6; i++) twist_n_corr[i] = (C_COEFF_Np1 * x_np1_hat[i] + C_COEFF_N * twist_n[i] + C_COEFF_Nm1 * twist_nm1[i] + C_COEFF_Nm2 * twist_nm2[i]);
    SE3_TimeSpace(R_n, p_n, t_step, twist_n_corr, R_x_np1, p_x_np1) ;
    for (int i = 0; i < 3; ++i) out_x_np1[i+6] = p_x_np1[i];
    for (int i = 0; i < 9; ++i) out_x_np1[i+9] = R_x_np1[i];
#endif

    // copy to output
    for (int i=0; i<6; i++) {
        out_xdot_n[i] 	= xdot_n[i];
    }

}


// definitions needed for twist exponential calculation
#define EPS 1.0e-12   // the threshold for assuming ||u|| to be approximately 0, so that we should use pure translation equation


// calculate R_np1 and p_np1 analytically using twist exponential, without numerical integration
//   g_np1 = g_n * expm ( \hat{\xi}^b *h ),  where \xi^b= [ 0 0 1 u_n^T ]^T
//   g = [R p; 0 0 0 1];
template <typename adType>
void SE3_TimeSpace(adType in_R_n[9], adType in_p_n[3], adType h, adType in_twist_n[6], adType out_R_np1[9], adType out_p_np1[3]) {
    adType R_n[9], p_n[3], v_n[3], w_n[3];

    mCopy_AB<3*3>(in_R_n, R_n);
    mCopy_AB<3>(in_p_n, p_n);

    for (int i = 0; i < 3; ++i) {
        v_n[i] = in_twist_n[i];
        w_n[i] = in_twist_n[i+3];
    }

    // calculate R_np1 and p_np1 analytically, without numerical integration
    adType Rdelta[9], pdelta[3];
    adType umagsq, umagsqresp, umag, umagresp, unorm[3], delsumag, ImRuxvpuuTvds[3], Rnpd[3];
    umagsq = vNormSq<3>(w_n);

    adType p_dot[3];
    mMult_AB<3,3,1>(R_n, v_n, p_dot);

    if (umagsq < EPS) {
        mCopy_AB<3*3>(R_n, out_R_np1);
        out_p_np1[0] = p_n[0] + p_dot[0] * h;
        out_p_np1[1] = p_n[1] + p_dot[1] * h;
        out_p_np1[2] = p_n[2] + p_dot[2] * h;
    }
    else {
        umagsqresp = 1.0 / umagsq;
        umag = sqrt(umagsq);
        umagresp = 1.0 / umag;
        mMult_sA<3, 1>(umagresp, w_n, unorm);
        delsumag = h * umag;
        RodriguesExpanded(unorm, delsumag, Rdelta);
        mMult_AB<3, 3, 3>(R_n, Rdelta, out_R_np1);

        adType what[9], wxv[3], ImRdelta[9], ImRdeltawxv[3] , wwTv[3], wwTvh[3];
        wHat(w_n, what);
        mMult_AB<3,3,1>(what, v_n, wxv);

        ImRdelta[0] = 1 - Rdelta[0]; ImRdelta[1] = -Rdelta[1]; ImRdelta[2] = -Rdelta[2];
        ImRdelta[3] =  -Rdelta[3]; ImRdelta[4] = 1 -Rdelta[4]; ImRdelta[5] = -Rdelta[5];
        ImRdelta[6] = -Rdelta[6]; ImRdelta[7] = -Rdelta[7]; ImRdelta[8] = 1 - Rdelta[8];

        mMult_AB<3,3,1>(ImRdelta, wxv, ImRdeltawxv);

        wwTv[0] = w_n[0]* w_n[0] * v_n[0] + w_n[0]*w_n[1] * v_n[1] + w_n[0] * w_n[2] * v_n[2];
        wwTv[1] = w_n[1]* w_n[0] * v_n[1] + w_n[1]*w_n[1] * v_n[1] + w_n[1] * w_n[2] * v_n[2];
        wwTv[2] = w_n[2]* w_n[0] * v_n[0] + w_n[2]*w_n[1] * v_n[1] + w_n[2] * w_n[2] * v_n[2];

        for (int i = 0; i < 3; ++i) {
            wwTvh[i] = wwTv[i] * h;
        }

        mAdd_AB<3, 1>(ImRdeltawxv, wwTvh, ImRuxvpuuTvds);
        mMult_sA<3, 1>(umagsqresp, ImRuxvpuuTvds, pdelta);

        mMult_AB<3, 3, 1>(R_n, pdelta, Rnpd);
        mAdd_AB<3, 1>(p_n, Rnpd, out_p_np1);
    }
}

#undef EPS

template <typename adType>
void DYNNLEquation(adType in_x[], adType out_y[], NLEqnParams<adType> Params, adType out_u0[NUM_FLEX_SEG*3]) {

    double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3];

    // output for time advance, not used in BVP, just placeholders
    adType m_L[NUM_ACT_SET][3], n_L[NUM_ACT_SET][3], ftip[3];
    // don't forget to scale parameters before passing to the CRMSolverIVP
    if (Params.ContactMode == ContactModeType::FREE_TIP) {
        for (int i = 0; i < NUM_ACT_SET; ++i) {
            for (int j = 0; j < 3; ++j) {
                m_L[i][j] = IVALUE_SCALE_M * in_x[j+ 6*i];
                n_L[i][j] = IVALUE_SCALE_N * in_x[j+3 + 6*i];
            }
        }
        for (int i = 0; i < 3; i++) {
            ftip[i] = Params.TipForce[i];  // for free-tip, this parameter is not given by the nonlinear equation solver, and hence, does not need to be scaled
        }
    }

//    std::cout << "m_L dyn: " << m_L[0] << " " << m_L[1] << " " <<m_L[2] << std::endl;
//    std::cout << "n_L dyn: " << n_L[0] << " " << n_L[1] << " " <<n_L[2] << std::endl;

    auto & update_m_L = Params.m_L;
    auto & update_n_L = Params.n_L;

    for (int i = 0; i < NUM_ACT_SET; ++i) {
        for (int j = 0; j < 3; ++j) {
            update_m_L[i][j] = m_L[i][j];
            update_n_L[i][j] = n_L[i][j];
        }
    }

    adType xf[NUM_STATES], MomentResidual[NUM_RESIDUAL], RL[NUM_ACT_SET][9], pL[NUM_ACT_SET][3], Rcoil[9], x_coil[NUM_COIL_STATES], out_x_coil[NUM_COIL_STATES];

    int localmin = 0;
    double u0_calc[NUM_FLEX_SEG][3], out_ftip[3];
    CRMShootingMethodBVP_DYN(Params, u0_calc, out_ftip, localmin);
//    std::cout << "u0_calc dyn: " << u0_calc[0] << " " << u0_calc[1] << " " <<u0_calc[2] << std::endl;

    //Solve the Initial Value Problem to calculate the shape of the catheter and Get P and R
    CRMSolverIVP_Core ( Params, u0_calc, ftip,  xf, MomentResidual, pL, RL, p_atLocMarkers);

    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
        for (int j = 0; j < 3; ++j) {
            out_u0[j+3*i] = u0_calc[i][j];
        }
    }
//    mCopy_AB<3>(u0_calc, out_u0);
    double actMass, actInertia[9], damping[6], residual_per_coil[6], n_L_p[3], m_L_p[3], MagMoment[3], muhat[9];
    double v1[3], v2[3], v3[3];
    double RESIDUAL[NUM_DYN_RESIDUAL];
//    std::cout << "MomentResidual: " << MomentResidual[0] << " " << MomentResidual[1] << " " << MomentResidual[2] <<  std::endl;
    for (int i = 0; i < NUM_ACT_SET; ++i) {
        for (int j = 0; j < 3; ++j) {
            x_coil[j] = Params.v_L_pre[i][j];
            x_coil[j+3] = Params.w_L_pre[i][j];
            x_coil[j+6] = Params.p_pre[i][j];
        }
        for (int j = 0; j < 9; ++j) {
            x_coil[j+9] = Params.R_pre[i][j];
        }

        actMass = Params.actMass[i];
        for (int j = 0; j < 9; ++j) {
            actInertia[j] = Params.actInertia[i][j];
        }
    //    std::cout << "in x_coil PL: " << x_coil[6] << " " << x_coil[7] << " " << x_coil[8] <<  std::endl;
        for (int j = 0; j < 6; ++j) {
            damping[j] = Params.damping[i][j];
        }

        for (int j = 0; j < 3; ++j) {
            MagMoment[j] = Params.MagMoment[i][j];
        }
        wHat(MagMoment,muhat);

//        for (int j = 0; j < 3; ++j) {
//            tau[j] = Tbcoil[i][j] - m_L[i][j];
//        }

        for (int j = 0; j < 3; ++j) {
            n_L_p[j] = n_L[i][j];
        }
        for (int j = 0; j < 3; ++j) {
            m_L_p[j] = m_L[i][j];
        }

        CoilDynamics(x_coil, n_L_p, Params.g, actMass, actInertia, damping, Params.DELTA_T, Params.B0, muhat,  m_L_p,out_x_coil);

        for (int j = 0; j < 3; ++j) {
            residual_per_coil[j] = pL[i][j] - out_x_coil[j+6];
        }

        for (int j = 0; j < 9; ++j) {
            Rcoil[j] = out_x_coil[j+9];
        }
        for (int j = 0; j < 3; ++j) {
            v1[j] = Rcoil[j*3] - RL[i][j*3];
            v2[j] = Rcoil[1+j*3] - RL[i][1+ j*3];
            v3[j] = Rcoil[2+j*3] - RL[i][2 + j*3];
        }
        residual_per_coil[3] = vNormSq<3>(v1);
        residual_per_coil[4] = vNormSq<3>(v2);
        residual_per_coil[5] = vNormSq<3>(v3);
        for (int j = 0; j < 3; ++j) residual_per_coil[j+3] = sqrt(residual_per_coil[j+3]);

        for (int j = 0; j < 3; ++j) {
            RESIDUAL[j + i*6] = residual_per_coil[j];
            RESIDUAL[j + 3+ i*6] = residual_per_coil[j+3];
        }
    }

    for (int i = 0; i < NUM_ACT_SET; ++i) {
        for (int j = 0; j < 3; ++j) {
            out_y[j+ 6*i] = RESIDUAL_SCALE_P * RESIDUAL[j+ 6*i]; // RESIDUAL_SCALE_F * WrenchResidual[i] ;
            out_y[j+3+ 6*i] = RESIDUAL_SCALE_R * RESIDUAL[j+3+6*i]; // RESIDUAL_SCALE_F * WrenchResidual[i] ;

        }
    }

}


template <typename adType>
void CRMShootingMethodBVP_DYN(	NLEqnParams<adType> in_Params, adType out_u0[NUM_FLEX_SEG][3], adType out_ftip[3], int& out_localmin) {

//    ContactModeType ContactMode = in_Params.ContactMode;
    int NLEq_Dim = NUM_RESIDUAL;
    // Scale parameters and call the nonlinear equation solver
    const double uscaleinv = 1.0 / IVALUE_SCALE_U;
    const double fscaleinv = 1.0 / IVALUE_SCALE_F;
    auto* initialguessscaled = new double [NLEq_Dim];
    auto* returnedparamscaled = new adType [NLEq_Dim];

    if (in_Params.ContactMode == ContactModeType::FREE_TIP) {
        for (int i = 0; i < NUM_FLEX_SEG; ++i) {
            for (int j = 0; j < 3; ++j) {
                initialguessscaled[j+3*i] = uscaleinv * in_Params.u0_initialguess[i][j];
            }
        }
        for (int i = 0; i < 3; i++) {
            in_Params.TipForce[i] = in_Params.TipForce[i];  // if the catheter is not in contact, the tip force specified within in_Params needs to be used; this would not be scaled as it is not changed by the solver
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
    TrustRegionDogleg(NLEq_Dim, x, residual, tol, info, wa, lwa, in_Params);
#else // undefined
    exit(1);
#endif
    localmin = (info == 1) ? 0 : (info - 1);
    for (int i = 0; i < NLEq_Dim; i++) returnedparamscaled[i] = x[i];
    delete[] residual;
    delete[] x;
    delete[] wa;

    if (in_Params.ContactMode == ContactModeType::FREE_TIP) {
        for (int i = 0; i < NUM_FLEX_SEG; ++i) {
            for (int j = 0; j < 3; ++j) {
                out_u0[i][j] = IVALUE_SCALE_U * returnedparamscaled[j+3*i];
            }
        }
        for (int i = 0; i < 3; i++) {
            out_ftip[i] = in_Params.TipForce[i];  // if it is free-tip, return the tip force specified within in_Params
        }
    }


    out_localmin = localmin;
    delete[] initialguessscaled;
    delete[] returnedparamscaled;

}

template <typename adType>
void DynamicsBVP(	CRMShootingMethodParams<adType> in_Params, double in_u0_initialguess[NUM_FLEX_SEG][3],
                              double in_mL_initialguess[NUM_ACT_SET][3], double in_nL_initialguess[NUM_ACT_SET][3], double in_ftip_initialguess[3],
                              adType out_u0[NUM_FLEX_SEG][3], adType out_mL[NUM_ACT_SET][3], adType out_nL[NUM_ACT_SET][3], adType out_ftip[3], int& out_localmin) {

    ContactModeType ContactMode = in_Params.ContactMode;
    int NLEq_Dim;  // Dimension of the Nonlinear Equation to Solve
    NLEq_Dim = NUM_DYN_RESIDUAL;

    // Call CRMSolverIVP_Prep, to pre-process parameters
//    adType x_0[NUM_STATES], u0_calc[NUM_FLEX_SEG*3];

    auto* x_0 = new adType[NUM_STATES];
    auto* u0_calc = new adType[NUM_FLEX_SEG*3];

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
    NLEqnParams<adType> DYNNLEParams;

    CRMSolverIVP_Prep(x_0, in_Params.IntegrationStepSize,
                      in_Params.Li, in_Params.dlambdainv,
                      in_Params.SegEndLambdas, in_Params.LocMarkerLambdas,
                      in_Params.K, in_Params.Kinv, in_Params.ustar,
                      in_Params.MagMoment, in_Params.fcumlambda,
                      in_Params.B0, in_Params.g, in_Params.actMass, in_Params.actInertia, in_Params.damping, in_Params.DELTA_T,
                      in_Params.v_L_pre, in_Params.w_L_pre, in_Params.p_pre, in_Params.R_pre,
                      in_mL_initialguess, in_nL_initialguess,true, FinalValueOnly,DYNNLEParams);

    DYNNLEParams.ContactMode = ContactMode;
    mCopy_AB<3>(in_Params.TipConstraintPoint, DYNNLEParams.TipConstraintPoint);
    mCopy_AB<3>(in_ftip_initialguess, DYNNLEParams.ftip_initialguess);

//    mCopy_AB<3>(in_u0_initialguess, DYNNLEParams.u0_initialguess);
    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
        for (int j = 0; j < 3; ++j) {
            DYNNLEParams.u0_initialguess[i][j] = in_u0_initialguess[i][j];
        }
    }

    // Scale parameters and call the nonlinear equation solver
//    const double uscaleinv = 1.0 / IVALUE_SCALE_U;
    const double nscaleinv = 1.0 / IVALUE_SCALE_N;
    const double mscaleinv = 1.0 / IVALUE_SCALE_M;

    const double fscaleinv = 1.0 / IVALUE_SCALE_F;
    auto* initialguessscaled = new double [NLEq_Dim];
    auto* returnedparamscaled = new adType [NLEq_Dim];

    if (ContactMode == ContactModeType::FREE_TIP) {
        for (int i = 0; i < NUM_ACT_SET; ++i) {
            for (int j = 0; j < 3; j++) {
                initialguessscaled[j+ 6*i] = mscaleinv * in_mL_initialguess[i][j];
                initialguessscaled[j+3 + 6*i] = nscaleinv * in_nL_initialguess[i][j];
            }
        }
        for (int i = 0; i < 3; i++) {
            DYNNLEParams.TipForce[i] = in_Params.TipForce[i];  // if the catheter is not in contact, the tip force specified within in_Params needs to be used; this would not be scaled as it is not changed by the solver
        }
    }

    int localmin = 0;

    adType* x = new adType[NLEq_Dim]; // we will create a new variable here and not use initial guess scaled since truss-region-dogleg algorithm uses the same variable for both input and output

    adType* residual = new adType[NLEq_Dim];
    int info;
    int lwa = (NLEq_Dim * (3 * NLEq_Dim + 13)) / 2;
    double tol = 0.00001;
    adType* wa = new adType [lwa];
    for (int i = 0; i < NLEq_Dim; i++) x[i] = initialguessscaled[i];
#if defined( TRUSTREGION )
    TrustRegionDogleg_dyn(NLEq_Dim, x, residual, tol, info, wa, lwa, DYNNLEParams, u0_calc);
#else // undefined
    exit(1);
#endif
    localmin = (info == 1) ? 0 : (info - 1);
    for (int i = 0; i < NLEq_Dim; i++) returnedparamscaled[i] = x[i];
    delete[] residual;
    delete[] x;
    delete[] wa;

    for (int i = 0; i < NUM_ACT_SET; ++i) {
        for (int j = 0; j < 3; j++) {
            out_mL[i][j] = IVALUE_SCALE_M * returnedparamscaled[j+ 6*i];
            out_nL[i][j] = IVALUE_SCALE_N * returnedparamscaled[j+3 + 6*i];
        }
    }
    for (int i = 0; i < 3; i++) {
        out_ftip[i] = in_Params.TipForce[i];  // if it is free-tip, return the tip force specified within in_Params
    }

    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
        for (int j = 0; j < 3; ++j) {
            out_u0[i][j] = u0_calc[j+ 3*i];
        }
    }


    delete[] x_0;
    delete[] u0_calc;

    out_localmin = localmin;
    delete[] initialguessscaled;
    delete[] returnedparamscaled;

}


template <typename adType>
void DYNSolverIVP(	CRMShootingMethodParams<adType> in_Params, adType in_u0[NUM_FLEX_SEG][3],
                      adType in_mL[NUM_ACT_SET][3], adType in_nL[NUM_ACT_SET][3], adType in_ftip[3],
                      bool in_FinalValueOnly, adType out_x_N[NUM_STATES], adType out_coil_state[NUM_COIL_STATES],
                      double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]){

    adType x_0[NUM_STATES];
    for (int i = 0; i < NUM_STATES; i++) {
        if (i < 3) x_0[i] = in_u0[0][i];
        else if (i < 12) x_0[i] = in_Params.R0[i - 3];
        else if (i < 15) x_0[i] = in_Params.p0[i - 12];
    }

    CRMIVPCoreParams<adType> CoreParams;
    adType u_0[NUM_FLEX_SEG][3], n_L[NUM_ACT_SET][3], m_L[NUM_ACT_SET][3], ftip[3];
    // We need to pass u_0 as input argument as the values in CoreParams will be overriden with the values provided in the input arguments - functionality needed for solving Boundary Value Problems (BVP)
    // copy to local variable
    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
        for (int j = 0; j < 3; j++) u_0[i][j] = in_u0[i][j];
    }
    for (int i = 0; i < NUM_ACT_SET; ++i) {
        for (int j = 0; j < 3; ++j) {
            n_L[i][j] = in_nL[i][j];
            m_L[i][j] = in_mL[i][j];
        }
    }
    for (int i = 0; i < 3; i++) ftip[i] = in_ftip[i];

    CRMSolverIVP_Prep(x_0, in_Params.IntegrationStepSize,
                      in_Params.Li, in_Params.dlambdainv,
                      in_Params.SegEndLambdas, in_Params.LocMarkerLambdas,
                      in_Params.K, in_Params.Kinv, in_Params.ustar,
                      in_Params.MagMoment, in_Params.fcumlambda,
                      in_Params.B0, in_Params.g, in_Params.actMass, in_Params.actInertia, in_Params.damping, in_Params.DELTA_T,
                      in_Params.v_L_pre, in_Params.w_L_pre, in_Params.p_pre, in_Params.R_pre,
                      m_L, n_L, true, in_FinalValueOnly,CoreParams);

    adType pcoil[NUM_ACT_SET][3], Rcoil[NUM_ACT_SET][9], MomentResidual[3], x_coil[NUM_COIL_STATES];
    adType MagMoment[3] ,muhat[9];

    CRMSolverIVP_Core ( CoreParams, u_0, ftip, out_x_N, MomentResidual, pcoil, Rcoil, out_p_atLocMarkers);

    for (int i = 0; i < 3; ++i) {
        x_coil[i] = CoreParams.v_L_pre[0][i];
        x_coil[i+3] = CoreParams.w_L_pre[0][i];
        x_coil[i+6] = CoreParams.p_pre[0][i];
    }
    for (int i = 0; i < 9; ++i) {
        x_coil[i+9] = CoreParams.R_pre[0][i];
    }

    double actMass = CoreParams.actMass[0];
    double actInertia[9];
    for (int j = 0; j < 9; ++j) {
        actInertia[j] = CoreParams.actInertia[0][j];
    }
    double damping[6];
    for (int i = 0; i < 6; ++i) {
        damping[i] = CoreParams.damping[0][i];
    }

//    double tau[3];
//    for (int i = 0; i < 3; ++i) {
//        tau[i] = Tbcoil[0][i] - m_L[0][i];
//    }
    for (int i = 0; i <3 ; ++i) {
        MagMoment[i] = CoreParams.MagMoment[0][i];
    }
    wHat(MagMoment,muhat);


    CoilDynamics(x_coil, n_L[0], CoreParams.g, actMass, actInertia, damping, CoreParams.DELTA_T, in_Params.B0, muhat, m_L[0], out_coil_state);

}

template <typename adType>
void rotationMatrixToEulerAngles(adType in_R[9], adType out_v[3]){
    double sy = sqrt(in_R[0] * in_R[0] + in_R[3] * in_R[3] );
    bool singular = sy < 1e-6; // If

    double x, y, z;
    if (!singular){
        x = atan2(in_R[7] , in_R[8]);
        y = atan2(-in_R[6], sy);
        z = atan2(in_R[3], in_R[0]);
    }else{
        std::cout << "Singular Rotation matrix............" << std::endl;
        x = atan2(-in_R[5], in_R[4]);
        y = atan2(-in_R[6], sy);
        z = 0;
    }

    out_v[0] = x; out_v[1] = y; out_v[2] = z;

}