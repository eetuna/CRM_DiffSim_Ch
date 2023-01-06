#pragma once
#include "CRM_BVPIVP_APIDeclarations.hpp"
#include "numerical/minpack_Declarations.hpp"
#include "numerical/minpack_DYN.hpp"

#define t_step 0.0005
#include <math.h>
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
void CoilIntegrad(adType in_twist[6], adType in_n[3], adType g[3], adType R[9], adType actMass, adType actInertia[9], adType damping[6], adType in_tau[3], adType twistdot[6]){

    adType v[3], w[3], n_L[3], tau[3];
    adType RTg[3], w_v[3], w_hat[9], inertiaw[3], w_inertia_w[3], diff_tau_w[3];
    adType vdot[3], wdot[3];

     mMult_ATB<3,3,1>(R, g, RTg);

     for (int i = 0; i < 3; ++i) {
        v[i] = in_twist[i];
        w[i] = in_twist[i+3];
    }
     for (int i = 0; i < 3; ++i) {
         n_L[i] = in_n[i];
         tau[i] = in_tau[i];
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


template <typename adType>
void CoilDynamics( adType in_coil_state[NUM_COIL_STATES], adType in_n[3], adType g[3],
                  adType actMass, adType actInertia[9], adType damping[6], double DELTA_T, adType in_tau[3], adType out_coil_state[NUM_COIL_STATES]){

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
    adType nL[3], tau[3];

    for (int i = 0; i < 3; ++i) nL[i] = in_n[i];
    for (int i = 0; i < 3; ++i) tau[i] = in_tau[i];

    // initialize the iteration items
    for (int i = 0; i < NUM_COIL_STATES; i++) {
            x_n[i] = in_coil_state[i];
    }


    int N = ceil(DELTA_T / t_step);

    for (int idx=0; idx<N; idx++) {

        if (idx<3) {  // RK2 initialization steps
            RK2_coildyn(x_n, nL, g,  actMass, actInertia, damping,tau, x_np1, xdot_n );
        }
        else { 		 // ABM4 steps
            ABM4_coildyn(	x_n, xdot_nm1, xdot_nm2, xdot_nm3, x_nm1, x_nm2, x_nm3,
                             nL, g,  actMass, actInertia, damping, tau,x_np1, xdot_n);
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

//        std::cout << "in coil loop:  " << std::endl;
//        for (int i = 0; i < NUM_COIL_STATES; ++i) {
//            std::cout << "x_n: " << x_n[i] << std::endl;
//        }
    }

    // Copy final value
    for (int i=0; i<NUM_COIL_STATES; i++) {
        out_coil_state[i]=x_n[i];
    }

}



//[x_np1, xdot_n] = RK2_step(x_n, t_n, h, Integrand)
template <typename adType>
void RK2_coildyn(adType in_x_n[NUM_COIL_STATES], adType in_n[3], adType g[3],  adType actMass, adType actInertia[9], adType damping[6], adType in_tau[3],
                 adType out_x_np1[NUM_COIL_STATES], adType out_xdot_n[6] ) {


    adType nL[3], tau[3], twist_n[6];       // from input
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
        tau[i] = in_tau[i];
    }

//    std::cout << "input v_n: " << twist_n[0] << " " << twist_n[1] << " " << twist_n[2] <<  std::endl;
//    std::cout << "input w_n: " << twist_n[3] << " " << twist_n[4] << " " << twist_n[5] <<  std::endl;

    //RK2_STEP_STEP1:
    CoilIntegrad(twist_n, nL, g, R_n, actMass, actInertia, damping,tau, xdot_n);
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
    CoilIntegrad(x_n_p_k1o2, nL, g, R_np1half, actMass, actInertia, damping,tau, k2oh);

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
                   adType in_n[3], adType g[3],  adType actMass, adType actInertia[9], adType damping[6], adType in_tau[3],
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
//        x_n[i] = in_x_n[i];
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

    for (int i = 0; i < 3; ++i) {
        nL[i] = in_n[i];
        tau[i] = in_tau[i];
    }

    adType R_np1_hat[9], p_np1_hat[3];

    //ABM4_STEP_STEP1:
    CoilIntegrad(twist_n, nL, g, R_n, actMass, actInertia, damping,tau, xdot_n);

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
    CoilIntegrad(x_np1_hat, nL, g, R_np1_hat, actMass, actInertia, damping,tau, xdot_np1_hat);

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
//    std::cout << " w_n: " << w_n[0] << " " << w_n[1] << " " << w_n[2] <<  std::endl;


    // calculate R_np1 and p_np1 analytically, without numerical integration
    adType Rdelta[9], pdelta[3];
    adType umagsq, umagsqresp, umag, umagresp, unorm[3], delsumag, ImRuxvpuuTvds[3], Rnpd[3];
    umagsq = vNormSq<3>(w_n);

    adType p_dot[3];
    mMult_AB<3,3,1>(R_n, v_n, p_dot);
//    std::cout << "v_n: " << v_n[0] << " " << v_n[1] << " " << v_n[2] <<  std::endl;
//    std::cout << "p_dot: " << p_dot[0] << " " << p_dot[1] << " " << p_dot[2] <<  std::endl;

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
//        std::cout << "pdelta: " << pdelta[0] << " " << pdelta[1] << " " << pdelta[2] <<  std::endl;

        mMult_AB<3, 3, 1>(R_n, pdelta, Rnpd);
        mAdd_AB<3, 1>(p_n, Rnpd, out_p_np1);
    }
}

#undef EPS

template <typename adType>
void DYNNLEquation(adType in_x[], adType out_y[], DYNNLEqnParams<adType> Params, adType out_u0[NUM_FLEX_SEG][3]) {

//    double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3];

    // output for time advance, not used in BVP, just placeholders
    adType m_L[3], n_L[3], ftip[3];
    // don't forget to scale parameters before passing to the CRMSolverIVP
    if (Params.ContactMode == ContactModeType::FREE_TIP) {
        for (int i = 0; i < 3; i++) {
            m_L[i] = IVALUE_SCALE_M * in_x[i];
            n_L[i] = IVALUE_SCALE_N * in_x[i+3];

            ftip[i] = Params.TipForce[i];  // for free-tip, this parameter is not given by the nonlinear equation solver, and hence, does not need to be scaled
        }
    }
    else { // FIXED_TIP
        for (int i = 0; i < 3; i++) {
            m_L[i] = IVALUE_SCALE_M * in_x[i];
            ftip[i] = IVALUE_SCALE_F * in_x[i + 9];
        }
    }

//    std::cout << "m_L dyn: " << m_L[0] << " " << m_L[1] << " " <<m_L[2] << std::endl;
//    std::cout << "n_L dyn: " << n_L[0] << " " << n_L[1] << " " <<n_L[2] << std::endl;

    auto & update_m_L = Params.m_L;
    mCopy_AB<3>(m_L, update_m_L);

    auto & update_n_L = Params.n_L;
    mCopy_AB<3>(n_L, update_n_L);


    adType MomentResidual[3], Tbcoil[3], RL[9], pL[3], Rcoil[9], x_coil[NUM_COIL_STATES], out_x_coil[NUM_COIL_STATES];

    int localmin = 0;
    double u0_calc[3], out_ftip[3];

    StateVector<adType> xi, xf, xi_new;

    mCopy_AB<3>(Params.xi + 0, xi._p);
    mCopy_AB<9>(Params.xi + 3, xi._R);
    mCopy_AB<3>(Params.u0_initialguess, xi._u);


//    double moment_residual[3];

    int ind_s = Params.StartSegmentIndex;
    CRMShootingMethodBVP_DYN(Params, m_L, n_L,xi, ind_s, false, u0_calc, xf, MomentResidual, out_ftip, localmin);
//    std::cout << "u0_calc dyn: " << u0_calc[0] << " " << u0_calc[1] << " " <<u0_calc[2] << std::endl;

    //Pass the state to the rigid link and Get P and R
//    CRMSolverIVP_Core( Params, u0_calc, ftip,  xf, MomentResidual, Tbcoil, pL, RL, p_atLocMarkers);

    ind_s +=1;
    CRMRigidLinkPass(MomentResidual, ind_s, Params, xf, xi_new,
                     Tbcoil, pL, RL);

    ind_s +=1;
    double u0_cal_1[3], deltau[3],  m_L_u[3];
    double m_L_[3] = {0,0,0};
    double n_L_[3] = {0,0,0};
    CRMShootingMethodBVP_DYN(Params, m_L_, n_L_, xi_new, ind_s, true, u0_cal_1, xf, MomentResidual, out_ftip, localmin);
    int  fsegno=ind_s>>1; // i/2, flexible segment no

    mSub_AB<3, 1>(u0_cal_1, Params.ustar[fsegno], deltau);
    mMult_AB<3, 3, 1>(Params.K[fsegno], deltau, m_L_u);
//    std::cout << "u0_cal_free dyn: " << u0_cal_free[0] << " " << u0_cal_free[1] << " " <<u0_cal_free[2] << std::endl;
//    std::cout << "ustar dyn: " << Params.ustar[fsegno][0] << " " << Params.ustar[fsegno][1] << " " <<Params.ustar[fsegno][2] << std::endl;
//
//    std::cout << "m_L_u dyn: " <<m_L_u[0] << " " << m_L_u[1] << " " <<m_L_u[2] << std::endl;

    //    for (int i = 0; i < NUM_STATES; ++i) {
//        std::cout << " x_N: " << x_N[i] << std::endl;
//    }
//    std::cout << "RL dyn: " << RL[0] << " " << RL[1] << " " <<RL[2] << std::endl;
//    std::cout << "RL dyn: " << RL[3] << " " << RL[4] << " " <<RL[5] << std::endl;
//    std::cout << "RL dyn: " << RL[6] << " " << RL[7] << " " <<RL[8] << std::endl;

//    double u0_calc_0[3];
//    mCopy_AB<3>(u0_calc, u0_calc_0);
    for (int i = 0; i < 3; ++i) {
        out_u0[0][i] = u0_calc[i];
    }
    for (int i = 0; i < 3; ++i) {
        out_u0[1][i] = u0_cal_1[i];
    }



    for (int i = 0; i < 3; ++i) {
        x_coil[i] = Params.v_L_pre[i];
        x_coil[i+3] = Params.w_L_pre[i];
        x_coil[i+6] = Params.p_pre[i];
    }
    for (int i = 0; i < 9; ++i) {
        x_coil[i+9] = Params.R_pre[i];
    }

    double actMass = Params.actMass[0];
    double actInertia[9];
    for (int j = 0; j < 9; ++j) {
        actInertia[j] = Params.actInertia[0][j];
    }

    double tau[3], temp[3];
    mSub_AB<3, 1>(Tbcoil, m_L_u, temp);
    mSub_AB<3, 1>(temp, m_L, tau);

    CoilDynamics(x_coil, n_L, Params.g, actMass, actInertia, Params.damping, Params.DELTA_T, tau, out_x_coil);

    double RESIDUAL[NUM_DYN_RESIDUAL];
//    std::cout << "pL: " << pL[0] << " " << pL[1] << " " << pL[2] <<  std::endl;
//    std::cout << "out_x_coil pL: " << out_x_coil[6] << " " << out_x_coil[7] << " " << out_x_coil[8] <<  std::endl;

    for (int i = 0; i < 3; ++i) {
        RESIDUAL[i] = pL[i] - out_x_coil[i+6];
    }

    //rotation residual
    for (int i = 0; i < 9; ++i) {
        Rcoil[i] = out_x_coil[i+9];
    }


    // calculate vector norm of the rotation matrices
    double v1[3], v2[3], v3[3];
    for (int i = 0; i < 3; ++i) {
        v1[i] = Rcoil[i*3] - RL[i*3];
        v2[i] = Rcoil[1+ i*3] - RL[1+ i*3];
        v3[i] = Rcoil[2 + i*3] - RL[2 + i*3];
    }

    RESIDUAL[3] = vNormSq<3>(v1);
    RESIDUAL[4] = vNormSq<3>(v2);
    RESIDUAL[5] = vNormSq<3>(v3);
    for (int i = 0; i < 3; ++i) RESIDUAL[i+3] = sqrt(RESIDUAL[i+3]);

////    euler angles difference
//    double v1[3], v2[3];
//    rotationMatrixToEulerAngles(Rcoil, v1);
//    rotationMatrixToEulerAngles(RL, v2);
//    for (int i = 0; i < 3; ++i) {
//        RESIDUAL[i+3] = v1[i] - v2[i];
//    }

//    std::cout << "p diff: " << RESIDUAL[0] << " " << RESIDUAL[1] << " " << RESIDUAL[2] <<  std::endl;
//    std::cout << "R diff: " << RESIDUAL[3] << " " << RESIDUAL[4] << " " << RESIDUAL[5] <<  std::endl;

    // don't forget to scale parameters before returning to the nonlinear equation solver
    for (int i = 0; i < 3; i++) {
        out_y[i] = RESIDUAL_SCALE_P * RESIDUAL[i]; // RESIDUAL_SCALE_F * WrenchResidual[i] ;
        out_y[i+3] = RESIDUAL_SCALE_R * RESIDUAL[i+3]; // RESIDUAL_SCALE_F * WrenchResidual[i] ;
    }

}


template <typename adType>
void CRMShootingMethodBVP_DYN(	DYNNLEqnParams<adType> in_Params, const double m_L[3], const double n_L[3], const StateVector<adType>& xi_i, int ind_s, bool LastLinkFlexible,  adType out_u0[3],
                                  StateVector<adType> &xf, adType MomentResidual[3], adType out_ftip[3], int& out_localmin) {

//    ContactModeType ContactMode = in_Params.ContactMode;
    int NLEq_Dim = NUM_RESIDUAL;
    // Scale parameters and call the nonlinear equation solver
    const double uscaleinv = 1.0 / IVALUE_SCALE_U;
    const double fscaleinv = 1.0 / IVALUE_SCALE_F;
    auto* initialguessscaled = new double [NLEq_Dim];
    auto* returnedparamscaled = new adType [NLEq_Dim];

    if (in_Params.ContactMode == ContactModeType::FREE_TIP) {
        for (int i = 0; i < 3; i++) {
            initialguessscaled[i] = uscaleinv * in_Params.u0_initialguess[i];
            in_Params.TipForce[i] = in_Params.TipForce[i];  // if the catheter is not in contact, the tip force specified within in_Params needs to be used; this would not be scaled as it is not changed by the solver
        }
    }else { // FIXED_TIP
        for (int i = 0; i < 3; i++) {
            initialguessscaled[i] = uscaleinv * in_Params.u0_initialguess[i];
            initialguessscaled[i+3] = fscaleinv * in_Params.ftip_initialguess[i];
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
#if defined( COIL_TRUSTREGION )
    TrustRegionDogleg_fleseg(NLEq_Dim, x, residual, tol, info, wa, lwa, m_L, n_L, LastLinkFlexible, xi_i, ind_s, in_Params);
#else // undefined
    exit(1);
#endif

    localmin = (info == 1) ? 0 : (info - 1);
    for (int i = 0; i < NLEq_Dim; i++) returnedparamscaled[i] = x[i];
    delete[] residual;
    delete[] x;
    delete[] wa;

    if (in_Params.ContactMode == ContactModeType::FREE_TIP) {
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
//    std::cout << "out_u0: " << out_u0[0] << " " << out_u0[1] << " " << out_u0[2] <<  std::endl;

    out_localmin = localmin;
    delete[] initialguessscaled;
    delete[] returnedparamscaled;

    double l_zero[3];

    CRMSolverIVP_PropagateBCThroughFlexibleLink( xi_i, ind_s, in_Params.SegBounds, in_Params.SegSteps, in_Params.InsertedLength,
                                                 in_Params.dlambdainv, in_Params.K, in_Params.Kinv, l_zero, in_Params.ustar, in_Params.fcumlambda, in_Params.n_L,
                                                 in_Params.m_L, LastLinkFlexible, in_Params.FinalValueOnly, in_Params.LocMarkers, in_Params.NextLocMarker,
                                                 in_Params.p_atLocMarkers, xf, MomentResidual );
//    mCopy_AB<3>(xf._p, x_N + 0 );
//    mCopy_AB<9>(xf._R, x_N + 3 );
//    mCopy_AB<3>(xf._u,x_N + 12);

}

template <typename adType>
void DynamicsBVP(	CRMShootingMethod_DYNParams<adType> in_Params, double in_u0_initialguess[3],
                              double in_mL_initialguess[3], double in_nL_initialguess[3], double in_ftip_initialguess[3],
                              adType out_u0[3], adType out_mL[3], adType out_nL[3], adType out_ftip[3], int& out_localmin) {



    double u0_calc[NUM_FLEX_SEG][3];


    ContactModeType ContactMode = in_Params.ContactMode;
    int NLEq_Dim;  // Dimension of the Nonlinear Equation to Solve
    NLEq_Dim = NUM_DYN_RESIDUAL;

    // Call CRMSolverIVP_Prep, to pre-process parameters
    adType x_0[NUM_STATES];
    for (int i = 0; i < 3; i++) x_0[i] = in_Params.p0[i];
    for (int i = 0; i < 9; i++) x_0[i + 3] = in_Params.R0[i];
    for (int i = 0; i < 3; i++) x_0[i + 3 + 9] = in_u0_initialguess[i];

    bool FinalValueOnly = true;
    DYNNLEqnParams<adType> DYNNLEParams;

    CRMSolverIVP_DYNPrep(x_0, in_Params.IntegrationStepSize,
                      in_Params.Li, in_Params.dlambdainv,
                      in_Params.SegEndLambdas, in_Params.LocMarkerLambdas,
                      in_Params.K, in_Params.Kinv, in_Params.ustar,
                      in_Params.CoilAlignmentTurnAreaMatrix,
                      in_Params.MagMoment, in_Params.fcumlambda,
                      in_Params.B0,in_Params.g, in_Params.actMass, in_Params.actInertia, in_Params.damping, in_Params.DELTA_T,
                      in_Params.v_L_pre, in_Params.w_L_pre, in_Params.p_pre, in_Params.R_pre,
                      in_mL_initialguess, in_nL_initialguess,FinalValueOnly,DYNNLEParams);

    DYNNLEParams.ContactMode = ContactMode;
    mCopy_AB<3>(in_Params.TipConstraintPoint, DYNNLEParams.TipConstraintPoint);
    mCopy_AB<3>(in_u0_initialguess, DYNNLEParams.u0_initialguess);
    mCopy_AB<3>(in_ftip_initialguess, DYNNLEParams.ftip_initialguess);


    // Scale parameters and call the nonlinear equation solver
//    const double uscaleinv = 1.0 / IVALUE_SCALE_U;
    const double nscaleinv = 1.0 / IVALUE_SCALE_N;
    const double mscaleinv = 1.0 / IVALUE_SCALE_M;

    const double fscaleinv = 1.0 / IVALUE_SCALE_F;
    auto* initialguessscaled = new double [NLEq_Dim];
    auto* returnedparamscaled = new adType [NLEq_Dim];

    if (ContactMode == ContactModeType::FREE_TIP) {
        for (int i = 0; i < 3; i++) {
            initialguessscaled[i] = mscaleinv * in_mL_initialguess[i];
            initialguessscaled[i+3] = nscaleinv * in_nL_initialguess[i];

            DYNNLEParams.TipForce[i] = in_Params.TipForce[i];  // if the catheter is not in contact, the tip force specified within in_Params needs to be used; this would not be scaled as it is not changed by the solver
        }
    }else { // FIXED_TIP
        for (int i = 0; i < 3; i++) {
            initialguessscaled[i] = mscaleinv * in_mL_initialguess[i];
            initialguessscaled[i+3] = fscaleinv * in_ftip_initialguess[i];
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
#if defined( COIL_TRUSTREGION )
    TrustRegionDogleg_dyn(NLEq_Dim, x, residual, tol, info, wa, lwa, DYNNLEParams, u0_calc);
#else // undefined
    exit(1);
#endif
    localmin = (info == 1) ? 0 : (info - 1);
    for (int i = 0; i < NLEq_Dim; i++) returnedparamscaled[i] = x[i];
    delete[] residual;
    delete[] x;
    delete[] wa;

    for (int i = 0; i < 3; ++i) {
        out_u0[i] = u0_calc[0][i];
    }

    if (ContactMode == ContactModeType::FREE_TIP) {
        for (int i = 0; i < 3; i++) {
            out_mL[i] = IVALUE_SCALE_M * returnedparamscaled[i];
            out_nL[i] = IVALUE_SCALE_N * returnedparamscaled[i+3];

            out_ftip[i] = in_Params.TipForce[i];  // if it is free-tip, return the tip force specified within in_Params
        }
    }


//    std::cout << "out_u0: " << out_u0[0] << " " << out_u0[1] << " " << out_u0[2] <<  std::endl;
//    std::cout << "out_mL: " << out_mL[0] << " " << out_mL[1] << " " << out_mL[2] <<  std::endl;
//    std::cout << "out_nL: " << out_nL[0] << " " << out_nL[1] << " " << out_nL[2] <<  std::endl;

    out_localmin = localmin;
    delete[] initialguessscaled;
    delete[] returnedparamscaled;

}


template <typename adType>
void DYNSolverIVP(	CRMShootingMethod_DYNParams<adType> in_Params, adType in_u0[3],
                      adType in_mL[3], adType in_nL[3], adType in_ftip[3],
                      bool in_FinalValueOnly,
                      adType out_x_N[NUM_STATES], adType out_coil_state[NUM_COIL_STATES],
                      double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]){

    adType x_0[NUM_STATES];
    for (int i = 0; i < 3; i++) x_0[i] = in_Params.p0[i];
    for (int i = 0; i < 9; i++) x_0[i + 3] = in_Params.R0[i];
    for (int i = 0; i < 3; i++) x_0[i + 3 + 9] = in_u0[i];

    CRMIVPCore_DYNParams<adType> CoreParams;
    adType u_0[3], n_L[3], m_L[3], ftip[3];
    // We need to pass u_0 as input argument as the values in CoreParams will be overriden with the values provided in the input arguments - functionality needed for solving Boundary Value Problems (BVP)
    // copy to local variable
    for (int i = 0; i < 3; i++) u_0[i] = in_u0[i];
    for (int i = 0; i < 3; i++) n_L[i] = in_nL[i];
    for (int i = 0; i < 3; i++) m_L[i] = in_mL[i];
    for (int i = 0; i < 3; i++) ftip[i] = in_ftip[i];

    CRMSolverIVP_DYNPrep(x_0, in_Params.IntegrationStepSize,
                      in_Params.Li, in_Params.dlambdainv,
                      in_Params.SegEndLambdas, in_Params.LocMarkerLambdas,
                      in_Params.K, in_Params.Kinv, in_Params.ustar,
                      in_Params.CoilAlignmentTurnAreaMatrix,
                      in_Params.MagMoment, in_Params.fcumlambda,
                      in_Params.B0,in_Params.g, in_Params.actMass, in_Params.actInertia, in_Params.damping, in_Params.DELTA_T,
                      in_Params.v_L_pre, in_Params.w_L_pre, in_Params.p_pre, in_Params.R_pre,
                      m_L, n_L, in_FinalValueOnly, CoreParams);


    adType Tbcoil[3], pcoil[3], Rcoil[9], MomentResidual[3], x_coil[NUM_COIL_STATES];

    StateVector<adType> x_N;
    CRMSolverIVP_Core ( CoreParams, in_u0, in_mL, in_nL, ftip, x_N, MomentResidual, Tbcoil, pcoil, Rcoil, out_p_atLocMarkers);

    for (int i = 0; i < 3; ++i) {
        x_coil[i] = CoreParams.v_L_pre[i];
        x_coil[i+3] = CoreParams.w_L_pre[i];
        x_coil[i+6] = CoreParams.p_pre[i];
    }
    for (int i = 0; i < 9; ++i) {
        x_coil[i+9] = CoreParams.R_pre[i];
    }

    double actMass = CoreParams.actMass[0];
    double actInertia[9];
    for (int j = 0; j < 9; ++j) {
        actInertia[j] = CoreParams.actInertia[0][j];
    }

    double tau[3];
    mSub_AB<3, 1>(Tbcoil, m_L, tau);

    CoilDynamics(x_coil, n_L, CoreParams.g, actMass, actInertia, CoreParams.damping, CoreParams.DELTA_T, tau, out_coil_state);

    // NO NEED FOR CRMSolverIVP_Return
    mCopy_AB<3>(x_N._p, out_x_N + 0);
    mCopy_AB<9>(x_N._R, out_x_N + 3);
    mCopy_AB<3>(x_N._u, out_x_N + 3 + 9);
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


