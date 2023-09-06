#pragma once
#include "CRMBVP_Defs.hpp"
#include "CRMIVP_Defs.hpp"
#include "Numerical/minpack_DYN.hpp"

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
    adType RscTB0[3], RTg[3], w_v[3], w_hat[9], inertiaw[3], w_inertia_w[3], diff_tau_w[3];
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
                  adType actMass, adType actInertia[9], adType damping[6], double DELTA_T, adType in_B0[3], adType in_muhat[9],
                  adType in_mL[3], adType out_coil_state[NUM_COIL_STATES], adType out_xdot_n[6]){

//    for (int i = 0; i < 3; ++i) {
//        std::cout << "in_n: " << in_n[i] << std::endl;
//    }
//    for (int i = 0; i < 3; ++i) {
//        std::cout << "in_tau: " << in_tau[i] << std::endl;
//    }

//
//        std::cout << "R_pre: " << std::endl;
//        std::cout <<  in_coil_state[9] << " " << in_coil_state[10]<< " " << in_coil_state[11] <<  std::endl;
//        std::cout <<  in_coil_state[12]<< " " << in_coil_state[13]<< " " << in_coil_state[14] <<  std::endl;
//        std::cout << in_coil_state[15]<< " " << in_coil_state[16] << " " << in_coil_state[17] <<  std::endl;
//


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
//        std::cout << "input v_n: " << x_n[0] << " " << x_n[1] << " " << x_n[2] <<  std::endl;
//         std::cout << "input w_n: " << x_n[3] << " " << x_n[4] << " " << x_n[5] <<  std::endl;
//

    int N = ceil(DELTA_T / t_step);

    for (int idx=0; idx<N; idx++) {

        if (idx<3) {  // RK2 initialization steps
            RK2_coildyn(x_n, nL, g,  actMass, actInertia, damping, B0, muhat, mL, x_np1, xdot_n );
        }
        else { 		 // ABM4 steps
            ABM4_coildyn(	x_n, xdot_nm1, xdot_nm2, xdot_nm3, x_nm1, x_nm2, x_nm3,
                             nL, g,  actMass, actInertia, damping,  B0, muhat, mL,x_np1, xdot_n);
            if ( isnan(x_n[0]) ) {
                std::cout << "Coil integration Unbounded!! " << std::endl;
//                exit( 3 );
            }
        }

//        std::cout << "input v_n: " << xdot_n[0] << " " << xdot_n[1] << " " << xdot_n[2] <<  std::endl;
//         std::cout << "input w_n: " << xdot_n[3] << " " << xdot_n[4] << " " << xdot_n[5] <<  std::endl;

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
    for (int i = 0; i < 6; ++i) {
        out_xdot_n[i] = xdot_n[i];
    }
//
//    std::cout << "RL_now: " << std::endl;
//    std::cout <<  x_n[9] << " " << x_n[10]<< " " << x_n[11] <<  std::endl;
//    std::cout <<  x_n[12]<< " " << x_n[13]<< " " << x_n[14] <<  std::endl;
//    std::cout << x_n[15]<< " " << x_n[16] << " " << x_n[17] <<  std::endl;
//



//    std::cout << "out_coil_state v_n: " << out_coil_state[0] << " " << out_coil_state[1] << " " << out_coil_state[2] <<  std::endl;
//    std::cout << "out_coil_state w_n: " << out_coil_state[3] << " " << out_coil_state[4] << " " << out_coil_state[5] <<  std::endl;
//    std::cout << "out_coil_state v dot: " << out_xdot_n[0] << " " << out_xdot_n[1] << " " << out_xdot_n[2] <<  std::endl;
//    std::cout << "out_coil_state w dot: " << out_xdot_n[3] << " " << out_xdot_n[4] << " " << out_xdot_n[5] <<  std::endl;
//    std::cout << " ---------  " <<  std::endl;
}



//[x_np1, xdot_n] = RK2_step(x_n, t_n, h, Integrand)
template <typename adType>
void RK2_coildyn(adType in_x_n[NUM_COIL_STATES], adType in_n[3], adType g[3],  adType actMass, adType actInertia[9], adType damping[6], adType in_B0[3], adType in_muhat[9], adType in_mL[3],
                 adType out_x_np1[NUM_COIL_STATES], adType out_xdot_n[6] ) {


    adType nL[3], B0[3], muhat[9], mL[3],  twist_n[6];       // from input
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


    //RK2_STEP_STEP1:
    CoilIntegrad(twist_n, nL, g, R_n, actMass, actInertia, damping,B0, muhat, mL, xdot_n);
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
    CoilIntegrad(x_n_p_k1o2, nL, g, R_np1half, actMass, actInertia, damping,B0, muhat, mL,  k2oh);

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
    adType twist_n[6], R_n[9], p_n[3], nL[3];				// variables used in analytical calculation
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
    CoilIntegrad(twist_n, nL, g, R_n, actMass, actInertia, damping,B0, muhat, mL, xdot_n);

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
    CoilIntegrad(x_np1_hat, nL, g, R_np1_hat, actMass, actInertia, damping,B0, muhat, mL, xdot_np1_hat);

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
void DYNNLEquation(adType in_x[], adType out_y[], NLEqnParams<adType> Params, adType out_u0[3], adType out_tau[NUM_ACT_SET*3]) {

    // output for time advance, not used in BVP, just placeholders
    adType m_L[NUM_ACT_SET][3], n_L[NUM_ACT_SET][3], n_0[3];
    // don't forget to scale parameters before passing to the CRMSolverIVP
    for (int j= 0; j < NUM_ACT_SET; ++j) {
        for (int i = 0; i < 3; i++) {
            m_L[j][i] = IVALUE_SCALE_M * in_x[i + j*6];
            n_L[j][i] = IVALUE_SCALE_N * in_x[i+ j*6 + 3];
        }
    }

//
//    for (int i = 0; i < NUM_ACT_SET; ++i) {
//        std::cout << "nL dyn: " << n_L[i][0] << " " << n_L[i][1] << " " << n_L[i][2] << std::endl;
//        std::cout << "ML dyn: " << m_L[i][0] << " " << m_L[i][1] << " " << m_L[i][2] << std::endl;
//    }


    for (int i = 0; i < 3; i++) {
        n_0[i] = Params.TipForce[i];  // for free-tip, this parameter is not given by the nonlinear equation solver, and hence, does not need to be scaled
    }

//    for (int i = 0; i < NUM_ACT_SET; ++i) {
//        std::cout << "m_L dyn: " << m_L[i][0] << " " << m_L[i][1] << " " <<m_L[i][2] << std::endl;
//        std::cout << "n_L dyn: " << n_L[i][0] << " " << n_L[i][1] << " " <<n_L[i][2] << std::endl;
//    }

    adType x_coil[NUM_ACT_SET][NUM_COIL_STATES], out_x_coil[NUM_ACT_SET][NUM_COIL_STATES];
    adType MagMoment[NUM_ACT_SET][3], muhat[NUM_ACT_SET][9], actMass[NUM_ACT_SET], actInertia[NUM_ACT_SET][9];

    double mu[3], muhattemp[9];


    for (int j = 0; j < NUM_ACT_SET; ++j) {
        for (int i = 0; i < 3; ++i) {
            x_coil[j][i] = Params.v_L_pre[j][i];
            x_coil[j][i+3] = Params.w_L_pre[j][i];
            x_coil[j][i+6] = Params.p_pre[j][i];
        }
        for (int i = 0; i < 9; ++i) {
            x_coil[j][i+9] = Params.R_pre[j][i];
        }

        actMass[j] = Params.actMass[j];
        for (int i = 0; i < 9; ++i) {
            actInertia[j][i] = Params.actInertia[j][i];
        }

        for (int i = 0; i <3 ; ++i) {
            MagMoment[j][i] = Params.MagMoment[j][i];
        }
        for (int i = 0; i < 3; ++i) {
            mu[i] =  MagMoment[j][i];
        }
        wHat(mu,muhattemp);
        for (int i = 0; i < 9; ++i) {
            muhat[j][i] = muhattemp[i];
        }
    }


    auto & K = Params.K;
    auto & Kinv = Params.Kinv;
    auto & ustar = Params.ustar;
    auto & SegBounds = Params.SegBounds;

    double net_mL[3], tau[3], K2invResidual[3], u_t[3], du[3];
//    double n_0[3] = {0,0,0}; //

    double tau_0[3] = {0.0,0.0,0.0};
    double p_t[3], R_t[9], u_tau[3], p_[3], R_[9];
    int fsegi, actno, actseg, actno_mn;
    double u_L[3], u_f[3], p_f[3], R_f[9];
    double out_xdot[6], residual[NUM_ACT_SET][6];

    double p_L[3], R_L[9];
    double v1[3], v2[3], v3[3], v_val[3];

    // root configurations for residual
    double p_d[3], R_d[9];
    for (int i = 0; i < 3; ++i) {
        p_d[i] = Params.xi[3+9+i];
    }
    for (int i = 0; i < 9; ++i) {
        R_d[i] = Params.xi[3+i];
    }

    double RigidSegmentLength;
    double net_nL[3];  // placeholder for force applied on the flexible segment

    for (int segi = NUM_SEGMENTS-1; segi >=0; --segi) { // starting from the last segment
        if ( segi%2 == 0 ) {
            fsegi = segi>>1;

            if(segi == NUM_SEGMENTS-1){ //last segment is flexible
                // assume free tip no torque tau_0 = [0,0,0]
                mMult_AB<3,3,1>( Kinv[fsegi], tau_0, K2invResidual );
                mAdd_AB<3,1>( ustar[fsegi], K2invResidual, u_t );

                for (int i = 0; i < 3; ++i) {
                    p_t[i] = Params.xf[3+9+i];
                }
                for (int i = 0; i < 9; ++i) {
                    R_t[i] = Params.xf[3+i];
                }
                CRMFlexible_IVP_Back ( segi, p_t, R_t, Params,u_t , n_0,
                                       u_tau, p_, R_);

                mSub_AB<3,1>( u_tau, ustar[fsegi], du);
                mMult_AB<3,3,1>( K[fsegi], du, tau ); // moment at upper side of the coil

                actno = fsegi - 1;

                mSub_AB<3,1>( m_L[actno], tau, net_mL);  // the net torque applied to the downwards coil

                //copy to the output for forward calculation
                for (int i = 0; i < 3; ++i) {
                    out_tau[i+actno*3] = tau[i];
                }

            }else{ // the non-free tip flexible segments with coils on top

                actno = fsegi;
                mMult_AB<3,3,1>( Kinv[fsegi], m_L[actno], K2invResidual );
                mAdd_AB<3,1>( ustar[fsegi], K2invResidual, u_L );

                actseg = segi + 1; //should be the upper actuator

                //calculate starting positios and roations given last coil states
                double RigidSegmentLength	=	SegBounds[actseg+1]-SegBounds[actseg];	// how far we need to move along the length of the rigid segment to reach the next flexible segment
                for (int i = 0; i < 9; ++i) {
                    R_L[i] = out_x_coil[actno][i+9];
                }
                for (int i = 0; i < 3; ++i) {
                    p_L[i] = out_x_coil[actno][i+6] - R_L[ i*3 +2 ]*RigidSegmentLength * 0.5;
                }

                // nL is the force applied to the flexible segment
                CRMFlexible_IVP_Back ( segi, p_L, R_L, Params,u_L , n_L[actno],
                                       u_f, p_f, R_f);

                // we now compute the terms needed for the actuator below
                actno = fsegi - 1;
//                std::cout << "u_f: " << u_f[0] << " " << u_f[1] << " " << u_f[2] <<  std::endl;
//                std::cout << "p_f: " << p_f[0] << " " << p_f[1] << " " << p_f[2] <<  std::endl;

                if(actno > -1){ //if there is a coil linked below, we need to calculate the net torque again
                    mSub_AB<3,1>( u_f, ustar[fsegi], du);
                    mMult_AB<3,3,1>( K[fsegi], du, tau ); // calculate moment at upper side of the coil

                    mSub_AB<3,1>( m_L[actno], tau, net_mL); // net torque applied at the downwards coil
                    //copy to the output for forward calculation

                    for (int i = 0; i < 3; ++i) {
                        out_tau[i+actno*3] = tau[i];
                    }
                }

            }
        }
        else{

            actno = (segi -1 )>>1;

//            std::cout << "n_L dyn: " << n_L[actno][0] << " " << n_L[actno][1] << " " << n_L[actno][2] << std::endl;
//            std::cout << "net_mL dyn: " << net_mL[0] << " " << net_mL[1] << " " << net_mL[2] << std::endl;

            if(actno< NUM_ACT_SET-1){
                mSub_AB<3,1>( n_L[actno], n_L[actno+1], net_nL);
            }else{
                mSub_AB<3,1>( n_L[actno], n_0, net_nL);
            }

            // nL is the force applied to the flexible segment, that is negated in the calculation
            CoilDynamics(x_coil[actno], net_nL, Params.g, actMass[actno], actInertia[actno],
                         Params.damping[actno], Params.DELTA_T, Params.B0,
                         muhat[actno], net_mL, out_x_coil[actno], out_xdot);

//            std::cout << "RL now: " << std::endl;
//            std::cout <<  out_x_coil[actno][9] << " " << out_x_coil[actno][10]<< " " << out_x_coil[actno][11] <<  std::endl;
//            std::cout <<  out_x_coil[actno][12]<< " " << out_x_coil[actno][13]<< " " << out_x_coil[actno][14] <<  std::endl;
//            std::cout << out_x_coil[actno][15]<< " " << out_x_coil[actno][16] << " " << out_x_coil[actno][17] <<  std::endl;


//            for (int i = 0; i < NUM_COIL_STATES; ++i) {
//                std::cout << "out_xdot" << out_xdot[i] << std::endl;
//            }
//            std::cout << "------------------------------" << std::endl;
            // compute the residual against the last flexible segment above
            if(actno<NUM_ACT_SET-1){
                 RigidSegmentLength	=	SegBounds[segi+1]-SegBounds[segi];	// how far we need to move along the length of the rigid segment to reach the next flexible segment

                //indices for the residual
                actno_mn = actno + 1;

                for (int i = 0; i < 9; ++i) {
                    R_L[i] = out_x_coil[actno][i+9];
                }
                for (int i = 0; i < 3; ++i) {
                    p_L[i] = out_x_coil[actno][i+6] + R_L[ i*3 +2 ]*RigidSegmentLength * 0.5;
                }
                for (int i = 0; i < 3; ++i) {
                    residual[actno_mn][i] = p_f[i] - p_L[i];
                }
                for (int i = 0; i < 3; ++i) {
                    v1[i] = R_f[i*3] - R_L[i*3];
                    v2[i] = R_f[1+ i*3] - R_L[1+ i*3];
                    v3[i] = R_f[2 + i*3]  - R_L[2 + i*3];
                }

                v_val[0] = vNormSq<3>(v1);
                v_val[1] = vNormSq<3>(v2);
                v_val[2] = vNormSq<3>(v3);
                for (int i = 0; i < 3; ++i) residual[actno_mn][i+3] = sqrt(v_val[i]);
            }

        }

    }

    mCopy_AB<3>(u_f, out_u0);

    // last segment
    actno = 0;
    for (int i = 0; i < 3; ++i) {
        residual[actno][i] = p_f[i] - p_d[i];
    }
    // calculate vector norm of the rotation matrices
    for (int i = 0; i < 3; ++i) {
        v1[i] = R_f[i*3] - R_d[i*3];
        v2[i] = R_f[1+ i*3] - R_d[1+ i*3];
        v3[i] = R_f[2 + i*3]  - R_d[2 + i*3];
    }

    v_val[0] = vNormSq<3>(v1);
    v_val[1] = vNormSq<3>(v2);
    v_val[2] = vNormSq<3>(v3);
    for (int i = 0; i < 3; ++i) residual[actno][i+3] = sqrt(v_val[i]);

//    for (int i = 0; i < NUM_ACT_SET; ++i) {
//        std::cout << "residual p: " << residual[i][0] << " " << residual[i][1] << " " <<residual[i][2] << std::endl;
//        std::cout << "residual R: " << residual[i][3] << " " << residual[i][4] << " " <<residual[i][5] << std::endl;
//    }
//
//    std::cout << " --------------------------------- " << std::endl;

    for (int i = 0; i < NUM_ACT_SET; ++i) {
        for (int j = 0; j < 3; ++j) {
            out_y[j + i*6] = RESIDUAL_SCALE_P * residual[i][j];
        }
        for (int j = 3; j < 6; ++j) {
            out_y[j + i*6] = RESIDUAL_SCALE_R * residual[i][j]; // 10000000000
        }
    }
//
//    std::cout << "out_y: " << out_y[0] << " " << out_y[1] << " " << out_y[2] <<  std::endl;
//    std::cout << "out_y: " << out_y[3] << " " << out_y[4] << " " << out_y[5] <<  std::endl;
//    std::cout << "out_y: " << out_y[6] << " " << out_y[7] << " " << out_y[8] <<  std::endl;
//    std::cout << "out_y: " << out_y[9] << " " << out_y[10] << " " << out_y[11] <<  std::endl;
}

void CRMIVP_DYN(	 CRMIVPCoreParams<double> CoreParams, const double in_u0[3], const double in_p0[3], const double in_R0[9],
                     const double in_mL[NUM_ACT_SET][3], const double in_nL[NUM_ACT_SET][3], const double in_tau[NUM_ACT_SET][3], const double in_ftip[3],
                     double out_coil_state[NUM_ACT_SET][NUM_COIL_STATES],  double out_u_new[3], double out_p_new[3], double out_R_new[9],
                     double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]) {

    double x_coil[NUM_ACT_SET][NUM_COIL_STATES];
    double MagMoment[NUM_ACT_SET][3], muhat[NUM_ACT_SET][9], actMass[NUM_ACT_SET], actInertia[NUM_ACT_SET][9];


//    for (int i = 0; i < NUM_ACT_SET; ++i) {
//        std::cout << "nL dyn: " << in_nL[i][0] << " " << in_nL[i][1] << " " << in_nL[i][2] << std::endl;
//        std::cout << "ML dyn: " << in_mL[i][0] << " " << in_mL[i][1] << " " << in_mL[i][2] << std::endl;
//    }
    double mu[3], muhattemp[9];
    for (int j = 0; j < NUM_ACT_SET; ++j) {
        for (int i = 0; i < 3; ++i) {
            x_coil[j][i] = CoreParams.v_L_pre[j][i];
            x_coil[j][i+3] = CoreParams.w_L_pre[j][i];
            x_coil[j][i+6] = CoreParams.p_pre[j][i];
        }
        for (int i = 0; i < 9; ++i) {
            x_coil[j][i+9] = CoreParams.R_pre[j][i];
        }

        actMass[j] = CoreParams.actMass[j];
        for (int i = 0; i < 9; ++i) {
            actInertia[j][i] = CoreParams.actInertia[j][i];
        }

        for (int i = 0; i <3 ; ++i) {
            MagMoment[j][i] = CoreParams.MagMoment[j][i];
        }
        for (int i = 0; i < 3; ++i) {
            mu[i] =  MagMoment[j][i];
        }
        wHat(mu,muhattemp);
        for (int i = 0; i < 9; ++i) {
            muhat[j][i] = muhattemp[i];
        }
    }

    double u_new[3], p_new[3], R_new[9];

    int   actno, fsegip1;	// actuator no, flexible segment before, flexible segment after
    double RscTB0[3], Tb[3], Residual[3], inertiaw[3];
    auto & K = CoreParams.K;
    auto & Kinv = CoreParams.Kinv;
    auto & ustar = CoreParams.ustar;
    auto & B0 = CoreParams.B0;
    auto & SegBounds = CoreParams.SegBounds;
    double w[3], w_dot[3], w_hat[9], w_inertia_w[3], inertia_wdot[3], out_xdot[6], w_terms[3], Tb_ml[3], K2invResidual[3];

    double u_0[3], p0[3], R0[9], n_L[NUM_ACT_SET][3], m_L[NUM_ACT_SET][3], tau[NUM_ACT_SET][3];

    mCopy_AB<NUM_ACT_SET, 3>(in_nL, n_L);
    mCopy_AB<NUM_ACT_SET, 3>(in_mL, m_L);
    mCopy_AB<NUM_ACT_SET, 3>(in_tau, tau);

    for (int i = 0; i < 3; ++i) {
        u_0[i] = in_u0[i];
        p0[i] = in_p0[i];
    }
    for (int i = 0; i < 9; ++i) {
        R0[i] = in_R0[i];
    }

    double net_nL[3], net_mL[3];

    double RigidSegmentLength;

    double n_f[3];

    double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3];
    double update_coil_state[NUM_COIL_STATES];
    //forward pass
    for (int SegmentIndex = 0; SegmentIndex < NUM_SEGMENTS; ++SegmentIndex) {
        if(SegmentIndex %2 == 0){

            actno = SegmentIndex >> 1;

            if(SegmentIndex == NUM_SEGMENTS -1){ // if this is free tip flexible segment
                for (int i = 0; i < 3; ++i) {
                    n_f[i] = in_ftip[i];
                }
//                n_f[0] = n_f[1] = n_f[2] = 0.0;
            }else{
                for (int i = 0; i < 3; ++i) {
                    n_f[i] =  n_L[actno][i];
                }
            }
            CRMFlexForward_pass (  SegmentIndex, p0, R0,  CoreParams, u_0 , n_f, u_new, p_new, R_new, p_atLocMarkers);
        }else{
            RigidSegmentLength=(SegBounds[SegmentIndex+1]-SegBounds[SegmentIndex]);

            actno=(SegmentIndex-1)>>1;		// actuator no
            fsegip1=actno+1;


            if(actno< NUM_ACT_SET-1){
                mSub_AB<3,1>( n_L[actno], n_L[actno+1], net_nL);
            }else{
                mSub_AB<3,1>( n_L[actno], in_ftip, net_nL);
            }

            mSub_AB<3,1>( m_L[actno], tau[actno], net_mL); // net torque applied at the downwards coil


            CoilDynamics(x_coil[actno], net_nL, CoreParams.g, actMass[actno], actInertia[actno], CoreParams.damping[actno], CoreParams.DELTA_T,
                         CoreParams.B0, muhat[actno], net_mL, update_coil_state, out_xdot);

            for (int i = 0; i < 3; ++i) {
                w[i] = update_coil_state[i+3];
                w_dot[i] = out_xdot[i+3];
            }
            wHat(w,w_hat);

            mMult_ATB<3,3,1>(&(update_coil_state[9]),B0,RscTB0);
            mMult_AB<3,3,1>(muhat[actno],RscTB0,Tb);

            mMult_AB<3,3,1>(actInertia[actno], w, inertiaw);
            mMult_AB<3,3,1>(w_hat, inertiaw, w_inertia_w);
            mMult_AB<3,3,1>(actInertia[actno], w_dot, inertia_wdot);
            mAdd_AB<3,1>(inertia_wdot, w_inertia_w, w_terms);

//            mSub_AB<3,1>( Tb, m_L[actno], Tb_ml);
//            mSub_AB<3,1>( Tb_ml, w_terms, Residual);

            mMult_AB<3,3,1>( Kinv[fsegip1], tau[actno], K2invResidual );
            mAdd_AB<3,1>( ustar[fsegip1], K2invResidual, u_0 );

            for (int i = 0; i < 3; ++i) {
                p0[i] = update_coil_state[i+6] + 0.5*RigidSegmentLength * update_coil_state[9+ i*3+2 ];
            }

            for (int i = 0; i < 9; ++i) {
                R0[i] = update_coil_state[i+9];
            }

            for (int i = 0; i < NUM_COIL_STATES; ++i) {
                out_coil_state[actno][i] = update_coil_state[i];
            }

        }
    }

    // return tip state
    for (int i = 0; i < 3; ++i) {
        out_u_new[i] = u_new[i];
        out_p_new[i] = p_new[i];
    }
    for (int i = 0; i < 9; ++i) {
        out_R_new[i] = R_new[i];
    }

    // Copy marker locations to the output
    if (!CoreParams.FinalValueOnly) {
        for (int i=0; i<NUM_LOCALIZATION_MARKERS; i++) {
            for (int j=0; j<3; j++) {
                out_p_atLocMarkers[i][j]=p_atLocMarkers[(NUM_LOCALIZATION_MARKERS-1)-i][j];
            }
        }
    }

}

template <typename adType>
void DynamicsBVP(	CRMShootingMethodParams<adType> in_Params, const double xf[NUM_STATES],
                              double in_mL_initialguess[NUM_ACT_SET][3], double in_nL_initialguess[NUM_ACT_SET][3], double in_ftip_initialguess[3],
                              adType out_u0[3], adType out_mL[NUM_ACT_SET][3], adType out_nL[NUM_ACT_SET][3], adType out_tau[NUM_ACT_SET][3], adType out_ftip[3], int& out_localmin) {

    ContactModeType ContactMode = in_Params.ContactMode;
    int NLEq_Dim;  // Dimension of the Nonlinear Equation to Solve
    NLEq_Dim = NUM_DYN_RESIDUAL;

    // Call CRMSolverIVP_Prep, to pre-process parameters
    adType x_0[NUM_STATES];
    for (int i = 0; i < NUM_STATES; i++) {
        if (i < 3) {
            x_0[i] = 0.0; //placeholder
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

    for (int i = 0; i < NUM_STATES; ++i) {
        DYNNLEParams.xf[i] = xf[i];
    }

    // Scale parameters and call the nonlinear equation solver
//    const double uscaleinv = 1.0 / IVALUE_SCALE_U;
    const double nscaleinv = 1.0 / IVALUE_SCALE_N;
    const double mscaleinv = 1.0 / IVALUE_SCALE_M;

//    const double fscaleinv = 1.0 / IVALUE_SCALE_F;
    auto* initialguessscaled = new double [NLEq_Dim];
    auto* returnedparamscaled = new double [NLEq_Dim];

    for (int j = 0; j < NUM_ACT_SET; ++j) {
        for (int i = 0; i < 3; i++) {
            initialguessscaled[i + j*6] = mscaleinv * in_mL_initialguess[j][i];
            initialguessscaled[i + j*6 +3] = nscaleinv * in_nL_initialguess[j][i];
        }
    }

    for (int i = 0; i < 3; ++i) {
        DYNNLEParams.TipForce[i] = in_Params.TipForce[i];  // if the catheter is not in contact, the tip force specified within in_Params needs to be used; this would not be scaled as it is not changed by the solver
    }

    int localmin = 0;

    double* x = new double[NLEq_Dim]; // we will create a new variable here and not use initial guess scaled since truss-region-dogleg algorithm uses the same variable for both input and output

    double* residual = new double[NLEq_Dim];
    int info;
    int lwa = (NLEq_Dim * (3 * NLEq_Dim + 13)) / 2;
    double tol = 0.00001;
    double* wa = new double [lwa];
    for (int i = 0; i < NLEq_Dim; i++) x[i] = initialguessscaled[i];

//    int REPS = 100;
//    // Get starting timepoint
//    auto start = high_resolution_clock::now();
//    for (int cnt = 0; cnt < REPS; cnt++)
//        DYNNLEquation(x, returnedparamscaled, DYNNLEParams, out_u0);
//    auto stop = high_resolution_clock::now();
//    auto duration = duration_cast<microseconds>(stop - start);
//    std::cout << std::endl << "Average time taken by DYNNLEquation Solution in " << REPS << " repetitions: " << duration.count() / REPS << " microseconds" << std::endl;

    double tau[NUM_ACT_SET*3];

#if defined( TRUSTREGION )
        TrustRegionDogleg_dyn(NLEq_Dim, x, residual, tol, info, wa, lwa, DYNNLEParams, out_u0, tau);
#else // undefined
    exit(1);
#endif

    localmin = (info == 1) ? 0 : (info - 1);
    for (int i = 0; i < NLEq_Dim; i++) returnedparamscaled[i] = x[i];
    delete[] residual;
    delete[] x;
    delete[] wa;

    for (int j = 0; j < NUM_ACT_SET; ++j) {
        for (int i = 0; i < 3; i++) {
            out_mL[j][i] = IVALUE_SCALE_M * returnedparamscaled[i + j * 6];
            out_nL[j][i] = IVALUE_SCALE_N * returnedparamscaled[i + j * 6 + 3];
        }

        for (int i = 0; i < 3; ++i) {
            out_tau[j][i] = tau[i+j*3];
        }
    }


//    for (int i = 0; i < NUM_ACT_SET; ++i) {
//        std::cout << "out_mL: " << out_mL[i][0] << " " << out_mL[i][1] << " " <<out_mL[i][2] << std::endl;
//        std::cout << "out_nL: " << out_nL[i][0] << " " << out_nL[i][1] << " " <<out_nL[i][2] << std::endl;
//        std::cout << "out_tau: " << out_tau[i][0] << " " << out_tau[i][1] << " " <<out_tau[i][2] << std::endl;
//    }

    for (int i = 0; i < 3; ++i) {
        out_ftip[i] = in_Params.TipForce[i];  // if it is free-tip, return the tip force specified within in_Params
    }

    out_localmin = localmin;
    delete[] initialguessscaled;
    delete[] returnedparamscaled;

}


template <typename adType>
void DYNSolverIVP(	CRMShootingMethodParams<adType> in_Params, adType in_u0[3],
                      adType in_mL[NUM_ACT_SET][3], adType in_nL[NUM_ACT_SET][3], double in_tau[NUM_ACT_SET][3], adType in_ftip[3],
                      bool in_FinalValueOnly,
                      adType out_x_N[NUM_STATES], adType out_coil_state[NUM_ACT_SET][NUM_COIL_STATES],
                      double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]){

    adType x_0[NUM_STATES];
    for (int i = 0; i < NUM_STATES; i++) {
        if (i < 3) x_0[i] = in_u0[i];
        else if (i < 12) x_0[i] = in_Params.R0[i - 3];
        else if (i < 15) x_0[i] = in_Params.p0[i - 12];
    }

    CRMIVPCoreParams<adType> CoreParams;
    adType u_0[3], n_L[NUM_ACT_SET][3], m_L[NUM_ACT_SET][3], tau[NUM_ACT_SET][3], ftip[3];
    // We need to pass u_0 as input argument as the values in CoreParams will be overriden with the values provided in the input arguments - functionality needed for solving Boundary Value Problems (BVP)
    // copy to local variable
    for (int i = 0; i < 3; i++) u_0[i] = in_u0[i];
    for (int i = 0; i < 3; i++) ftip[i] = in_ftip[i];

    for (int j = 0; j < NUM_ACT_SET; ++j) {
        for (int i = 0; i < 3; i++) n_L[j][i] = in_nL[j][i];
        for (int i = 0; i < 3; i++) m_L[j][i] = in_mL[j][i];
        for (int i = 0; i < 3; i++) tau[j][i] = in_tau[j][i];

    }
    CRMSolverIVP_Prep(x_0, in_Params.IntegrationStepSize,
                      in_Params.Li, in_Params.dlambdainv,
                      in_Params.SegEndLambdas, in_Params.LocMarkerLambdas,
                      in_Params.K, in_Params.Kinv, in_Params.ustar,
                      in_Params.MagMoment, in_Params.fcumlambda,
                      in_Params.B0, in_Params.g, in_Params.actMass, in_Params.actInertia, in_Params.damping, in_Params.DELTA_T,
                      in_Params.v_L_pre, in_Params.w_L_pre, in_Params.p_pre, in_Params.R_pre,
                      m_L, n_L, true, in_FinalValueOnly,CoreParams);

    double u_new[3], p_new[3], R_new[9];
    double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3];
    CRMIVP_DYN(	 CoreParams, u_0, in_Params.p0, in_Params.R0, m_L, n_L, tau,ftip, out_coil_state, u_new, p_new, R_new, p_atLocMarkers);

    for (int i = 0; i < NUM_STATES; ++i) {
        if(i<3){
            out_x_N[i] = u_new[i];
        }else if(i<3+9){
            out_x_N[i] = R_new[i-3];
        }else{
            out_x_N[i] = p_new[i-3-9];
        }
    }

    // Copy marker locations to the output
    if (!in_FinalValueOnly) {
        for (int i=0; i<NUM_LOCALIZATION_MARKERS; i++) {
            for (int j=0; j<3; j++) {
                out_p_atLocMarkers[i][j]=p_atLocMarkers[(NUM_LOCALIZATION_MARKERS-1)-i][j];
            }
        }
    }

}

