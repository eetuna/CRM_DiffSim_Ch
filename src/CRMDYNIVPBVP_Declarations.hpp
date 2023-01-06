#pragma once


#define COIL_TRUSTREGION

#define t_step 0.0005
#define NUM_DYN_RESIDUAL 6 // residual dim for m + n
#define NUM_COIL_STATES (6+3+9) // v[3], w[3], p[3] ,R[9]

#define IVALUE_SCALE_M  1.0 //
#define IVALUE_SCALE_N	1.0		// the variable used in Nonlinear Solver is multiplied with this scale to calculate ftip (tip force) that will be used in IVP
#define RESIDUAL_SCALE_R	1.0 //(10.0)			// the residual for tip position error coming out of the IVP will be multiplied with this scale to return to the Nonlinear Solver
#define NUM_RESIDUAL 3 // moment residual for CRM BVP solver

#define out_Dim 3 + 3 + 3 + 9 + 9					// vL + wL + pL + RL + u,m,n



template <typename adType>
struct CRMIVPCore_DYNParams : CRMIVPCoreParams<adType>{
    double 	g[3];		                                //  Gravity vector (in spatial coordinates)
    double v_L_pre[3];                                  // The linear velocity at the coil (L) at the previous time period
    double w_L_pre[3];                                  // The angular velocity at the coil (L) at the previous time period
    double actMass[NUM_ACT_SET];
    double actInertia[NUM_ACT_SET][9];

    double p_pre[3];
    double R_pre[9];
    double m_L[3];
    double n_L[3];
    double damping[6]; //v and w of the coil
    double DELTA_T;
};

template <typename adType>
struct CRMShootingMethod_DYNParams : CRMShootingMethodParams<adType> {
    double 	g[3];		                                //  Gravity vector (in spatial coordinates)
    double v_L_pre[3];                                  // The linear velocity at the coil (L) at the previous time period
    double w_L_pre[3];                                  // The angular velocity at the coil (L) at the previous time period
    double actMass[NUM_ACT_SET];
    double actInertia[NUM_ACT_SET][9];

    double p_pre[3];
    double R_pre[9];
    double damping[6]; //v and w of the coil
    double DELTA_T;
};



template <typename adType>
struct DYNNLEqnParams : CRMIVPCore_DYNParams<adType> {
    ContactModeType ContactMode;			// Enumerated type defining catheter contact mode.  ContactMode == FREE_TIP if the catheter is not in contact with a surface, FIXED_TIP if catheter tip is constrained to be at TipContraintPoint
    double			TipConstraintPoint[3];	// The spatial coordinates of the point where the catheter tip is constrained to be (used if ContactMode == FIXED_TIP)
    adType			TipForce[3];			// External point force (in spatial coordinates) applied at the tip of the catheter (\lambda = 0) (used if ContactMode == FREE_TIP)
    // terms used in the dynamics
    double	u0_initialguess[3];
    double	ftip_initialguess[3];
};

template <typename adType>
struct CRMDYNFKData {
    CRMCatheterModelParams* CathParams;		// physical parameters of the catheter
    CatheterConfiguration*	CathConfig;		// catheter configuration in space
    ContactModeType ContactMode;			// contact mode
    double	TipConstraintPoint[3];			// spatial coordinates of the point that the catheter tip is constrained to - only used if ContactMode == ContactModeType::FIXED_TIP
    adType	TipForce[3];					// force applied at the tip of the catheter in spatial coordinates - only used if ContactMode == ContactModeType::FREE_TIP
    double	u0_initialguess[3];				// initial guess for the curvature at the base of the catheter
    double	mL_initialguess[3];				// initial guess for the curvature at the base of the catheter
    double	nL_initialguess[3];				// initial guess for the curvature at the base of the catheter
    double	ftip_initialguess[3];			// initial guess for the tip force in spatial coordinates - only used if ContactMode == ContactModeType::FIXED_TIP
    double	IntegrationStepSize;			// length stepsize used in numerical integration performed as part of IVP calculations (part of BVP)
    bool	FinalValueOnly;					// Flag used to indicate if only final value is returned (true) or if Marker Locations are returned as well (false)
    double	(*ReportedMarkerPos)[NUM_LOCALIZATION_MARKERS][3];  // output which returns the calculated marker locations
    double actInertia[NUM_ACT_SET][9];
    double v_L_pre[3];
    double w_L_pre[3];
    double pL_pre[3];
    double RL_pre[9];
    double damping[6];
    double DELTA_T;
};


/*****************************************************************
 * * * BVP functions with moment & internal forces * * *
 *****************************************************************/


template <typename adType>
void CRMConstructShootingMethodParamSet(	CRMCatheterModelParams CathParams, CatheterConfiguration CathConfig,
                                            adType InsertionLength, adType ActuationCurrents[NUM_ACT_SET][3],
                                            ContactModeType ContactMode,
                                            double TipConstraintPoint[3], adType TipForce[3],
                                            double IntegrationStepSize, double in_ActInertia[NUM_ACT_SET][9],
                                            double in_v_L_pre[3], double in_w_L_pre[3], double in_p_pre[3], double in_R_pre[9], double in_damping[6], double in_DELTA_T,
                                            CRMShootingMethod_DYNParams<adType> &ShootingParams);

/**
 * override the original ELQ to account for internal moment and force
 * @tparam adType
 * @param in_x
 * @param out_y
 * @param Params
 */
template <typename adType>
void CRM_DYNNLEquation(adType in_x[], adType out_y[], const double m_L[3], const double n_L[3], bool LastLinkFlexible, const StateVector<adType>& xi_i, int ind_s, DYNNLEqnParams<adType> Params);


/*********************************************************
 * * * IVP functions with moment & internal forces * * *
 *********************************************************/

/**
 * RUn the IVP of CRM with dynamics parameters for getting initial p and R at coil
 * @tparam adType
 * @param in_Params
 * @param in_u0
 * @param in_ftip
 * @param in_FinalValueOnly
 * @param out_x_N
 * @param out_MomentResidual
 * @param Tbcoil
 * @param pcoil
 * @param Rcoil
 * @param out_p_atLocMarkers
 */
template <typename adType>
void CRMSolverIVP(	CRMShootingMethod_DYNParams<adType> in_Params,
                      adType in_u0[3], adType in_mL[3], adType in_nL[3],
                      adType in_ftip[3], bool in_FinalValueOnly,
                      adType out_x_N[NUM_STATES], adType out_MomentResidual[3],
                      adType Tbcoil[3], adType pcoil[3], adType Rcoil[9],
                      double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]) ;


template <typename adType>
void CRMSolverIVP_DYNPrep( adType in_x_0[NUM_STATES], double in_IntegrationStepSize,
                           adType in_Li, double in_dlambdainv,
                           double in_SegEndLambdas[NUM_SEGMENTS], const double in_LocMarkerLambdas[NUM_LOCALIZATION_MARKERS],
                           double in_K[NUM_FLEX_SEG][9], double in_Kinv[NUM_FLEX_SEG][9], double in_ustar[NUM_FLEX_SEG][3],
                           double in_CoilAlignmentTurnAreaMatrix[NUM_ACT_SET][9],
                           adType in_MagMoment[NUM_ACT_SET][3], double in_fcumlambda[NUM_FCUM_LAMBDA+1][3],
                           double in_B0[3], double in_g[3], const double in_actMass[NUM_ACT_SET], double in_actInertia[NUM_ACT_SET][9], double in_damping[6], double in_delta_t,
                           double in_v_L_pre[3], double in_w_L_pre[3], double in_p_pre[3], double in_R_pre[9], double in_mL[3], double in_nL[3],
                           bool in_FinalValueOnly, CRMIVPCore_DYNParams<adType> &out_CoreParams	) ;

/**
 * IVP solver with internal moment and force included, used in thr coil dynamics IVP and BVP functions
 * @tparam adType
 * @param in_params
 * @param in_u
 * @param in_ftip
 * @param out_x_N
 * @param out_MomentResidual
 * @param out_Tbcoil
 * @param out_pcoil
 * @param out_Rcoil
 * @param out_p_atLocMarkers
 */
template <typename adType>
void CRMSolverIVP_Core ( const CRMIVPCore_DYNParams<adType>& in_params,
                         adType in_u0[3], adType in_m[3], adType in_nL[3], adType in_ftip[3],
                         StateVector<adType>& out_x_N, adType out_MomentResidual[3],adType out_Tbcoil[3], adType out_pcoil[3], adType out_Rcoil[9],
                         double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3] );

/**
 * Override the original CRMSolverIVP_PropagateBCThroughRigidLink function with internal moment term added
 * @tparam adType
 * @param Residual_i
 * @param RigidSegmentLength
 * @param MagMoment
 * @param B0
 * @param m_L
 * @param ustar_ip1
 * @param Kinv_ip1
 * @param xf_i
 * @param xi_ip1
 * @param Residual_ip1
 * @param out_Tbcoil
 * @param out_pcoil
 * @param out_Rcoil
 */
template<typename adType>
void CRMSolverIVP_PropagateBCThroughRigidLink( const adType Residual_i[3],
                                               const adType RigidSegmentLength, const adType MagMoment[3],
                                               const double B0[3], const double ustar_ip1[3], const double Kinv_ip1[9],
                                               const StateVector<adType>& xf_i, StateVector<adType>& xi_ip1, adType out_Tbcoil[3], adType out_pcoil[3], adType out_Rcoil[9]);

/**
 *
 * @tparam adType
 * @param flex_residual
 * @param ind_s
 * @param SegBounds
 * @param MagMoment
 * @param B0
 * @param ustar
 * @param Kinv
 * @param xf_i
 * @param xi_o
 * @param out_Tbcoil
 * @param out_pcoil
 * @param out_Rcoil
 */
template<typename adType>
void CRMRigidLinkPass(const double flex_residual[3], int ind_s, const CRMIVPCore_DYNParams<adType>& in_params,
                      const StateVector<adType>& xf_i, StateVector<adType>& xi_o, adType out_Tbcoil[3], adType out_pcoil[3], adType out_Rcoil[9]  );

/**
 * Integrate through one flexible segment & calculate the tip moment residual
 * @tparam adType
 * @param xi_i
 * @param ind_s
 * @param SegBounds
 * @param SegSteps
 * @param InsertedLength
 * @param dlambdainv
 * @param K
 * @param Kinv
 * @param l_zero
 * @param ustar
 * @param fcumlambda
 * @param n_L
 * @param m_L
 * @param LastFlexibleLink
 * @param FinalValueOnly
 * @param LocMarkers
 * @param NextLocMarker
 * @param p_atLocMarkers
 * @param xf_o
 * @param MomentResidual
 */
template<typename adType>
void CRMSolverIVP_PropagateBCThroughFlexibleLink(  const StateVector<adType>& xi_i, const int ind_s, const double SegBounds[NUM_SEGMENTS+1], const int SegSteps[NUM_FLEX_SEG],
                                                   const double InsertedLength, const double dlambdainv, const double K[NUM_FLEX_SEG][9],const double Kinv[NUM_FLEX_SEG][9],
                                                   const double l_zero[3], const double ustar[NUM_FLEX_SEG][3],const double fcumlambda[NUM_FCUM_LAMBDA+1][3], const double n_L[3],
                                                   const double m_L[3], const bool LastFlexibleLink, const bool FinalValueOnly, const double LocMarkers[NUM_LOCALIZATION_MARKERS], int &NextLocMarker,
                                                   double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3],
                                                   StateVector<adType>& xf_o, double MomentResidual[3] );

/********************************************
 * * * Coil dynamics functions  * * *
 ********************************************/

/**
 * Dynamics forward kinematics function
 * @tparam adType
 * @param in_x
 * @param Params
 * @param localmin
 * @return
 */
template <typename adType>
adTypeVector<adType> CRM_DynamicsFK_FreeSpace(const adTypeVector<adType> in_x, CRMDYNFKData<adType> Params, int& localmin);

/**
 * Shooting of the flexible segment in the dynamics
 * @tparam adType
 * @param in_Params
 * @param out_u0
 * @param out_ftip
 * @param out_localmin
 */
template <typename adType>
void CRMShootingMethodBVP_DYN(	DYNNLEqnParams<adType> in_Params,  const double m_L[3], const double n_L[3], const StateVector<adType>& xi_i, int ind_s, bool LastLinkFlexible,
                                  adType out_u0[3], StateVector<adType> &xf,  adType MomentResidual[3], adType out_ftip[3], int& out_localmin);


/**
 * Integration function for coil integrad
 * @tparam adType
 * @param in_coil_state
 * @param in_n
 * @param g
 * @param actMass
 * @param actInertia
 * @param damping
 * @param DELTA_T
 * @param in_tau
 * @param out_coil_state
 */
template <typename adType>
void CoilDynamics( adType in_coil_state[NUM_COIL_STATES], adType in_n[3], adType g[3],
                   adType actMass, adType actInertia[9], adType damping[6], double DELTA_T, adType in_tau[3], adType out_coil_state[NUM_COIL_STATES]);


template <typename adType>
void CoilIntegrad(adType in_twist[6], adType in_n[3], adType g[3], adType R[9], adType actMass, adType actInertia[9], adType damping[6], adType in_tau[3], adType twistdot[6]);

template <typename adType>
void SE3_TimeSpace(adType in_R_n[9], adType in_p_n[3], adType h, adType in_twist_n[6], adType out_R_np1[9], adType out_p_np1[3]) ;

template <typename adType>
void RK2_coildyn(adType in_x_n[NUM_COIL_STATES], adType in_n[3], adType g[3],  adType actMass, adType actInertia[9], adType damping[6], adType in_tau[3],
                 adType out_x_np1[NUM_COIL_STATES], adType out_xdot_n[6] ) ;

template <typename adType>
void ABM4_coildyn(	adType in_x_n[NUM_COIL_STATES],adType in_xdot_nm1[6], adType in_xdot_nm2[6], adType in_xdot_nm3[6],
                      adType in_x_nm1[NUM_COIL_STATES], adType in_x_nm2[NUM_COIL_STATES], adType in_x_nm3[NUM_COIL_STATES],
                      adType in_n[3], adType g[3],  adType actMass, adType actInertia[9], adType damping[6], adType in_tau[3],
                      adType out_x_np1[NUM_COIL_STATES], adType out_xdot_n[6]);

/**
 * Dynamics equation used in the trust region function
 * @tparam adType
 * @param in_x
 * @param out_y
 * @param Params
 * @param out_u0 : returns the optimized initial curvature
 */
template <typename adType>
void DYNNLEquation(adType in_x[], adType out_y[], DYNNLEqnParams<adType> Params, adType out_u0[3]);

/**
 * BVP functions optimizes the u, m, and n for each flexible + coil segment
 * @tparam adType
 * @param in_Params
 * @param in_u0_initialguess
 * @param in_mL_initialguess
 * @param in_nL_initialguess
 * @param in_ftip_initialguess
 * @param out_u0
 * @param out_mL
 * @param out_nL
 * @param out_ftip
 * @param out_localmin
 */
template <typename adType>
void DynamicsBVP(CRMShootingMethod_DYNParams<adType> in_Params, double in_u0_initialguess[3],
                 double in_mL_initialguess[3], double in_nL_initialguess[3], double in_ftip_initialguess[3],
                 adType out_u0[3], adType out_mL[3], adType out_nL[3], adType out_ftip[3], int& out_localmin);


template <typename adType>
void DYNSolverIVP(	CRMShootingMethod_DYNParams<adType> in_Params, adType in_u0[3],
                      adType in_mL[3], adType in_nL[3], adType in_ftip[3],
                      bool in_FinalValueOnly,
                      adType out_x_N[NUM_STATES], adType out_coil_state[NUM_COIL_STATES],
                      double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]);
