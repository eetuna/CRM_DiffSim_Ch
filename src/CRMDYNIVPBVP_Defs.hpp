#pragma once

template <typename adType>
void CRMSolverIVP(	CRMShootingMethod_DYNParams<adType> in_Params,
                      adType in_u0[3], adType in_mL[3], adType in_nL[3],
                      adType in_ftip[3], bool in_FinalValueOnly,
                      adType out_x_N[NUM_STATES], adType out_MomentResidual[3],
                      adType Tbcoil[3], adType pcoil[3], adType Rcoil[9],
                      double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]) {

    adType x_0[NUM_STATES];
    for (int i = 0; i < 3; i++) x_0[i] = in_Params.p0[i];
    for (int i = 0; i < 9; i++) x_0[i + 3] = in_Params.R0[i];
    for (int i = 0; i < 3; i++) x_0[i + 3 + 9] = in_u0[i];

    CRMIVPCore_DYNParams<adType> CoreParams;
    adType u_0[3], ftip[3];

    StateVector<adType> x_N;

    adType m_L[3] = {0.0, 0.0, 0.0};
    adType n_L[3] = {0.0, 0.0, 0.0};

    CRMSolverIVP_DYNPrep(x_0, in_Params.IntegrationStepSize,
                      in_Params.Li, in_Params.dlambdainv,
                      in_Params.SegEndLambdas, in_Params.LocMarkerLambdas,
                      in_Params.K, in_Params.Kinv, in_Params.ustar,
                      in_Params.CoilAlignmentTurnAreaMatrix,
                      in_Params.MagMoment, in_Params.fcumlambda,
                      in_Params.B0,in_Params.g, in_Params.actMass, in_Params.actInertia, in_Params.damping, in_Params.DELTA_T,
                      in_Params.v_L_pre, in_Params.w_L_pre, in_Params.p_pre, in_Params.R_pre,
                      m_L, n_L, in_FinalValueOnly,CoreParams);


    // copy to local variable
    for (int i = 0; i < 3; i++) ftip[i] = in_ftip[i];

    // We need to pass u_0 as input argument as the values in CoreParams will be overriden with the values provided in the input arguments - functionality needed for solving Boundary Value Problems (BVP)
//    for (int i = 0; i < 3; i++) u_0[i] = in_u0[i];

    CRMSolverIVP_Core(CoreParams, in_u0, in_mL, in_nL, ftip, x_N, out_MomentResidual, Tbcoil, pcoil, Rcoil, out_p_atLocMarkers);

    // NO NEED FOR CRMSolverIVP_Return
    mCopy_AB<3>(x_N._p, out_x_N + 0);
    mCopy_AB<9>(x_N._R, out_x_N + 3);
    mCopy_AB<3>(x_N._u, out_x_N + 3 + 9);

}


template <typename adType>
void CRMSolverIVP_DYNPrep ( adType in_x_0[NUM_STATES], double in_IntegrationStepSize,
                         adType in_Li, double in_dlambdainv,
                         double in_SegEndLambdas[NUM_SEGMENTS], const double in_LocMarkerLambdas[NUM_LOCALIZATION_MARKERS],
                         double in_K[NUM_FLEX_SEG][9], double in_Kinv[NUM_FLEX_SEG][9], double in_ustar[NUM_FLEX_SEG][3],
                         double in_CoilAlignmentTurnAreaMatrix[NUM_ACT_SET][9],
                         adType in_MagMoment[NUM_ACT_SET][3], double in_fcumlambda[NUM_FCUM_LAMBDA+1][3],
                         double in_B0[3], double in_g[3], const double in_actMass[NUM_ACT_SET], double in_actInertia[NUM_ACT_SET][9], double in_damping[6], double in_delta_t,
                         double in_v_L_pre[3], double in_w_L_pre[3], double in_p_pre[3], double in_R_pre[9], double in_mL[3], double in_nL[3],
                         bool in_FinalValueOnly, CRMIVPCore_DYNParams<adType> &out_CoreParams	) {
    // Process the incoming parameters (including changing from distal-proximal order to proximal-distal order)
    //   and package them to be passed to CRMSolverIVP_Core
    // Everything is copied to local variables (inside out_CoreParams), therefore, this function
    //   can be used to transfer data from global memory to device memory for subsequent computations
    // The packaged data can be used to call CRMSolverIVP_Core multiple times by changing
    //   only u[0..2] components of xi --- other parameters should not change

    // local copies for variables accessed out-of-order
    double SegEndLambdas[NUM_SEGMENTS];
    mCopy_AB<NUM_SEGMENTS>(in_SegEndLambdas,SegEndLambdas);

    // create aliases for variables in CoreParams
    auto & xi = out_CoreParams.xi;

    auto & SegBounds = out_CoreParams.SegBounds;
    auto & SegSteps = out_CoreParams.SegSteps;

    auto & InsertedLength = out_CoreParams.InsertedLength;
    auto & dlambdainv = out_CoreParams.dlambdainv;
    auto & K = out_CoreParams.K;
    auto & Kinv = out_CoreParams.Kinv;
    auto & ustar = out_CoreParams.ustar;
    auto & CoilAlignmentTurnAreaMatrix = out_CoreParams.CoilAlignmentTurnAreaMatrix;
    auto & fcumlambda = out_CoreParams.fcumlambda;
//	auto & ftip = out_CoreParams.ftip;
    auto & FinalValueOnly = out_CoreParams.FinalValueOnly;
    auto & LocMarkers = out_CoreParams.LocMarkers;

    auto & B0 = out_CoreParams.B0;
    auto & MagMoment = out_CoreParams.MagMoment;

    auto & StartSegmentIndex = out_CoreParams.StartSegmentIndex;
    auto & NextLocMarker = out_CoreParams.NextLocMarker;
    auto & p_atLocMarkers = out_CoreParams.p_atLocMarkers;

    // local variables
    double IntegrationStepSize=in_IntegrationStepSize;
    double DeltaSInv=1.0/IntegrationStepSize;
    adType tempadType;

    // process parameters as needed and copy into CoreParams
    mCopy_AB<NUM_STATES>(in_x_0,xi);		// Initial value of the state for the next segment to be integrated
    //  States are packed p[0..2],R[0..8],u[0..2]  (R: 3x3 matrix stored in row major order R11 R12 R13 R21 R22 R23 R31 R32 R33)
    mCopy_AB<3>(in_B0,B0);					//  B0 field vector of the MRI scanner (in spatial coordinates)

    dlambdainv=in_dlambdainv;				// reciprocal of dlambda (lambda stepsize used in discretizing fcumlambda)
    FinalValueOnly=in_FinalValueOnly;		// Flag used to indicate if only final value (xf) is returned (true) or if Marker Locations are returned as well (false)
    InsertedLength=in_Li;					// Inserted Length (length of the catheter from the entry point to the tip)
    // If the catheter is inserted more than the length of the catheter, clamp it to catheter length
    if (InsertedLength > SegEndLambdas[NUM_SEGMENTS - 1]) InsertedLength = SegEndLambdas[NUM_SEGMENTS - 1];

    // WE ARE GOING TO REORDER SEGMENT AND ACTUATOR UNITS SO THAT THEY ARE ORDERED FROM THE INSERTION POINT TO THE TIP
    //   I.E., SWITCH TO PROXIMAL TO DISTAL ORDERING

    StartSegmentIndex=NUM_SEGMENTS-1;		// Index of the segment where the integration to solve IVP will start -- the segment located at the entry point; note that segment indices start at 0
    SegBounds[NUM_SEGMENTS]=InsertedLength;	// End s value of last segment is s=InsertedLength
    // Entry point has a value of s=0 -- we may not simulate full length of the most proximal segment in the chamber
    for (int i=0; i<NUM_SEGMENTS; i++) {
        tempadType =InsertedLength-SegEndLambdas[i];
        if (tempadType >0.0) {
            SegBounds[(NUM_SEGMENTS-1)-i]= tempadType;
            StartSegmentIndex--;				// Integration will start at the previous segment
        }
        else {
            SegBounds[(NUM_SEGMENTS-1)-i]=0.0;	// This segment boundary is still inside the sheath
        }
    }

    for (int i=0; i<NUM_FLEX_SEG; i++) { 		// flexible catheter segment
        // Calculate the number of integration steps based on the given IntegrationStepSize
        SegSteps[i]=int(ceil( dVal(SegBounds[2*i+1]-SegBounds[2*i])*DeltaSInv ));
    }

    NextLocMarker=NUM_LOCALIZATION_MARKERS;
    if (!FinalValueOnly) {
        // While changing the order and converting from lambda to s,
        //   also find the index of the first localization marker after the entry point
        //   and assign (extrapolated) locations to markers which are still inside the sheath
        for (int i=0; i<NUM_LOCALIZATION_MARKERS; i++) {
            tempadType =InsertedLength-in_LocMarkerLambdas[i];
            LocMarkers[(NUM_LOCALIZATION_MARKERS-1)-i]= dVal(tempadType);
            if (tempadType >0.0) NextLocMarker--;
            else {   // and assign (extrapolated) locations to markers which are still inside the sheath
                LocMarkerUpdate(p_atLocMarkers[(NUM_LOCALIZATION_MARKERS - 1) - i], xi, tempadType);
            }
        }
    }

    // be careful - order is reversed in K, Kinv, ustar, CoilAlignmentTurnAreaMatrix, and MagMoment
    for (int i=0; i<NUM_FLEX_SEG; i++) {
        for (int j=0; j<9; j++) {
            K[(NUM_FLEX_SEG-1)-i][j]		=	in_K[i][j];
            Kinv[(NUM_FLEX_SEG-1)-i][j]		=	in_Kinv[i][j];
        }
        for (int j=0; j<3; j++) {
            ustar[(NUM_FLEX_SEG-1)-i][j]	=	in_ustar[i][j];
        }
    }
    for (int i=0; i<NUM_ACT_SET; i++) {
        for (int j=0; j<3; j++) {
            MagMoment[(NUM_ACT_SET-1)-i][j]	=	in_MagMoment[i][j];
        }
        mCopy_AB<9>(in_CoilAlignmentTurnAreaMatrix[i], CoilAlignmentTurnAreaMatrix[(NUM_ACT_SET - 1) - i]);
    }

    mCopy_ABm<(NUM_FCUM_LAMBDA+1),3>(in_fcumlambda,fcumlambda);
    //mCopy_AB<3>(in_ftip, ftip);

    /**
     * The params for dynamics
     * **/
    auto & g = out_CoreParams.g;
    mCopy_AB<3>(in_g, g);               // gravitational vector
    auto & v_L_pre = out_CoreParams.v_L_pre;
    auto & w_L_pre = out_CoreParams.w_L_pre;
    mCopy_AB<3>(in_v_L_pre, v_L_pre);               // linear velocity vector
    mCopy_AB<3>(in_w_L_pre, w_L_pre);               // angular velocity vector

    auto & p_pre = out_CoreParams.p_pre;
    auto & R_pre = out_CoreParams.R_pre;
    mCopy_AB<3>(in_p_pre, p_pre);
    mCopy_AB<9>(in_R_pre, R_pre);

    auto & actMass = out_CoreParams.actMass;
    auto & actInertia = out_CoreParams.actInertia;
    for (int i = 0; i < NUM_ACT_SET; ++i) {
        actMass[i] = in_actMass[i];
        for (int j = 0; j < 9; ++j) {
            actInertia[i][j] = in_actInertia[i][j];
        }
    }
    auto & m_L = out_CoreParams.m_L;
    mCopy_AB<3>(in_mL, m_L);

    auto & n_L = out_CoreParams.n_L;
    mCopy_AB<3>(in_nL, n_L);

    auto & damping = out_CoreParams.damping;
    mCopy_AB<6>(in_damping, damping);

    auto & delta_t = out_CoreParams.DELTA_T;
    delta_t = in_delta_t;

}


template <typename adType>
void CRMSolverIVP_Core ( const CRMIVPCore_DYNParams<adType>& in_params,
                         adType in_u[3], adType in_mL[3], adType in_nL[3],  adType in_ftip[3],
                         StateVector<adType>& out_x_N, adType out_MomentResidual[3],adType out_Tbcoil[3], adType out_pcoil[3], adType out_Rcoil[9],
                         double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3] ){

//    /**
//     * Parameters for dynamics
//     */
//    auto & m_L = in_params.m_L;
//    auto & n_L = in_params.n_L;

//    // convenience definitions
//    const double Identity3x3[9] = { 1,0,0,0,1,0,0,0,1 };
//    const double Zero3[3] = { 0,0,0 };

    double l_zero[3] = { 0.0,0.0,0.0 };
    adType RigidSegmentLength;		// Length of the rigid segment - intermediate variable
    //adType muhat[9];
    //adType RscTB0[3], Tb[3], deltau1[3], K1deltau1[3], K2invResidual[3]; // intermediate variables
    int   actno, fsegi, fsegip1;	// actuator no, flexible segment before, flexible segment after
    bool  LastSegmentIsRigid = true;	// Flag indicating if the last segment processed is rigid (true) or not (false)

    // we need to copy in_ftip to local variable
    adType ftip[3];
    mCopy_AB<3>(in_ftip, ftip);

    // we will copy anything we will access more than once (or write to) to local variables
    int   StartSegmentIndex = in_params.StartSegmentIndex;
    int	  NextLocMarker = in_params.NextLocMarker;
    int	  InitialLocMarker = NextLocMarker;
    // for others, we will create aliases
    auto& SegBounds = in_params.SegBounds;
    auto& SegSteps = in_params.SegSteps;
    auto& InsertedLength = in_params.InsertedLength;
    auto& dlambdainv = in_params.dlambdainv;
    auto& K = in_params.K;
    auto& Kinv = in_params.Kinv;
    auto& ustar = in_params.ustar;
    auto& fcumlambda = in_params.fcumlambda;
    auto& FinalValueOnly = in_params.FinalValueOnly;
    auto& LocMarkers = in_params.LocMarkers;
    auto& B0 = in_params.B0;
    auto& MagMoment = in_params.MagMoment;
    auto& CoilAlignmentTurnAreaMatrix = in_params.CoilAlignmentTurnAreaMatrix;

    StateVector<adType> xi;  		//  Initial value of the state for the next segment to be integrated
    auto& xf = out_x_N;					//  Final value of the state for the last segment integrated
    auto& Residual = out_MomentResidual;	// Residual at the catheter tip -- will be returned
    double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3];	// Positions at markers (ordered proximal to distal)

    // we will also copy the already filled entries of p_atLocMarkers from params
    for (int i = 0; i < InitialLocMarker; i++) {
        for (int j = 0; j < 3; j++) {
            p_atLocMarkers[i][j] = in_params.p_atLocMarkers[i][j];
        }
    }

    // initial conditions for the regular states
    // we need to copy xi from in_params to the local variable and update it with u[0..2] specified in in_u
    mCopy_AB<3>(in_params.xi + 0, xi._p);
    mCopy_AB<9>(in_params.xi + 3, xi._R);
    mCopy_AB<3>(in_u, xi._u);


    // IMPORTANT NOTE: most proximal segment is assumed to be always flexible
    //    and the flexible and rigid segments are assumed to be alternating
    //    most distal segment can be flexible or rigid

    // If the starting segment is a rigid segment, then we will need to move initial conditions to the start of the next flexible segment and change the StartingSegment to that segment
    if ( StartSegmentIndex % 2 == 1 ) {
        RigidSegmentLength	=	SegBounds[StartSegmentIndex+1]-SegBounds[StartSegmentIndex];	// how far we need to move along the length of the rigid segment to reach the next flexible segment
        for (int i = 0; i < 3; i++) { 										// u[0..2] and R[0..8] remain the same
            xi._p[i] = xi._p[i] + RigidSegmentLength * xi._R[i * 3 + 2];	// p[0..2] will translate along the z direction of the R matrix (3rd column)
        }
        StartSegmentIndex++;  // move the start segment to the next segment
        if (!FinalValueOnly) {
            // are there any localization markers?  If so, calculate their positions
            CalculateLocMarkers(NextLocMarker, xf, LocMarkers, SegBounds[StartSegmentIndex], p_atLocMarkers);  //	xf has already been updated
        }
    }

    // In case there is nothing to integrate or the rigid segment is the only segment (it is both the start and the end segment),
    //   in which case, the main integration loop will not execute,
    //   we need to have valid return values for xf and
    //   Residual is initialized to the 0 vector and residual is initialized to the 0 vector
    xf = xi;
    for (int i = 0; i < 3; i++) Residual[i] = 0.0;

    adType cur_residual[3];
    for (int i = 0; i < 3; i++) cur_residual[i] = 0.0;

    bool lastFlexibleLink = false;
    // Integrate each of the remaining segments
    for (int i=StartSegmentIndex; i<NUM_SEGMENTS; i++){
            if ( i%2 == 0 ) {  // Flexible Segment
            LastSegmentIsRigid=false;
             int   fsegno=i>>1; // i/2, flexible segment no

             CRMSolverIVP_PropagateBCThroughFlexibleLink(  xi, i, SegBounds, SegSteps, InsertedLength, dlambdainv, K, Kinv,
                                                               l_zero, ustar, fcumlambda, in_nL, in_mL,
                                                               lastFlexibleLink, FinalValueOnly, LocMarkers, NextLocMarker,
                                                               p_atLocMarkers,xf, cur_residual);

        }else {  // Need to do actuation/rigid segment calculations to transfer Initial Conditions to next flexible segment

            mCopy_AB<3>(cur_residual, Residual);
            lastFlexibleLink = true;

            CRMRigidLinkPass(cur_residual, i, in_params, xi, xf,
                                 out_Tbcoil, out_pcoil, out_Rcoil);

            if (!FinalValueOnly) {
                // are there any localization markers?  If so, calculate their positions
                CalculateLocMarkers(NextLocMarker, xf, LocMarkers, SegBounds[i + 1], p_atLocMarkers);
            }

        }
        xi = xf; // The calculated final values will be the initial value of the next iteration
    }

    // Copy marker locations to the output
    //    note that the order is being reversed
    if (!FinalValueOnly) {
        for (int i=0; i<NUM_LOCALIZATION_MARKERS; i++) {
            for (int j=0; j<3; j++) {
                out_p_atLocMarkers[i][j]=p_atLocMarkers[(NUM_LOCALIZATION_MARKERS-1)-i][j];
            }
        }
    }

    // no need to copy the final values of the state x_f and residual to the output
    // since we have created an alias
    //out_x_N = xf;
    //mCopy_AB<3>(Residual, out_MomentResidual);

}


template<typename adType>
void CRMSolverIVP_PropagateBCThroughRigidLink( const adType Residual_i[3],
                                               const adType RigidSegmentLength, const adType MagMoment[3],
                                               const double B0[3], const double ustar_ip1[3], const double Kinv_ip1[9],
                                               const StateVector<adType>& xf_i, StateVector<adType>& xi_ip1, adType out_Tbcoil[3], adType out_pcoil[3], adType out_Rcoil[9]) {

    adType muhat[9], RscTB0[3], Tb[3], K2invResidual[3]; // intermediate variables

    // p
    for (int j = 0; j < 3; j++) {
        xi_ip1._p[j] = xf_i._p[j] + xf_i._R[j * 3 + 2] * RigidSegmentLength;
    }
    // R
    mCopy_AB<9>(xf_i._R, xi_ip1._R);
    // u
    // Tb=\mu_c \cross R_sc^T B0,s
    wHat<adType>(MagMoment, muhat);
    mMult_ATB<3, 3, 1>(xf_i._R, B0, RscTB0);
    mMult_AB<3, 3, 1>(muhat, RscTB0, Tb);
    // u2=u2star + ( K2inv K1 (u1 - u1star ) - K2inv Tb )

//    mSub_AB<3,1>( m_L, Residual_i , Residual_ip1);

    for (int j = 0; j < 3; ++j) {
        out_Tbcoil[j] = Tb[j];
    }

    for (int j=0; j<3; j++) {
        out_pcoil[j]=xf_i._p[j] + xf_i._R[j * 3 + 2] * RigidSegmentLength * 0.5;
    }
    mCopy_AB<9>(xf_i._R, out_Rcoil);

    mMult_AB<3, 3, 1>(Kinv_ip1, Residual_i, K2invResidual);
    mAdd_AB<3, 1>(ustar_ip1, K2invResidual, xi_ip1._u);

}


template<typename adType>
void CRMSolverIVP_PropagateBCThroughFlexibleLink(  const StateVector<adType>& xi_i, const int ind_s, const double SegBounds[NUM_SEGMENTS+1], const int SegSteps[NUM_FLEX_SEG],
                                                   const double InsertedLength, const double dlambdainv, const double K[NUM_FLEX_SEG][9],const double Kinv[NUM_FLEX_SEG][9],
                                                   const double l_zero[3], const double ustar[NUM_FLEX_SEG][3],const double fcumlambda[NUM_FCUM_LAMBDA+1][3], const double n_L[3],
                                                   const double m_L[3], const bool LastFlexibleLink, const bool FinalValueOnly, const double LocMarkers[NUM_LOCALIZATION_MARKERS], int &NextLocMarker,
                                                   double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3],
                                                   StateVector<adType>& xf_o, double MomentResidual[3] ){
    // Prepare the CRMIntegrand Parameters
    int   fsegno; 					// flexible segment no
    adType h;						// integration stepsize
    adType deltau[3];

    fsegno=ind_s>>1; // i/2, flexible segment no

    // Calculate the actual stepsize, based on the number of steps
    h=(SegBounds[ind_s+1]-SegBounds[ind_s])/(SegSteps[fsegno]*1.0);


    // Integrate
    ABM4( 	xi_i, SegBounds[ind_s], SegSteps[fsegno], h,
             InsertedLength, dlambdainv, K[fsegno], Kinv[fsegno], l_zero, ustar[fsegno],
             fcumlambda, n_L,
             FinalValueOnly, LocMarkers, NextLocMarker,
             xf_o, p_atLocMarkers);

    if(!LastFlexibleLink){
        double freetip_residual[3];
        mSub_AB<3, 1>(xf_o._u, ustar[fsegno], deltau);
        mMult_AB<3, 3, 1>(K[fsegno], deltau, freetip_residual);

        mSub_AB<3,1>( m_L, freetip_residual , MomentResidual);
    }else{
        // residual = K (u1 - u1star)
        mSub_AB<3, 1>(xf_o._u, ustar[fsegno], deltau);
        mMult_AB<3, 3, 1>(K[fsegno], deltau, MomentResidual);
    }

}

template<typename adType>
void CRMRigidLinkPass(const double flex_residual[3], int ind_s, const CRMIVPCore_DYNParams<adType>& in_params,
                      const StateVector<adType>& xf_i, StateVector<adType>& xi_o, adType out_Tbcoil[3], adType out_pcoil[3], adType out_Rcoil[9]  ){
    int   actno, fsegip1;	// actuator no, flexible segment before, flexible segment after
int  fsegi;
    double RigidSegmentLength=(in_params.SegBounds[ind_s+1]-in_params.SegBounds[ind_s]);
    actno = (ind_s - 1) >> 1;		// actuator no
    fsegi = actno;				// index of flexible segment before (immediately proximal to) the rigid link
    fsegip1 = actno + 1;		// index of flexible segment after (immediately distal to) the rigid link

    CRMSolverIVP_PropagateBCThroughRigidLink(flex_residual, RigidSegmentLength, in_params.MagMoment[actno], in_params.B0,
                                             in_params.ustar[fsegip1], in_params.Kinv[fsegip1], xf_i, xi_o,
                                             out_Tbcoil, out_pcoil, out_Rcoil);


};



template <typename adType>
void CRMConstructShootingMethodParamSet(	CRMCatheterModelParams CathParams, CatheterConfiguration CathConfig,
                                            adType InsertionLength, adType ActuationCurrents[NUM_ACT_SET][3],
                                            ContactModeType ContactMode,
                                            double TipConstraintPoint[3], adType TipForce[3],
                                            double IntegrationStepSize, double in_ActInertia[NUM_ACT_SET][9],
                                            double in_v_L_pre[3], double in_w_L_pre[3], double in_p_pre[3], double in_R_pre[9], double in_damping[6], double in_DELTA_T,
                                            CRMShootingMethod_DYNParams<adType> &ShootingParams) {
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

        ShootingParams.actMass[i] = CathParams.ActMass[i];
        for (int j = 0; j < 9; ++j) {
            ShootingParams.actInertia[i][j] = in_ActInertia[i][j];
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
    ShootingParams.IntegrationStepSize = IntegrationStepSize;
    mCopy_AB<9>(CathConfig.R0, ShootingParams.R0);
    mCopy_AB<3>(CathConfig.p0, ShootingParams.p0);

    ShootingParams.ContactMode = ContactMode;
    mCopy_AB<3>(TipConstraintPoint, ShootingParams.TipConstraintPoint);
    mCopy_AB<3>(TipForce, ShootingParams.TipForce);

    mCopy_AB<3>(CathConfig.g, ShootingParams.g);
    mCopy_AB<3>(in_v_L_pre, ShootingParams.v_L_pre);
    mCopy_AB<3>(in_w_L_pre, ShootingParams.w_L_pre);
    mCopy_AB<3>(in_p_pre, ShootingParams.p_pre);
    mCopy_AB<9>(in_R_pre, ShootingParams.R_pre);
    mCopy_AB<6>(in_damping, ShootingParams.damping);
    ShootingParams.DELTA_T = in_DELTA_T;
}

template <typename adType>
void CRM_DYNNLEquation(adType in_x[], adType out_y[], const double m_L[3], const double n_L[3], bool LastLinkFlexible, const StateVector<adType>& xi_i, int ind_s, DYNNLEqnParams<adType> Params) {

    StateVector<adType> x_N;
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
//    // We will only call the IVP_Core, since preprocessing is already done
//    double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3];
//    adType Tbcoil[3], pcoil[3], Rcoil[9];
//    CRMSolverIVP_Core(Params, u_0, ftip, x_N, MomentResidual, Tbcoil, pcoil, Rcoil, p_atLocMarkers);
//
////    std::cout <<"CRMSolverIVP_Core residual: " << MomentResidual[0] << " "  <<MomentResidual[1] <<  " " << MomentResidual[2] << std::endl;

    double l_zero[3]= {0,0,0};

    //u_0 curvature is the only one get updated
    mCopy_AB<3>(u_0, xi_i._u);

    CRMSolverIVP_PropagateBCThroughFlexibleLink( xi_i, ind_s, Params.SegBounds, Params.SegSteps, Params.InsertedLength,
            Params.dlambdainv, Params.K, Params.Kinv, l_zero, Params.ustar, Params.fcumlambda, n_L,
            m_L, LastLinkFlexible, Params.FinalValueOnly, Params.LocMarkers, Params.NextLocMarker,
            Params.p_atLocMarkers, x_N, MomentResidual );

//    std::cout <<"flexible residual: " << MomentResidual[0] << " "  <<MomentResidual[1] <<  " " << MomentResidual[2] << std::endl;

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

template <typename adType>
adTypeVector<adType> CRM_DynamicsFK_FreeSpace(const adTypeVector<adType> in_x, CRMDYNFKData<adType> Params, int& localmin) {

//    int X_Dim = NUM_ACT_SET * 3 + 1;	// Dimension of the Input 3 x NUM_ACT_SET + 1 for inserted length

    adTypeVector<adType> out_y(out_Dim);

    adType ActuationCurrents[NUM_ACT_SET][3];
    adType InsertedLength;
    adType ftip_calc[3];
    adType xf[NUM_STATES];
    CRMShootingMethod_DYNParams<adType> BVPParams;
    auto& CathParams = *(Params.CathParams);
    auto& CathConfig = *(Params.CathConfig);

    for (int i = 0; i < NUM_ACT_SET; i++)	for (int j = 0; j < 3; j++)	ActuationCurrents[i][j] = in_x(i * 3 + j);
    InsertedLength = in_x(NUM_ACT_SET * 3);


    CRMConstructShootingMethodParamSet(CathParams, CathConfig, InsertedLength, ActuationCurrents,Params.ContactMode,
                                       Params.TipConstraintPoint, Params.TipForce, Params.IntegrationStepSize, Params.actInertia,
                                       Params.v_L_pre, Params.w_L_pre,  Params.pL_pre, Params.RL_pre,
                                       Params.damping, Params.DELTA_T, BVPParams);

    double out_u0[3], out_nL[3], out_mL[3];
    DynamicsBVP(BVPParams, Params.u0_initialguess, Params.mL_initialguess, Params.nL_initialguess, Params.ftip_initialguess,
                out_u0, out_mL, out_nL, ftip_calc, localmin);

//    std::cout << "out_u0: " << out_u0[0] << " " << out_u0[1] << " " << out_u0[2] <<  std::endl;
//    std::cout << "out_mL: " << out_mL[0] << " " << out_mL[1] << " " << out_mL[2] <<  std::endl;
//    std::cout << "out_nL: " << out_nL[0] << " " << out_nL[1] << " " << out_nL[2] <<  std::endl;

    double x_coil[NUM_COIL_STATES];
    DYNSolverIVP(BVPParams, out_u0, out_mL, out_nL, ftip_calc,
                 true, xf, x_coil,*(Params.ReportedMarkerPos));

    for (int i = 0; i < NUM_COIL_STATES; ++i) {
        out_y(i) = x_coil[i];
    }
    for (int i = 0; i < 9; ++i) {
        if(i<3){
            out_y(i+NUM_COIL_STATES) = out_u0[i];
        }else if(i<6){
            out_y(i+NUM_COIL_STATES) = out_mL[i-3];
        }else if(i<9){
            out_y(i+NUM_COIL_STATES) = out_nL[i-6];
        }
    }

    return (out_y);
}
