#pragma once
#include <cmath>

template <typename adType>
void CRMSolverIVP(	CRMShootingMethodParams<adType> in_Params,
					adType in_u0[NUM_FLEX_SEG][3], adType in_ftip[3],
                     bool in_FinalValueOnly,
					adType out_x_N[NUM_STATES], adType out_MomentResidual[NUM_RESIDUAL], adType pcoil[NUM_ACT_SET][3], adType Rcoil[NUM_ACT_SET][9],
					double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]){

	adType x_0[NUM_STATES];
	for (int i = 0; i < NUM_STATES; i++) {
		if (i < 3) x_0[i] = in_u0[0][i];
		else if (i < 12) x_0[i] = in_Params.R0[i - 3];
		else if (i < 15) x_0[i] = in_Params.p0[i - 12];
	}

	CRMIVPCoreParams<adType> CoreParams;
	adType u_0[NUM_FLEX_SEG][3], ftip[3];

    adType m_L[NUM_ACT_SET][3] = {0.0, 0.0, 0.0}; //This is the static case, where no input m/n is involved in the trustregion
    adType n_L[NUM_ACT_SET][3] = {0.0, 0.0, 0.0};

    CRMSolverIVP_Prep(x_0, in_Params.IntegrationStepSize,
                      in_Params.Li, in_Params.dlambdainv,
                      in_Params.SegEndLambdas, in_Params.LocMarkerLambdas,
                      in_Params.K, in_Params.Kinv, in_Params.ustar,
                      in_Params.MagMoment, in_Params.fcumlambda,
                      in_Params.B0, in_Params.g, in_Params.actMass, in_Params.actInertia, in_Params.damping, in_Params.DELTA_T,
                      in_Params.v_L_pre, in_Params.w_L_pre, in_Params.p_pre, in_Params.R_pre,
                      m_L, n_L, false, in_FinalValueOnly,CoreParams);

	// copy to local variable
	for (int i = 0; i < 3; i++) ftip[i] = in_ftip[i];

	// We need to pass u_0 as input argument as the values in CoreParams will be overriden with the values provided in the input arguments - functionality needed for solving Boundary Value Problems (BVP)
    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
        for (int j = 0; j < 3; j++) u_0[i][j] = in_u0[i][j];
    }

    CRMSolverIVP_Core ( CoreParams, u_0,ftip, out_x_N, out_MomentResidual, pcoil, Rcoil, out_p_atLocMarkers);

}

template <typename adType>
void CRMSolverIVP_Prep ( adType in_x_0[NUM_STATES], double in_IntegrationStepSize,
						 adType in_Li, double in_dlambdainv,
						 double in_SegEndLambdas[NUM_SEGMENTS], double in_LocMarkerLambdas[NUM_LOCALIZATION_MARKERS],
						 double in_K[NUM_FLEX_SEG][9], double in_Kinv[NUM_FLEX_SEG][9], double in_ustar[NUM_FLEX_SEG][3],
						 adType in_MagMoment[NUM_ACT_SET][3], double in_fcumlambda[NUM_FCUM_LAMBDA+1][3],
                         double in_B0[3], double in_g[3], const double in_actMass[NUM_ACT_SET], double in_actInertia[NUM_ACT_SET][9], double in_damping[NUM_ACT_SET][6], double in_delta_t,
                         double in_v_L_pre[NUM_ACT_SET][3], double in_w_L_pre[NUM_ACT_SET][3], double in_p_pre[NUM_ACT_SET][3],
                         double in_R_pre[NUM_ACT_SET][9], double in_mL[NUM_ACT_SET][3], double in_nL[NUM_ACT_SET][3],
                         bool in_DYNCORE, bool in_FinalValueOnly, CRMIVPCoreParams<adType> &out_CoreParams) {

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
    auto & h0 = out_CoreParams.h0;

	auto & InsertedLength = out_CoreParams.InsertedLength;
	auto & dlambdainv = out_CoreParams.dlambdainv;
	auto & K = out_CoreParams.K;
	auto & Kinv = out_CoreParams.Kinv;
	auto & ustar = out_CoreParams.ustar;
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
											//  States are packed u[0..2],R[0..8],p[0..2]  (R: 3x3 matrix stored in row major order R11 R12 R13 R21 R22 R23 R31 R32 R33)
	mCopy_AB<3>(in_B0,B0);				//  B0 field vector of the MRI scanner (in spatial coordinates)

	dlambdainv=in_dlambdainv;				// reciprocal of dlambda (lambda stepsize used in discretizing fcumlambda)
	FinalValueOnly=in_FinalValueOnly;		// Flag used to indicate if only final value (xf) is returned (true) or if Marker Locations are returned as well (false)
	InsertedLength=in_Li;					// Inserted Length (length of the catheter from the entry point to the tip)
	// If the catheter is inserted more than the length of the catheter, clamp it to catheter length
	if (InsertedLength > SegEndLambdas[NUM_SEGMENTS - 1]) InsertedLength = SegEndLambdas[NUM_SEGMENTS - 1];
//    std::cout << "IVP InsertedLength: " << InsertedLength << std::endl;

	// WE ARE GOING TO REORDER SEGMENT AND ACTUATOR UNITS SO THAT THEY ARE ORDERED FROM THE INSERTION POINT TO THE TIP
	// *  I.E., SWITCH TO PROXIMAL TO DISTAL ORDERING

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
        h0[i] = dVal((SegBounds[2*i+1]-SegBounds[2*i]))/(SegSteps[i]*1.0);
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

	// be careful - order is reversed in K, Kinv, ustar, and MagMoment
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
	}

	mCopy_AB<(NUM_FCUM_LAMBDA+1),3>(in_fcumlambda,fcumlambda);
//	mCopy_AB<3>(in_ftip, ftip);

	/*
	 * The params for dynamics
	 * */
    auto & g = out_CoreParams.g;
    mCopy_AB<3>(in_g, g);               // gravitational vector
    auto & v_L_pre = out_CoreParams.v_L_pre;
    auto & w_L_pre = out_CoreParams.w_L_pre;
//    mCopy_AB<3>(in_v_L_pre, v_L_pre);               // linear velocity vector
//    mCopy_AB<3>(in_w_L_pre, w_L_pre);               // angular velocity vector

    auto & p_pre = out_CoreParams.p_pre;
    auto & R_pre = out_CoreParams.R_pre;
//    mCopy_AB<3>(in_p_pre, p_pre);
//    mCopy_AB<9>(in_R_pre, R_pre);
    auto & damping = out_CoreParams.damping;
//    mCopy_AB<6>(in_damping, damping);
    auto & m_L = out_CoreParams.m_L;
    auto & n_L = out_CoreParams.n_L;

    auto & actMass = out_CoreParams.actMass;
    auto & actInertia = out_CoreParams.actInertia;
    for (int i = 0; i < NUM_ACT_SET; ++i) {
        actMass[i] = in_actMass[i];
        for (int j = 0; j < 9; ++j) {
            actInertia[i][j] = in_actInertia[i][j];
            R_pre[i][j] = in_R_pre[i][j];
        }
        for (int j = 0; j < 3; ++j) {
            v_L_pre[i][j] =   in_v_L_pre[i][j];
            w_L_pre[i][j] =   in_w_L_pre[i][j];
            p_pre[i][j] =   in_p_pre[i][j];
        }
        for (int j = 0; j < 6; ++j) {
            damping[i][j] = in_damping[i][j];
        }
        for (int j = 0; j < 3; ++j) {
            m_L[i][j] = in_mL[i][j];
            n_L[i][j] = in_nL[i][j];
        }

    }

    auto & dyn_signal = out_CoreParams.DYN_Core;
    dyn_signal = in_DYNCORE;

    auto & delta_t = out_CoreParams.DELTA_T;
    delta_t = in_delta_t;

}


template <typename adType>
void CRMSolverIVP_Core ( CRMIVPCoreParams<adType> in_params,
						 adType in_u[NUM_FLEX_SEG][3], adType in_ftip[3],
						 adType out_x_N[NUM_STATES], adType out_Residual[NUM_RESIDUAL], adType out_pcoil[NUM_ACT_SET][3], adType out_Rcoil[NUM_ACT_SET][9],
						 double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]){

//    std::cout << "in_u: " << in_u[0][0] << " " <<  in_u[0][1] << " " <<  in_u[0][2] << std::endl;
//    std::cout << "in_u: " << in_u[1][0] << " " <<  in_u[1][1] << " " <<  in_u[1][2] << std::endl;

    /*
     * Other params for dynamics
     * */
    auto & DYN_SIGNAL = in_params.DYN_Core;
    auto & in_m_L = in_params.m_L;
    auto & in_n_L = in_params.n_L;
    double m_L[NUM_FLEX_SEG][3], n_L[NUM_FLEX_SEG][3];
    for (int i = 0; i < NUM_ACT_SET; ++i) { // NUM_ACT_SET <= NUM_FLEX_SEG
        for (int j = 0; j < 3; ++j) {
            m_L[i][j] = in_m_L[i][j];
            n_L[i][j] = in_n_L[i][j];
        }
    }
    if ( NUM_ACT_SET < NUM_FLEX_SEG){ // last segment is flexible
        for (int i = 0; i < 3; ++i) {
            m_L[NUM_FLEX_SEG-1][i] = 0.0;
            n_L[NUM_FLEX_SEG-1][i] = 0.0;
        }
    }
//    std::cout << "in_u CORE: " << in_u[0] << " " << in_u[1] << " " <<in_u[2] << std::endl;
//    std::cout << "in_n CORE: " << n_L[0] << " " << n_L[1] << " " <<n_L[2] << std::endl;

    adType xi[NUM_STATES];  			//  Initial value of the state for the next segment to be integrated
	adType xf[NUM_STATES];			//  Final value of the state for the last segment integrated
	double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3];	// States at markers (ordered proximal to distal)

	bool loopcondition;				//  local variable for calculating a loop condition

	double l_zero[3]={0.0,0.0,0.0};
	adType RigidSegmentLength;		// Length of the rigid segment - intermediate variable
	adType muhat[9];
	adType RscTB0[3], Tb[3], deltau1[3], K1deltau1[3], K2invResidual[3]; // intermediate variables
	adType tempadType;				// intermediate variables
	adType Residual[NUM_FLEX_SEG][3];				// Residual at the catheter tip -- will be returned

	int   fsegno; 					// flexible segment no
	int   actno, fsegi, fsegip1;	// actuator no, flexible segment before, flexible segment after
	bool  LastSegmentIsRigid=true;	// Flag indicating if the last segment processed is rigid (true) or not (false)

    // we need to copy xi from in_params to the local variable and update it with u[0..2] specified in in_u
	for (int i=0; i<NUM_STATES; i++) xi[i]=in_params.xi[i];
	for (int i=0; i<3; i++) {
		xi[i]=in_u[0][i]; // THe first segment
	}

    // we need to copy in_ftip to local variable
	adType ftip[3];
	mCopy_AB<3,adType>(in_ftip, ftip);
	// we will copy anything we will access more than once (or write to) to local variables
	int   StartSegmentIndex=in_params.StartSegmentIndex;
	int	  NextLocMarker=in_params.NextLocMarker;
	int	  InitialLocMarker=NextLocMarker;
	// for others, we will create aliases
	auto & SegBounds = in_params.SegBounds;
	auto & SegSteps = in_params.SegSteps;
    auto & h0 = in_params.h0;

	auto & InsertedLength = in_params.InsertedLength;
	auto & dlambdainv = in_params.dlambdainv;
	auto & K = in_params.K;
	auto & Kinv = in_params.Kinv;
	auto & ustar = in_params.ustar;
	auto & fcumlambda = in_params.fcumlambda;
	auto & FinalValueOnly = in_params.FinalValueOnly;
	auto & LocMarkers = in_params.LocMarkers;
	auto & B0 = in_params.B0;
	auto & MagMoment = in_params.MagMoment;
	// we will also copy the already filled entries of p_atLocMarkers from params
	for (int i=0; i<InitialLocMarker; i++) {
		for (int j=0; j<3; j++) {
			p_atLocMarkers[i][j]=in_params.p_atLocMarkers[i][j];
		}
	}

    // IMPORTANT NOTE: most proximal segment is assumed to be always flexible
	//    and the flexible and rigid segments are assumed to be alternating
	//    most distal segment can be flexible or rigid

	// Integrate each of the segments
	// If the starting segment is a rigid segment, then we will need to move initial conditions to the start of the next flexible segment and change the StartingSegment to that segment
	if ( StartSegmentIndex % 2 == 1 ) {
		RigidSegmentLength	=	SegBounds[StartSegmentIndex+1]-SegBounds[StartSegmentIndex];	// how far we need to move along the length of the rigid segment to reach the next flexible segment
		for (int i=0; i<3; i++) { 										// u[0..2] and R[0..8] remain the same
			xi[i+9+3] = xi[i+9+3] + RigidSegmentLength * xi[3+i*3+2];	// p[0..2] will translate along the z direction of the R matrix (3rd column)
			Residual[0][i]=0.0;											// Residual is initialized to the 0 vector, just in case the rigid segment is the only segment (it is both the start and the end segment),
		}																//     In this case, the main integration loop will not execute, so, we need to have valid return values
        mCopy_AB<NUM_STATES>(xi,xf);									//     for both xf and Residual
		StartSegmentIndex++;  // move the start segment to the next segment
		if (!FinalValueOnly) {
			// are there any localization markers?  If so, calculate their positions
			loopcondition = ((NextLocMarker<NUM_LOCALIZATION_MARKERS) && (LocMarkers[NextLocMarker]<=SegBounds[StartSegmentIndex]));
			while ( loopcondition ) {  // remember StartSegmentIndex has already been incremented
				tempadType =LocMarkers[NextLocMarker]-SegBounds[StartSegmentIndex];  	// this would be a negative number
				LocMarkerUpdate(p_atLocMarkers[NextLocMarker], xi, tempadType);			//	xi has already been updated
				NextLocMarker++;
				loopcondition = ((NextLocMarker<NUM_LOCALIZATION_MARKERS) && (LocMarkers[NextLocMarker]<=SegBounds[StartSegmentIndex]));	}
		}
	}


    for (int i=StartSegmentIndex; i< NUM_SEGMENTS; i++){
		if ( i%2 == 0 ) {  // Flexible Segment
			LastSegmentIsRigid=false;
			// Prepare the CRMIntegrand Parameters
			fsegno=i>>1; // i/2, flexible segment no

            // Integrate
            ABM4( 	xi, SegBounds[i], SegSteps[fsegno], h0[fsegno],
                     InsertedLength, dlambdainv, K[fsegno], Kinv[fsegno], l_zero, ustar[fsegno], n_L[fsegno],
                     fcumlambda, ftip,
                     FinalValueOnly, LocMarkers, &NextLocMarker,
                     xf, p_atLocMarkers);

            if (i == NUM_SEGMENTS-1){ // last flexible segment
                mSub_AB<3,1>( &(xf[0]) , ustar[fsegno], deltau1);
                mMult_AB<3,3,1>( K[fsegno], deltau1, K1deltau1 );

                for (int j = 0; j < 3; ++j) {
                    Residual[fsegno][j] = K1deltau1[j];
                }
            }
        }
		else {  // Need to do actuation/rigid segment calculations to transfer Initial Conditions to next flexible segment
			LastSegmentIsRigid=true;
			RigidSegmentLength=(SegBounds[i+1]-SegBounds[i]);
			// R
			for (int j=0; j<9; j++) {
				xi[3+j]=xf[3+j];
			}
			// p
			for (int j=0; j<3; j++) {
				xi[3+9+j]=xf[3+9+j]+xf[3+j*3+2]*RigidSegmentLength;
			}

            // u
            actno=(i-1)>>1;		// actuator no
            fsegi=actno;
            fsegip1=actno+1;
            // Tb=\mu_c \cross R_sc^T B0,s
            wHat(MagMoment[actno],muhat);
            mMult_ATB<3,3,1>(&(xf[3]),B0,RscTB0);
            mMult_AB<3,3,1>(muhat,RscTB0,Tb);
//                std::cout << "Tb in core: " << Tb[0] << " " << Tb[1] << " " << Tb[2] << std::endl;

			// u2=u2star + ( K2inv K1 (u1 - u1star ) - K2inv Tb )
            mSub_AB<3,1>( &(xf[0]) , ustar[fsegi], deltau1);
            mMult_AB<3,3,1>( K[fsegi], deltau1, K1deltau1 );

            double res[3];
            if (!DYN_SIGNAL){
                mSub_AB<3,1>( K1deltau1 , Tb, res); // Quasi-static residual
            }else{
                mSub_AB<3,1>( m_L[fsegi], K1deltau1 , res);
            }
            for (int j = 0; j < 3; ++j) {
                Residual[fsegi][j] = res[j];
            }

            for (int j=0; j<3; j++) {
                out_pcoil[actno][j]=xf[3+9+j]+xf[3+j*3+2]*RigidSegmentLength * 0.5;
            }

            for (int j = 0; j < 9; ++j) {
                out_Rcoil[actno][j] = xf[j+3];
            }

            if (i < (NUM_SEGMENTS-1)) {	// we want to make sure that we are not at the last segment
//                mCopy_AB<3>( ustar[fsegip1], &(xi[0]) );
//                mMult_AB<3,3,1>( Kinv[fsegip1], Residual, K2invResidual );
//                mAdd_AB<3,1>( ustar[fsegip1], K2invResidual, &(xi[0]) );

                for (int j = 0; j < 3; ++j) {
                    xi[j] = in_u[fsegip1][j];
                }

            }else {	// otherwise, we are at the last segment, and we need to copy the R and p values calculated for xi to xf so that they can be returned
				for (int j=0; j<NUM_STATES; j++) {
					if (j<3) xf[j]=0.0;  // u[0..2] are assigned to zero
					else xf[j]=xi[j];
				}
			}

			if (!FinalValueOnly) {
				// are there any localization markers?  If so, calculate their positions
				loopcondition = ((NextLocMarker<NUM_LOCALIZATION_MARKERS) && (LocMarkers[NextLocMarker]<=SegBounds[i+1]));
				while ( loopcondition ) {
					tempadType =LocMarkers[NextLocMarker]-SegBounds[i+1];				// this would be a negative value
					LocMarkerUpdate(p_atLocMarkers[NextLocMarker], xi, tempadType);		//  xi has already been updated, it is the end point position
					NextLocMarker++;
					loopcondition = ((NextLocMarker<NUM_LOCALIZATION_MARKERS) && (LocMarkers[NextLocMarker]<=SegBounds[i+1]));
				}
			}

		}

	}


	// Copy marker locations to the output
	if (!FinalValueOnly) {
		for (int i=0; i<NUM_LOCALIZATION_MARKERS; i++) {
			for (int j=0; j<3; j++) {
				out_p_atLocMarkers[i][j]=p_atLocMarkers[(NUM_LOCALIZATION_MARKERS-1)-i][j];
			}
		}
	}

	// Copy the final values of the state x_f to the output
	mCopy_AB<NUM_STATES>(xf, out_x_N);

    // Copy the residual to the output
//	mCopy_AB<NUM_RESIDUAL>(Residual, out_Residual);

    for (int i = 0; i < NUM_FLEX_SEG; ++i) {
        for (int j = 0; j < 3; ++j) {
            out_Residual[j+3*i] = Residual[i][j];
        }
    }

//    std::cout << "Residual in core: " << Residual[0] << " " << Residual[1] << " " << Residual[2] << std::endl;

}

template <typename adType>
void LocMarkerUpdate(double p[3], adType xi[NUM_STATES], adType t) {
	p[0] = xi[12] + t * xi[5];
	p[1] = xi[13] + t * xi[8];
	p[2] = xi[14] + t * xi[11];
}


template <typename adType>
void CRMIntegrand (	adType s, adType x[NUM_STATES], adType Li, double dlambdainv,
					const double in_K[9], const double in_Kinv[9], const double in_l[3], const double in_ustar[3], adType in_nL[3],
                    double in_fcumlambda[NUM_FCUM_LAMBDA+1][3], adType in_ftip[3],
					adType xdot[NUM_INTEGRATION_STATES]) {

	adType Length;
	double deltalambdainv;
#ifdef ANALYTICAL_SE3_STEP
	adType u[3], R[9]; 	// p[3];   				we will not need this for analytical calculation
	adType udot[3]; 	// Rdot[9], pdot[3];	we will not need this for analytical calculation
#else
	adType u[3], R[9], p[3];
	adType udot[3], Rdot[9], pdot[3];
#endif
	adType K[9], Kinv[9], l[3], fcum[3], ustar[3], ustardot[3];  //  We are assuming Kdot=0.0 (K=const)

	// copy inputs and parameters to local variables
	for (int i=0; i<NUM_STATES; i++) {
		if (i<3) { 			// 0 <= i < 3
			u[i]=x[i];
		}
		else if (i<12) { 	// 3 <= i < 12
			R[i-3]=x[i];
		}
#ifndef ANALYTICAL_SE3_STEP
		// we will not need this for analytical calculation
		else { 				// 12 <= i < 15
			p[i - 12] = x[i];
		}
#endif
	}
	Length=Li;
	deltalambdainv=dlambdainv;
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			K[i*3+j]		=	in_K[i*3+j];
			Kinv[i*3+j]		=	in_Kinv[i*3+j];
		}
		l[i]			=	in_l[i];
		ustar[i]		=	in_ustar[i];
		ustardot[i]		=	0.0; //in_ustardot[i]; // we assume ustardot=0.0 since our rest shape model is piecewise constant curvature
	}

	// calculate interpolated value of fcum
	adType lambda=Length-s;
	adType ix=lambda*deltalambdainv;
	double ird_f=floor(dVal(ix)); // index for round down  -- doubleing point
	if (ird_f<0) ird_f=0;
	int ird=(int) ird_f;	//    integer index
	double iru_f=ceil(dVal(ix)); 	// index for round up  -- doubleing point
	if (iru_f>NUM_FCUM_LAMBDA) iru_f=NUM_FCUM_LAMBDA;
	int iru=(int) iru_f;	//    integer index
	adType ixmird=ix-ird;    // weight for interpolation
	adType irumix=iru-ix;	// weight for interpolation
	for (int i=0; i<3; i++) {
		fcum[i]		=	(in_fcumlambda[iru][i] * ixmird + in_fcumlambda[ird][i] * irumix);
	}

	// add the tip force to fcum
	for (int i = 0; i < 3; i++) {
		fcum[i]		+=	in_ftip[i];
	}
//    fcum[0] = fcum[1] = fcum[2] = 0.0;


    adType nL_spatial[3];
    mMult_AB<3,3,1>(R, in_nL, nL_spatial);
    // add the tip force to fcum
    for (int i = 0; i < 3; i++) {
        fcum[i]		+=	nL_spatial[i];
    }

//    std::cout << "fcum: " << fcum[0] << " " << fcum[1] << " " << fcum[2] << std::endl;

	// calculate u_hat
	adType u_hat[9];
	wHat(u,u_hat);

	// udot = ustardot - Kinv*((um*K+Kdot)*(u-ustar_s) + e3m*R'*intf + R'*l); % udot
	//
	//   e3hat*R' = [ -r12 -r22 -r32; r11 r21 r31; 0 0 0];
	adType e3hatRT[9];
	e3hatRT[0]=-R[1]; 	e3hatRT[1]=-R[4]; 	e3hatRT[2]=-R[7];
	e3hatRT[3]=R[0];	e3hatRT[4]=R[3];	e3hatRT[5]=R[6];
	e3hatRT[6]=0.0;	e3hatRT[7]=0.0;	e3hatRT[8]=0.0;
	adType e3hatRTfcum[3];
	mMult_AB<3,3,1>(e3hatRT, fcum, e3hatRTfcum);			//  e3m*R'*intf
//    std::cout << "e3hatRTfcum: " << e3hatRTfcum[0] << " " << e3hatRTfcum[1] << " " << e3hatRTfcum[2] << std::endl;

    adType RTl[3];
	mMult_ATB<3,3,1>(R,l,RTl);								// R'*l
	adType umustar[3];
	mSub_AB<3,1>(u,ustar,umustar);							// (u-ustar_s)
	adType Kumustar[3], uhatKumustar[3];
	mMult_AB<3,3,1>(K,umustar,Kumustar);
	mMult_AB<3,3,1>(u_hat,Kumustar,uhatKumustar); 		 	//(um*K+Kdot)*(u-ustar_s)  assuming Kdot=0
	adType sumterm[3];
	mAdd_ABC<3,1>(uhatKumustar,e3hatRTfcum,RTl,sumterm);	// ((um*K+Kdot)*(u-ustar_s) + e3m*R'*intf + R'*l)
	adType KinvSum[3];
	mMult_AB<3,3,1>(Kinv,sumterm,KinvSum);					// Kinv*((um*K+Kdot)*(u-ustar_s) + e3m*R'*intf + R'*l)
	mSub_AB<3,1>(ustardot,KinvSum,udot);					// udot = ustardot - Kinv*((um*K+Kdot)*(u-ustar_s) + e3m*R'*intf + R'*l);

#ifndef ANALYTICAL_SE3_STEP
	// we will not need these for analytical calculation
	// Rdot = R*u_hat
	mMult_AB<3, 3, 3>(R, u_hat, Rdot);
	// pdot = R*e3,
	for (int i = 0; i < 3; i++) {
		pdot[i] = R[i * 3 + 2];
	}
#endif

	//xdot(1:3) = ustardot - Kinv*((um*K+Kdot)*(u-ustar_s) + e3m*R'*intf + R'*l); % udot
	//xdot(4:12) = reshape(R*um,9,1); % Rdot   --- Note that matlab code reshapes in column major order while we are saving in row major order
	//xdot(13:15) = R*e3;             % pdot
	//  note: the sample code has matlab indexing starting from 1 to 15

    for (int i=0; i< 3; i++) {
#ifdef ANALYTICAL_SE3_STEP
		xdot[i] = udot[i];
#else
		if (i < 3) { // 0 <= i < 3
			xdot[i] = udot[i];
		}
		// we will not need these for analytical calculation
		else if (i < 12) { // 3 <= i < 12
			xdot[i] = Rdot[i - 3];
		}
		else if (i < 15) { // 12 <= i < 15
			xdot[i] = pdot[i - 12];
		}
#endif
	}

}



//
//
//	NUMERICAL INTEGRATION FUNCTIONS
//
//

//function [x_1toN] = ABM4(x_0, t_0, N, h, Integrand, initmethod)
template <typename adType>
void ABM4 (	adType in_x_0[NUM_STATES], adType t_0, int N, double h,
			adType Li, double dlambdainv, double in_K[9], double in_Kinv[9], double in_l[3], double in_ustar[3], adType in_nL[3],
            double in_fcumlambda[NUM_FCUM_LAMBDA+1][3], adType in_ftip[3], bool FinalValueOnly,
            double	in_LocMarkers[NUM_LOCALIZATION_MARKERS], int *inout_NextLocMarkerIdx,
               adType out_x_N[NUM_STATES], double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]  ) {

    adType x_nm3[NUM_STATES];
	adType x_nm2[NUM_STATES];
	adType x_nm1[NUM_STATES];
	adType x_n[NUM_STATES];
	adType x_np1[NUM_STATES];
	adType t_n;
	adType xdot_nm3[NUM_INTEGRATION_STATES];
	adType xdot_nm2[NUM_INTEGRATION_STATES];
	adType xdot_nm1[NUM_INTEGRATION_STATES];
	adType xdot_n[NUM_INTEGRATION_STATES];
	bool loopcondition;

	int NextLocMarkerIdx = *inout_NextLocMarkerIdx;
	double LocMarkers[NUM_LOCALIZATION_MARKERS];
	double delta_n_overh,delta_np1_overh;

//    adType fcum[3];

	// initialize the iteration items
	t_n = t_0;
	for (int i = 0; i < NUM_STATES; i++) {
		x_n[i] = in_x_0[i];
	}
	for (int i=0; i<NUM_LOCALIZATION_MARKERS; i++) {
		LocMarkers[i]=in_LocMarkers[i];
	}
//	for (int i = 0; i < NUM_INTEGRATION_STATES; i++) {
//		xdot_nm1[i]=0.0;
//		xdot_nm2[i]=0.0;
//		xdot_nm3[i]=0.0;
//	}


	for (int idx=0; idx<N; idx++) {

        if (idx<3) {  // RK2 initialization steps
			RK2_step(x_n, t_n, h, Li, dlambdainv, in_K, in_Kinv, in_l, in_ustar,  in_nL, in_fcumlambda, in_ftip, x_np1, xdot_n);
		}
		else { 		 // ABM4 steps
			ABM4_step(x_n, t_n, h, xdot_nm1, xdot_nm2, xdot_nm3, x_nm1, x_nm2, x_nm3,
                      Li, dlambdainv, in_K, in_Kinv, in_l, in_ustar, in_nL, in_fcumlambda, in_ftip, x_np1, xdot_n);
            if ( isnan(x_n[0]) ) {
                std::cout << "FLY ME TO THE MOON!! " << std::endl;
//                exit( 3 );
            }
		}

		// increment "time"
		t_n=t_n+h;

		if (!FinalValueOnly) {
			// are there any localization markers?  If so, calculate their positions
			// remember, t_n has already been incremented
			loopcondition = (NextLocMarkerIdx<NUM_LOCALIZATION_MARKERS ) && (LocMarkers[NextLocMarkerIdx]<=t_n);
			while ( loopcondition ) {
				delta_np1_overh=(dVal(t_n)-LocMarkers[NextLocMarkerIdx])/h;
				delta_n_overh=1.0-delta_np1_overh;
				out_p_atLocMarkers[NextLocMarkerIdx][0]=delta_np1_overh*dVal(x_n[12])+delta_n_overh*dVal(x_np1[12]);
				out_p_atLocMarkers[NextLocMarkerIdx][1]=delta_np1_overh*dVal(x_n[13])+delta_n_overh*dVal(x_np1[13]);
				out_p_atLocMarkers[NextLocMarkerIdx][2]=delta_np1_overh*dVal(x_n[14])+delta_n_overh*dVal(x_np1[14]);
				NextLocMarkerIdx++;
				loopcondition = (NextLocMarkerIdx<NUM_LOCALIZATION_MARKERS ) && (LocMarkers[NextLocMarkerIdx]<=t_n);
			}
		}

		// update the iteration items
		for (int i = 0; i < NUM_STATES; i++) {
			x_nm3[i] = x_nm2[i];
			x_nm2[i] = x_nm1[i];
			x_nm1[i] = x_n[i];
			x_n[i] = x_np1[i];
		}
		for (int i = 0; i < NUM_INTEGRATION_STATES; i++) {
			xdot_nm3[i]=xdot_nm2[i];
			xdot_nm2[i]=xdot_nm1[i];
			xdot_nm1[i]=xdot_n[i];
		}
	}

    // Copy final value
	for (int i=0; i<NUM_STATES; i++) {
		out_x_N[i]=x_n[i];
	}
	*inout_NextLocMarkerIdx=NextLocMarkerIdx;

}


// [x_np1, xdot_n, xdot_nm1, xdot_nm2] = ABM4_step(x_n, t_n, xdot_nm1, xdot_nm2, xdot_nm3, h, Integrand)
template <typename adType>
void ABM4_step(	adType in_x_n[NUM_STATES], adType t_n, double h,
				adType in_xdot_nm1[NUM_INTEGRATION_STATES], adType in_xdot_nm2[NUM_INTEGRATION_STATES], adType in_xdot_nm3[NUM_INTEGRATION_STATES],
				adType in_x_nm1[NUM_STATES], adType in_x_nm2[NUM_STATES], adType in_x_nm3[NUM_STATES],
				adType Li, double dlambdainv, double in_K[9], double in_Kinv[9], double in_l[3], double in_ustar[3], double in_nL[3], double in_fcumlambda[NUM_FCUM_LAMBDA+1][3], adType in_ftip[3],
                adType out_x_np1[NUM_STATES], adType out_xdot_n[NUM_INTEGRATION_STATES]) {

	const double P_COEFF_N=55.0/24.0, P_COEFF_Nm1=-59.0/24.0, P_COEFF_Nm2=37.0/24.0, P_COEFF_Nm3=-9.0/24.0;  // AB4 Predictor Coefficients
	const double C_COEFF_Np1=9.0/24.0, C_COEFF_N=19.0/24.0, C_COEFF_Nm1=-5.0/24.0, C_COEFF_Nm2=1.0/24.0;     // AM4 Corrector Coefficients
	adType x_n[NUM_STATES];       				// from input
	adType x_nm1[NUM_STATES];       				// from input
	adType x_nm2[NUM_STATES];       				// from input
	adType x_nm3[NUM_STATES];       				// from input
	adType u_n[3], R_n[9], p_n[3];				// variables used in analytical calculation
	adType x_np1_hat[NUM_STATES];    			// intermediate
	adType xdot_np1_hat[NUM_INTEGRATION_STATES];	// intermediate
	adType xdot_n[NUM_INTEGRATION_STATES];  		// for output
	adType xdot_nm1[NUM_INTEGRATION_STATES];  	// from input
	adType xdot_nm2[NUM_INTEGRATION_STATES];  	// from input
	adType xdot_nm3[NUM_INTEGRATION_STATES];  	// from input

    for (int i = 0; i < NUM_STATES; i++) {
        if (i<3)
            u_n[i]    =  in_x_n[i];
        else if (i<(9+3))
            R_n[i-3]  =  in_x_n[i];
        else if (i<(3+9+3))
            p_n[i-12] =  in_x_n[i];
        x_n[i] = in_x_n[i];
        x_nm1[i] = in_x_nm1[i];
        x_nm2[i] = in_x_nm2[i];
        x_nm3[i] = in_x_nm3[i];
    }
    for (int i = 0; i < NUM_INTEGRATION_STATES; i++) {
        xdot_nm1[i]=in_xdot_nm1[i];
        xdot_nm2[i]=in_xdot_nm2[i];
        xdot_nm3[i]=in_xdot_nm3[i];
    }


    //ABM4_STEP_STEP1:
	CRMIntegrand(t_n, x_n, Li, dlambdainv, in_K, in_Kinv, in_l, in_ustar, in_nL, in_fcumlambda, in_ftip, xdot_n);

    for (int i=0; i<3; i++) {
        x_np1_hat[i]    = x_n[i] + h * ( P_COEFF_N * xdot_n[i] + P_COEFF_Nm1 * xdot_nm1[i] + P_COEFF_Nm2 * xdot_nm2[i] + P_COEFF_Nm3 * xdot_nm3[i] );
    }
    for (int i=3; i<NUM_INTEGRATION_STATES; i++) {
        x_np1_hat[i+9+3]    = x_n[i+9+3] + h * ( P_COEFF_N * xdot_n[i] + P_COEFF_Nm1 * xdot_nm1[i] + P_COEFF_Nm2 * xdot_nm2[i] + P_COEFF_Nm3 * xdot_nm3[i] );
    }
#ifdef ANALYTICAL_SE3_STEP
	//      calculate R_np1_hat and p_np1_hat analytically, without numerical integration
	adType u_n_pred[3];
	for (int i = 0; i < 3; i++) u_n_pred[i] = (P_COEFF_N * x_n[i] + P_COEFF_Nm1 * x_nm1[i] + P_COEFF_Nm2 * x_nm2[i] + P_COEFF_Nm3 * x_nm3[i]);
	SE3_Analytical_Step(R_n, p_n, u_n_pred, h, x_np1_hat + 3 /*R_np1_hat*/, x_np1_hat + 12 /*p_np1_hat*/);
#endif
	//ABM4_STEP_STEP2:
	CRMIntegrand(t_n+h, x_np1_hat, Li, dlambdainv, in_K, in_Kinv, in_l, in_ustar, in_nL, in_fcumlambda, in_ftip, xdot_np1_hat);

    for (int i=0; i<3; i++) {
        out_x_np1[i]    = x_n[i] + h * ( C_COEFF_Np1 * xdot_np1_hat[i] + C_COEFF_N * xdot_n[i] + C_COEFF_Nm1 * xdot_nm1[i] + C_COEFF_Nm2 * xdot_nm2[i] );
    }
    for (int i=3; i<NUM_INTEGRATION_STATES; i++) {
        out_x_np1[i+9+3]    = x_n[i+9+3] + h * ( C_COEFF_Np1 * xdot_np1_hat[i] + C_COEFF_N * xdot_n[i] + C_COEFF_Nm1 * xdot_nm1[i] + C_COEFF_Nm2 * xdot_nm2[i] );
    }
#ifdef ANALYTICAL_SE3_STEP
	//      calculate R_np1 and p_np1 analytically, without numerical integration
	adType u_n_corr[3];
	for (int i = 0; i < 3; i++) u_n_corr[i] = (C_COEFF_Np1 * x_np1_hat[i] + C_COEFF_N * x_n[i] + C_COEFF_Nm1 * x_nm1[i] + C_COEFF_Nm2 * x_nm2[i]);
	SE3_Analytical_Step(R_n, p_n, u_n_corr, h, out_x_np1 + 3 /*R_np1*/, out_x_np1 + 12 /*p_np1*/);
#endif

	// copy to output
	for (int i=0; i<NUM_INTEGRATION_STATES; i++) {
		out_xdot_n[i] 	= xdot_n[i];
	}

}


//[x_np1, xdot_n] = RK2_step(x_n, t_n, h, Integrand)
template <typename adType>
void RK2_step(	adType in_x_n[NUM_STATES], adType t_n, double h,
				adType Li, double dlambdainv, double in_K[9], double in_Kinv[9], double in_l[3], double in_ustar[3], double in_nL[3], double in_fcumlambda[NUM_FCUM_LAMBDA+1][3], adType in_ftip[3],
                adType out_x_np1[NUM_STATES], adType out_xdot_n[NUM_INTEGRATION_STATES] ) {

    adType x_n[NUM_STATES];       // from input
	adType k1[NUM_INTEGRATION_STATES];
	adType k2oh[NUM_INTEGRATION_STATES];
	adType x_n_p_k1o2[NUM_STATES];
	adType xdot_n[NUM_INTEGRATION_STATES];

	adType u_n[3], R_n[9], p_n[3];				// variables used in analytical calculation
	for (int i = 0; i < NUM_STATES; i++) {
		if (i < 3)
			u_n[i] = in_x_n[i];
		else if (i < (9 + 3))
			R_n[i - 3] = in_x_n[i];
		else  if (i < (3+ 9 + 3))
			p_n[i - 12] = in_x_n[i];
        x_n[i] = in_x_n[i];
	}

	//RK2_STEP_STEP1:
	CRMIntegrand(t_n, x_n, Li, dlambdainv, in_K, in_Kinv, in_l, in_ustar, in_nL, in_fcumlambda, in_ftip, xdot_n);
	for (int i=0; i<3; i++) {
		k1[i] 			= h * xdot_n[i];
		x_n_p_k1o2[i] 	= x_n[i] + k1[i] * 0.5;
	}
    for (int i = 3; i < NUM_INTEGRATION_STATES; i++) { //twists
        k1[i] 			= h * xdot_n[i];
        x_n_p_k1o2[i+9+3] 	= x_n[i+9+3] + k1[i] * 0.5;
    }
#ifdef ANALYTICAL_SE3_STEP
	// we will calculate R_np1half and p_np1half analytically, without numerical integration
	adType R_np1half[9], p_np1half[3];
	SE3_Analytical_Step(R_n, p_n, u_n, h * 0.5, R_np1half, p_np1half);
	// copy these to x_n_p_k1o2
	for (int i = 3; i < 3+9+3; i++){ // In 3-15 indices
		if (i < (9 + 3))
			x_n_p_k1o2[i] = R_np1half[i - 3];	//R
		else
			x_n_p_k1o2[i] = p_np1half[i - 12];	//p
	}
#endif

	//RK2_STEP_STEP2:
	CRMIntegrand(t_n + h * 0.5, x_n_p_k1o2, Li, dlambdainv, in_K, in_Kinv, in_l, in_ustar, in_nL, in_fcumlambda, in_ftip, k2oh);

    for (int i = 0; i < 3; i++) {
        out_x_np1[i] = x_n[i] + h * k2oh[i];
    }
    for (int i = 3; i < NUM_INTEGRATION_STATES; i++) { //twists
        out_x_np1[i+9+3] = x_n[i+9+3] + h * k2oh[i];
    }

#ifdef ANALYTICAL_SE3_STEP
	// we will calculate R_np1 and p_np1 analytically, without numerical integration
	adType R_np1[9], p_np1[3];
	SE3_Analytical_Step(R_n, p_n, x_n_p_k1o2/*u_np1half*/, h, R_np1, p_np1);
	// copy these to the output state
	for (int i = 3; i < 3+9+3; i++) {
		if (i < (9 + 3))
			out_x_np1[i] = R_np1[i - 3];	//R
		else
			out_x_np1[i] = p_np1[i - 12];	//p
	}

#endif
	for (int i = 0; i < NUM_INTEGRATION_STATES; i++) {
		out_xdot_n[i] = xdot_n[i];
	}

}


// definitions needed for twist exponential calculation
#define EPS 1.0e-12   // the threshold for assuming ||u|| to be approximately 0, so that we should use pure translation equation


// calculate R_np1 and p_np1 analytically using twist exponential, without numerical integration
//   g_np1 = g_n * expm ( \hat{\xi}^b *h ),  where \xi^b= [ 0 0 1 u_n^T ]^T
//   g = [R p; 0 0 0 1];
template <typename adType>
void SE3_Analytical_Step(adType in_R_n[9], adType in_p_n[3], adType in_u_n[3], double h, adType out_R_np1[9], adType out_p_np1[3]) {
	adType R_n[9], p_n[3], u_n[3];

	mCopy_AB<3*3>(in_R_n, R_n);
	mCopy_AB<3>(in_p_n, p_n);
	mCopy_AB<3>(in_u_n, u_n);

	// calculate R_np1 and p_np1 analytically, without numerical integration
	adType Rdelta[9], pdelta[3];
	adType umagsq, umagsqresp, umag, umagresp, unorm[3], delsumag, uu3dels[3], ImRuxv[3], ImRuxvpuuTvds[3], Rnpd[3];
	umagsq = vNormSq<3>(u_n);
	if (umagsq < EPS) {
		mCopy_AB<3*3>(R_n, out_R_np1);
		out_p_np1[0] = p_n[0] + R_n[2] * h;
		out_p_np1[1] = p_n[1] + R_n[5] * h;
		out_p_np1[2] = p_n[2] + R_n[8] * h;
	}
	else {
		umagsqresp = 1.0 / umagsq;
		umag = sqrt(umagsq);
		umagresp = 1.0 / umag;
		mMult_sA<3, 1>(umagresp, u_n, unorm);
		delsumag = h * umag;
		RodriguesExpanded(unorm, delsumag, Rdelta);
		mMult_sA<3, 1>(u_n[2] * h, u_n, uu3dels);
		ImRuxv[0] = Rdelta[0 * 3 + 1] * u_n[0] - Rdelta[0 * 3 + 0] * u_n[1] + u_n[1];
		ImRuxv[1] = Rdelta[1 * 3 + 1] * u_n[0] - Rdelta[1 * 3 + 0] * u_n[1] - u_n[0];
		ImRuxv[2] = Rdelta[2 * 3 + 1] * u_n[0] - Rdelta[2 * 3 + 0] * u_n[1];
		mAdd_AB<3, 1>(ImRuxv, uu3dels, ImRuxvpuuTvds);
		mMult_sA<3, 1>(umagsqresp, ImRuxvpuuTvds, pdelta);
		mMult_AB<3, 3, 3>(R_n, Rdelta, out_R_np1);
		mMult_AB<3, 3, 1>(R_n, pdelta, Rnpd);
		mAdd_AB<3, 1>(p_n, Rnpd, out_p_np1);
	}
}

#undef EPS

//
//
//  Robotic Kinematics Related Functions
//
//

template <typename adType>
void wHat(adType in_w[3], adType out_what[9]) {  // what is stored as a 1-dim array in row major order

	out_what[0] = 0.0;
	out_what[1] = -in_w[2];
	out_what[2] = in_w[1];
	out_what[3] = in_w[2];
	out_what[4] = 0.0;
	out_what[5] = -in_w[0];
	out_what[6] = -in_w[1];
	out_what[7] = in_w[0];
	out_what[8] = 0.0;

}

template <typename adType>
void RodriguesFormula(adType in_w[3], adType in_theta, adType out_R[9]) {

	adType I3x3[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };

	adType w[3], what[9], whatsq[9], temp1[9], temp2[9], temp3[9];
	adType c = cos(in_theta);
	adType s = sin(in_theta);
	adType v = 1.0 - c;

	for (int i = 0; i < 3; i++) {
		w[i] = in_w[i];
	}

	wHat(w, what);
	mMult_sA<3, 3>(s, what, temp1);
	mMult_AB<3, 3, 3>(what, what, whatsq);
	mMult_sA<3, 3>(v, whatsq, temp2);
	mAdd_AB<3, 3>(I3x3, temp1, temp3);
	mAdd_AB<3, 3>(temp2, temp3, out_R);

}


template <typename adType>
void RodriguesExpanded(adType in_w[3], adType in_theta, adType out_R[9]) {

	adType w1 = in_w[0];
	adType w2 = in_w[1];
	adType w3 = in_w[2];
	adType c = cos(in_theta);
	adType s = sin(in_theta);
	adType v = 1.0 - c;
	adType w1s = w1 * s;
	adType w2s = w2 * s;
	adType w3s = w3 * s;
	adType w1w2v = w1 * w2 * v;
	adType w1w3v = w1 * w3 * v;
	adType w2w3v = w2 * w3 * v;
	out_R[0 * 3 + 0] = w1 * w1 * v + c;
	out_R[0 * 3 + 1] = w1w2v - w3s;
	out_R[0 * 3 + 2] = w1w3v + w2s;
	out_R[1 * 3 + 0] = w1w2v + w3s;
	out_R[1 * 3 + 1] = w2 * w2 * v + c;
	out_R[1 * 3 + 2] = w2w3v - w1s;
	out_R[2 * 3 + 0] = w1w3v - w2s;
	out_R[2 * 3 + 1] = w2w3v + w1s;
	out_R[2 * 3 + 2] = w3 * w3 * v + c;

}
