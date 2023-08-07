#pragma once
#include <cmath>


//
//
//	NUMERICAL INTEGRATION FUNCTIONS
//
//


template <typename adType>
void CRMFlexible_IVP_Back ( int SegmentIndex, const double in_p[3], const double in_R[9],  CRMIVPCoreParams<adType> in_params,
                       const double in_u[3], const double in_n_L[3],
                       adType out_u[3], double out_p[3], double out_R[9]){

//    std::cout << "in_u: " << in_u[0][0] << " " <<  in_u[0][1] << " " <<  in_u[0][2] << std::endl;
//    std::cout << "in_u: " << in_u[1][0] << " " <<  in_u[1][1] << " " <<  in_u[1][2] << std::endl;

    double n_L[3];
    for (int i = 0; i < 3; ++i) {
        n_L[i] = in_n_L[i];
    }
    adType xi[NUM_STATES];  			//  Initial value of the state for the next segment to be integrated
    adType xf[NUM_STATES];			//  Final value of the state for the last segment integrated
    adType h;

    double l_zero[3]={0.0,0.0,0.0};
    adType RigidSegmentLength;		// Length of the rigid segment - intermediate variable
//    adType tempadType;				// intermediate variables

    int   fsegno; 					// flexible segment no

    for (int i = 0; i < NUM_STATES; i++) {
        if (i < 3) {
            xi[i] = in_u[i];
        }
        else if (i < 12) {
            xi[i] = in_R[i - 3];
        }
        else if (i < 15) {
            xi[i] = in_p[i - 12];
        }
    }

    // we need to copy in_ftip to local variable
    adType ftip[3] = {0.0,0.0,0.0}; //placeholder
    // we will copy anything we will access more than once (or write to) to local variables
    int	  NextLocMarker=in_params.NextLocMarker;
    int	  InitialLocMarker=NextLocMarker;
    // for others, we will create aliases
    auto & SegBounds = in_params.SegBounds;
    auto & SegSteps = in_params.SegSteps;

    auto & InsertedLength = in_params.InsertedLength;
    auto & dlambdainv = in_params.dlambdainv;
    auto & K = in_params.K;
    auto & Kinv = in_params.Kinv;
    auto & ustar = in_params.ustar;
    auto & fcumlambda = in_params.fcumlambda;
    auto & FinalValueOnly = in_params.FinalValueOnly;
    auto & LocMarkers = in_params.LocMarkers;
    double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3];	// Positions at markers (ordered proximal to distal)
    // we will also copy the already filled entries of p_atLocMarkers from params
    for (int i = 0; i < InitialLocMarker; i++) {
        for (int j = 0; j < 3; j++) {
            p_atLocMarkers[i][j] = in_params.p_atLocMarkers[i][j];
        }
    }
    // Flexible Segment
    // Prepare the CRMIntegrand Parameters
    fsegno=SegmentIndex>>1; // i/2, flexible segment no
    h= -1* (SegBounds[SegmentIndex+1] - SegBounds[SegmentIndex] )/(SegSteps[fsegno]*1.0);

    // Integrate backwards
    ABM4( 	xi, SegBounds[SegmentIndex+1], SegSteps[fsegno], h,
             InsertedLength, dlambdainv, K[fsegno], Kinv[fsegno], l_zero, ustar[fsegno], n_L,
             fcumlambda, ftip,
             FinalValueOnly, LocMarkers, &NextLocMarker,
             xf, p_atLocMarkers);

    for (int j=0; j < 3; j++) {
        out_p[j]=xf[3+9+j];
    }
    for (int j = 0; j < 9; ++j) {
        out_R[j] = xf[j+3];
    }
    mCopy_AB<3>(&(xf[0]), out_u);

}




template <typename adType>
void CRMFlexForward_pass (  int SegmentIndex, const double in_p[3], const double in_R[9],  CRMIVPCoreParams<adType> in_params,
                            const double in_u[3], const double in_n_L[3],
                            adType out_u[3], double out_p[3], double out_R[9], double out_p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3]){

    double n_L[3];
    for (int i = 0; i < 3; ++i) {
        n_L[i] = in_n_L[i];
    }
    adType xi[NUM_STATES];  		//  Initial value of the state for the next segment to be integrated
    adType xf[NUM_STATES];			//  Final value of the state for the last segment integrated
    adType h;

    double l_zero[3]={0.0,0.0,0.0};
    adType RigidSegmentLength;		// Length of the rigid segment - intermediate variable
//    adType tempadType;				// intermediate variables

    int   fsegno; 					// flexible segment no

    for (int i = 0; i < NUM_STATES; i++) {
        if (i < 3) {
            xi[i] = in_u[i];
        }
        else if (i < 12) {
            xi[i] = in_R[i - 3];
        }
        else if (i < 15) {
            xi[i] = in_p[i - 12];
        }
    }

    // we need to copy in_ftip to local variable
    adType ftip[3] = {0.0,0.0,0.0}; //placeholder
    // we will copy anything we will access more than once (or write to) to local variables
    int	  NextLocMarker=in_params.NextLocMarker;
    int	  InitialLocMarker=NextLocMarker;
    // for others, we will create aliases
    auto & SegBounds = in_params.SegBounds;
    auto & SegSteps = in_params.SegSteps;

    auto & InsertedLength = in_params.InsertedLength;
    auto & dlambdainv = in_params.dlambdainv;
    auto & K = in_params.K;
    auto & Kinv = in_params.Kinv;
    auto & ustar = in_params.ustar;
    auto & fcumlambda = in_params.fcumlambda;
    auto & FinalValueOnly = in_params.FinalValueOnly;
    auto & LocMarkers = in_params.LocMarkers;
    double p_atLocMarkers[NUM_LOCALIZATION_MARKERS][3];	// Positions at markers (ordered proximal to distal)
    // we will also copy the already filled entries of p_atLocMarkers from params
    for (int i = 0; i < InitialLocMarker; i++) {
        for (int j = 0; j < 3; j++) {
            p_atLocMarkers[i][j] = in_params.p_atLocMarkers[i][j];
        }
    }
    // Prepare the CRMIntegrand Parameters
    fsegno=SegmentIndex>>1; // i/2, flexible segment no
    h=(SegBounds[SegmentIndex+1]-SegBounds[SegmentIndex])/(SegSteps[fsegno]*1.0);

    // Integrate
    ABM4( 	xi, SegBounds[SegmentIndex], SegSteps[fsegno], h,
             InsertedLength, dlambdainv, K[fsegno], Kinv[fsegno], l_zero, ustar[fsegno], n_L,
             fcumlambda, ftip,
             FinalValueOnly, LocMarkers, &NextLocMarker,
             xf, p_atLocMarkers);

    if (!FinalValueOnly) {
        for (int i=0; i<NUM_LOCALIZATION_MARKERS; i++) {
            for (int j=0; j<3; j++) {
                out_p_atLocMarkers[i][j]=p_atLocMarkers[(NUM_LOCALIZATION_MARKERS-1)-i][j];
            }
        }
    }
    // pass the outputs
    for (int j=0; j < 3; j++) {
        out_p[j]=xf[3+9+j];
    }
    for (int j = 0; j < 9; ++j) {
        out_R[j] = xf[j+3];
    }
    mCopy_AB<3>(&(xf[0]), out_u);

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