
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include "CRM.hpp"

namespace py = pybind11;
using namespace CRMCatheterModel;

// ============================================================
// Wrapper: Forward Kinematics
// ============================================================

Eigen::VectorXd forward_kinematics_wrapper(
    const Eigen::VectorXd &currents,
    double inserted_length,
    const std::string &config_file)
{
    // Load parameters from file
    CRMCatheterModelParams CathParams = Load_CRMCatheterModelParams(config_file.c_str());
    CatheterConfiguration CathConfig = Load_CatheterConfiguration(config_file.c_str());

    // Prepare input
    Eigen::VectorXd control_inputs(NUM_ACT_SET * 3 + 1);
    for (int i = 0; i < currents.size() - 1; i++)
    {
        control_inputs(i) = currents(i);
    }
    control_inputs(NUM_ACT_SET * 3) = inserted_length;

    // Prepare parameters
    CRMForwardKinematicsData FKParams;
    FKParams.CathParams = &CathParams;
    FKParams.CathConfig = &CathConfig;
    FKParams.ContactMode = ContactModeType::FREE_TIP;
    FKParams.FinalValueOnly = true;

    double dummy_energy = 0.0;
    int localmin = 0;

    // Allocate output arrays
    int output_dim = 15; // p[3] + R[9] + u[3]
    double *ReportedMarkerPos = new double[CathParams.no_locmarkers * 3];
    double *ReportedCoilPos = new double[NUM_ACT_SET * 3];
    double *ReportedCoilOrient = new double[NUM_ACT_SET * 9];

    FKParams.ReportedMarkerPos = (double (*)[3])ReportedMarkerPos;
    FKParams.ReportedCoilPos = (double (*)[3])ReportedCoilPos;
    FKParams.ReportedCoilOrient = (double (*)[9])ReportedCoilOrient;

    // Call forward kinematics
    Eigen::VectorXd result = CRM_ForwardKinematics(control_inputs, FKParams, dummy_energy, localmin);

    // Cleanup
    delete[] ReportedMarkerPos;
    delete[] ReportedCoilPos;
    delete[] ReportedCoilOrient;

    return result;
}

// ================================================
