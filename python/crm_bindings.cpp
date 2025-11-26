#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include "CRM.hpp"

namespace py = pybind11;
using namespace CRMCatheterModel;

// Simple wrapper for forward kinematics
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

PYBIND11_MODULE(crm_cpp, m)
{
    m.doc() = "MRI-Actuated Cosserat Rod Model - Python Bindings";

    // ============================================================
    // Expose Forward Kinematics
    // ============================================================
    m.def("forward_kinematics",
          &forward_kinematics_wrapper,
          py::arg("currents"),
          py::arg("inserted_length"),
          py::arg("config_file"),
          R"pbdoc(
            Compute forward kinematics of the catheter.
            
            Args:
                currents: (N,) array of actuation currents [A, A, A, ...]
                inserted_length: scalar inserted length [mm]
                config_file: path to catheter configuration file
                
            Returns:
                (15,) array [p_x, p_y, p_z, R_00, R_01, ..., R_22, u_x, u_y, u_z]
                where p is tip position (mm), R is 3x3 rotation matrix (row-major),
                u is tip curvature (1/mm)
        )pbdoc");

    // ============================================================
    // Expose ContactModeType enum
    // ============================================================
    py::enum_<ContactModeType>(m, "ContactMode")
        .value("FREE_TIP", ContactModeType::FREE_TIP)
        .value("FIXED_TIP", ContactModeType::FIXED_TIP)
        .export_values();

    // ============================================================
    // Expose utility function to check Eigen3 and constants
    // ============================================================
    m.def("get_num_actuators", []()
          { return NUM_ACT_SET; });
    m.def("get_num_states", []()
          { return NUM_STATES; });

    m.def("version_info", []()
          { return "CRM Catheter Model - Python Bindings v1.0\n"
                   "Cosserat Rod dynamics for MRI-actuated magnetic catheter\n"
                   "Compiled with pybind11"; });
}