/**
 * Python bindings for CRM catheter simulation.
 *
 * Wraps your existing C++ Cosserat rod solver for Python ML/RL use.
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "CRM.hpp"
#include "CRMDYN.hpp"
#include "CRM_BVPIVP_APIDeclarations.hpp"
#include "CRM_IVPJacobian.hpp"

namespace py = pybind11;

class ForwardKinematics
{
private:
    CRM::CatheterParams params_;

public:
    ForwardKinematics()
    {
        params_.LoadDefaultParameters();
    }

    py::dict compute(py::array_t<double> q, py::array_t<double> I)
    {
        auto q_buf = q.request();
        auto I_buf = I.request();

        if (I_buf.size != 6)
        {
            throw std::runtime_error("Expected 6 coil currents");
        }

        double *q_ptr = static_cast<double *>(q_buf.ptr);
        double *I_ptr = static_cast<double *>(I_buf.ptr);

        // Construct state vector for IVP solver
        int n_segments = q_buf.size / 3;
        Eigen::VectorXd state_vector(6 * n_segments + 6);
        state_vector.head(6).setZero();
        state_vector(2) = 1.0; // Initial orientation

        for (int i = 0; i < q_buf.size; ++i)
        {
            state_vector(6 + i) = q_ptr[i];
        }

        // Solve IVP
        CRM::IVPResult result;
        bool success = CRM::SolveIVP(params_, state_vector, I_ptr, result);

        if (!success)
        {
            throw std::runtime_error("Forward kinematics failed");
        }

        // Convert results to numpy arrays
        size_t n_points = result.positions.size();
        py::array_t<double> positions({n_points, 3});
        py::array_t<double> orientations({n_points, 9});

        auto pos_buf = positions.request();
        auto ori_buf = orientations.request();
        double *pos_ptr = static_cast<double *>(pos_buf.ptr);
        double *ori_ptr = static_cast<double *>(ori_buf.ptr);

        for (size_t i = 0; i < n_points; ++i)
        {
            // Position
            pos_ptr[i * 3 + 0] = result.positions[i](0);
            pos_ptr[i * 3 + 1] = result.positions[i](1);
            pos_ptr[i * 3 + 2] = result.positions[i](2);

            // Orientation (rotation matrix flattened)
            Eigen::Matrix3d R = result.rotations[i];
            for (int row = 0; row < 3; ++row)
            {
                for (int col = 0; col < 3; ++col)
                {
                    ori_ptr[i * 9 + row * 3 + col] = R(row, col);
                }
            }
        }

        py::dict output;
        output["positions"] = positions;
        output["orientations"] = orientations;
        output["success"] = success;

        return output;
    }

    py::array_t<double> get_tip_position(py::array_t<double> q, py::array_t<double> I)
    {
        auto result = compute(q, I);
        auto positions = result["positions"].cast<py::array_t<double>>();
        auto buf = positions.request();

        py::array_t<double> tip_pos(3);
        auto tip_buf = tip_pos.request();
        double *pos_ptr = static_cast<double *>(buf.ptr);
        double *tip_ptr = static_cast<double *>(tip_buf.ptr);

        size_t n_points = buf.shape[0];
        tip_ptr[0] = pos_ptr[(n_points - 1) * 3 + 0];
        tip_ptr[1] = pos_ptr[(n_points - 1) * 3 + 1];
        tip_ptr[2] = pos_ptr[(n_points - 1) * 3 + 2];

        return tip_pos;
    }

    void set_params(double length, double stiffness, double damping)
    {
        params_.catheter_length = length;
        params_.bending_stiffness = stiffness;
        params_.damping_coefficient = damping;
    }
};

class InverseKinematics
{
private:
    CRM::CatheterParams params_;

public:
    InverseKinematics()
    {
        params_.LoadDefaultParameters();
    }

    py::dict solve(
        py::array_t<double> target_position,
        py::object target_orientation,
        py::array_t<double> I,
        py::object q_init = py::none())
    {
        auto pos_buf = target_position.request();
        auto I_buf = I.request();

        if (pos_buf.size != 3 || I_buf.size != 6)
        {
            throw std::runtime_error("Invalid input dimensions");
        }

        // Target position
        Eigen::Vector3d target_pos;
        target_pos << *static_cast<double *>(pos_buf.ptr),
            *(static_cast<double *>(pos_buf.ptr) + 1),
            *(static_cast<double *>(pos_buf.ptr) + 2);

        double *I_ptr = static_cast<double *>(I_buf.ptr);

        // Solve BVP
        CRM::BVPResult result;
        bool success = CRM::SolveBVP(params_, target_pos, I_ptr, result);

        // Convert configuration to numpy
        py::array_t<double> q_result(result.q.size());
        auto q_buf = q_result.request();
        double *q_ptr = static_cast<double *>(q_buf.ptr);

        for (int i = 0; i < result.q.size(); ++i)
        {
            q_ptr[i] = result.q(i);
        }

        py::dict output;
        output["q"] = q_result;
        output["success"] = success;
        output["iterations"] = result.iterations;
        output["residual"] = result.residual_norm;

        return output;
    }
};

class DynamicsSimulator
{
private:
    CRMDYN::DynamicsParams params_;
    double dt_;

public:
    DynamicsSimulator(double dt = 0.001) : dt_(dt)
    {
        params_.LoadDefaultParameters();
    }

    py::array_t<double> step(py::array_t<double> state, py::array_t<double> control)
    {
        auto state_buf = state.request();
        auto ctrl_buf = control.request();

        if (ctrl_buf.size != 6)
        {
            throw std::runtime_error("Control must be 6D");
        }

        double *state_ptr = static_cast<double *>(state_buf.ptr);
        double *ctrl_ptr = static_cast<double *>(ctrl_buf.ptr);

        // Convert to Eigen
        Eigen::VectorXd current_state(state_buf.size);
        for (size_t i = 0; i < state_buf.size; ++i)
        {
            current_state(i) = state_ptr[i];
        }

        Eigen::VectorXd control_input(6);
        for (int i = 0; i < 6; ++i)
        {
            control_input(i) = ctrl_ptr[i];
        }

        // Integrate one step using your CRMDYN
        Eigen::VectorXd next_state = CRMDYN::IntegrateStep(
            current_state, control_input, dt_, params_);

        // Convert back to numpy
        py::array_t<double> result(next_state.size());
        auto result_buf = result.request();
        double *result_ptr = static_cast<double *>(result_buf.ptr);

        for (int i = 0; i < next_state.size(); ++i)
        {
            result_ptr[i] = next_state(i);
        }

        return result;
    }

    py::array_t<double> simulate(
        py::array_t<double> initial_state,
        py::array_t<double> control_sequence)
    {
        auto state_buf = initial_state.request();
        auto ctrl_buf = control_sequence.request();

        if (ctrl_buf.ndim != 2)
        {
            throw std::runtime_error("Control sequence must be 2D");
        }

        size_t n_steps = ctrl_buf.shape[0];
        size_t state_dim = state_buf.size;

        // Allocate trajectory storage
        py::array_t<double> trajectory({n_steps + 1, state_dim});
        auto traj_buf = trajectory.request();
        double *traj_ptr = static_cast<double *>(traj_buf.ptr);

        // Initial state
        double *state_ptr = static_cast<double *>(state_buf.ptr);
        std::copy(state_ptr, state_ptr + state_dim, traj_ptr);

        Eigen::VectorXd state(state_dim);
        for (size_t i = 0; i < state_dim; ++i)
        {
            state(i) = state_ptr[i];
        }

        double *ctrl_ptr = static_cast<double *>(ctrl_buf.ptr);

        // Simulate
        for (size_t t = 0; t < n_steps; ++t)
        {
            Eigen::VectorXd control(6);
            for (int j = 0; j < 6; ++j)
            {
                control(j) = ctrl_ptr[t * 6 + j];
            }

            state = CRMDYN::IntegrateStep(state, control, dt_, params_);

            for (size_t i = 0; i < state_dim; ++i)
            {
                traj_ptr[(t + 1) * state_dim + i] = state(i);
            }
        }

        return trajectory;
    }

    void set_timestep(double dt) { dt_ = dt; }
    double get_timestep() const { return dt_; }
};

class JacobianComputer
{
private:
    CRM::CatheterParams params_;

public:
    JacobianComputer()
    {
        params_.LoadDefaultParameters();
    }

    py::array_t<double> compute(py::array_t<double> q, py::array_t<double> I)
    {
        auto q_buf = q.request();
        auto I_buf = I.request();

        if (I_buf.size != 6)
        {
            throw std::runtime_error("Expected 6 coil currents");
        }

        double *q_ptr = static_cast<double *>(q_buf.ptr);
        double *I_ptr = static_cast<double *>(I_buf.ptr);
        int n_dof = q_buf.size;

        // Compute Jacobian using your implementation
        Eigen::MatrixXd J(6, n_dof);
        bool success = CRM::ComputeJacobian(params_, q_ptr, I_ptr, n_dof, J);

        if (!success)
        {
            throw std::runtime_error("Jacobian computation failed");
        }

        // Convert to numpy
        py::array_t<double> result({6, n_dof});
        auto result_buf = result.request();
        double *result_ptr = static_cast<double *>(result_buf.ptr);

        for (int i = 0; i < 6; ++i)
        {
            for (int j = 0; j < n_dof; ++j)
            {
                result_ptr[i * n_dof + j] = J(i, j);
            }
        }

        return result;
    }
};

PYBIND11_MODULE(crm_cpp, m)
{
    m.doc() = "C++ bindings for MRI-actuated catheter (CRM + CRMDYN)";

    py::class_<ForwardKinematics>(m, "ForwardKinematics")
        .def(py::init<>())
        .def("compute", &ForwardKinematics::compute,
             "Compute forward kinematics using IVP solver")
        .def("get_tip_position", &ForwardKinematics::get_tip_position,
             "Get tip position directly")
        .def("set_params", &ForwardKinematics::set_params,
             "Set catheter parameters");

    py::class_<InverseKinematics>(m, "InverseKinematics")
        .def(py::init<>())
        .def("solve", &InverseKinematics::solve,
             "Solve inverse kinematics using BVP solver");

    py::class_<DynamicsSimulator>(m, "DynamicsSimulator")
        .def(py::init<double>(), py::arg("dt") = 0.001)
        .def("step", &DynamicsSimulator::step,
             "Step dynamics forward by dt")
        .def("simulate", &DynamicsSimulator::simulate,
             "Simulate full trajectory")
        .def("set_timestep", &DynamicsSimulator::set_timestep)
        .def("get_timestep", &DynamicsSimulator::get_timestep);

    py::class_<JacobianComputer>(m, "JacobianComputer")
        .def(py::init<>())
        .def("compute", &JacobianComputer::compute,
             "Compute Jacobian matrix");

    m.attr("__version__") = "1.0.0";
}
