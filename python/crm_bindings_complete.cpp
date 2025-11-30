/**
 * @file crm_bindings_complete.cpp
 * @brief Complete PyBind11 bindings for CRM Catheter library
 * 
 * This file provides Python bindings for:
 * - CRM core classes (CatheterClass, IVPSolverInputs, BVPSolverInputs)
 * - Forward kinematics and Jacobian computation
 * - IVP and BVP solvers
 * - Dynamics simulation
 * - State vector manipulation
 * - Matrix operations
 */

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>

// CRM Headers
#include "CRM.hpp"
#include "CRM_MatrixOperations.hpp"
#include "CRM_StateVector_Definitions.hpp"
#include "CRM_BVPIVP_APIDeclarations.hpp"
#include "CRM_IVPJacobian.hpp"
#include "CRMDYN.hpp"
#include "CRMDYN_Numerical_Integration.hpp"

// Numerical methods
#include "minpack.hpp"

namespace py = pybind11;

// ============================================================
// HELPER FUNCTIONS FOR PYTHON CONVERSION
// ============================================================

/**
 * @brief Convert Python list to Eigen VectorXd
 */
Eigen::VectorXd list_to_vector(const py::list& lst) {
    Eigen::VectorXd vec(lst.size());
    for (size_t i = 0; i < lst.size(); ++i) {
        vec(i) = lst[i].cast<double>();
    }
    return vec;
}

/**
 * @brief Convert Eigen VectorXd to Python list
 */
py::list vector_to_list(const Eigen::VectorXd& vec) {
    py::list lst;
    for (int i = 0; i < vec.size(); ++i) {
        lst.append(vec(i));
    }
    return lst;
}

// ============================================================
// PYBIND11 MODULE DEFINITION
// ============================================================

PYBIND11_MODULE(crm_cpp, m) {
    m.doc() = "Python bindings for CRM Catheter - Cosserat Rod Model";

    // ========================================================
    // CATHETER CLASS
    // ========================================================
    py::class_<CatheterClass>(m, "CatheterClass")
        .def(py::init<>(), "Default constructor")
        
        // Physical properties
        .def_readwrite("L", &CatheterClass::L, 
            "Total length of catheter (m)")
        .def_readwrite("d", &CatheterClass::d, 
            "Diameter of catheter (m)")
        .def_readwrite("E", &CatheterClass::E, 
            "Young's modulus (Pa)")
        .def_readwrite("G", &CatheterClass::G, 
            "Shear modulus (Pa)")
        .def_readwrite("rho", &CatheterClass::rho, 
            "Density (kg/m³)")
        
        // Geometric properties
        .def_readwrite("I", &CatheterClass::I, 
            "Second moment of area (m⁴)")
        .def_readwrite("A", &CatheterClass::A, 
            "Cross-sectional area (m²)")
        .def_readwrite("J", &CatheterClass::J, 
            "Polar moment of inertia (m⁴)")
        
        // Stiffness properties
        .def_readwrite("K_se", &CatheterClass::K_se, 
            "Shear/extension stiffness matrix (3x3)")
        .def_readwrite("K_bt", &CatheterClass::K_bt, 
            "Bending/torsion stiffness matrix (3x3)")
        
        // Magnetic properties
        .def_readwrite("m", &CatheterClass::m, 
            "Magnetic moment vector (3x1, A·m²)")
        .def_readwrite("B0", &CatheterClass::B0, 
            "Background magnetic field (3x1, T)")
        
        // Discretization
        .def_readwrite("N", &CatheterClass::N, 
            "Number of discretization segments")
        
        // Methods
        .def("compute_stiffness", &CatheterClass::compute_stiffness,
            "Compute stiffness matrices from E, G, I, A, J")
        
        .def("__repr__", [](const CatheterClass& c) {
            return "<CatheterClass L=" + std::to_string(c.L) + 
                   "m, N=" + std::to_string(c.N) + ">";
        });

    // ========================================================
    // IVP SOLVER INPUTS
    // ========================================================
    py::class_<IVPSolverInputs>(m, "IVPSolverInputs")
        .def(py::init<>(), "Default constructor")
        
        .def_readwrite("p0", &IVPSolverInputs::p0,
            "Initial position (3x1)")
        .def_readwrite("R0", &IVPSolverInputs::R0,
            "Initial rotation matrix (3x3)")
        .def_readwrite("n0", &IVPSolverInputs::n0,
            "Initial internal force (3x1, N)")
        .def_readwrite("m0", &IVPSolverInputs::m0,
            "Initial internal moment (3x1, N·m)")
        
        .def_readwrite("distributed_load", &IVPSolverInputs::distributed_load,
            "External distributed load (6x1)")
        
        .def_readwrite("magnetic_field", &IVPSolverInputs::magnetic_field,
            "Magnetic field function: s -> B(s)")
        .def_readwrite("magnetic_gradient", &IVPSolverInputs::magnetic_gradient,
            "Magnetic gradient function: s -> ∇B(s)")
        
        .def("__repr__", [](const IVPSolverInputs&) {
            return "<IVPSolverInputs>";
        });

    // ========================================================
    // BVP SOLVER INPUTS
    // ========================================================
    py::class_<BVPSolverInputs>(m, "BVPSolverInputs")
        .def(py::init<>(), "Default constructor")
        
        .def_readwrite("p0", &BVPSolverInputs::p0,
            "Initial position (3x1)")
        .def_readwrite("R0", &BVPSolverInputs::R0,
            "Initial rotation matrix (3x3)")
        
        .def_readwrite("tip_load", &BVPSolverInputs::tip_load,
            "Tip load (force and moment, 6x1)")
        
        .def_readwrite("magnetic_field", &BVPSolverInputs::magnetic_field,
            "Magnetic field function")
        .def_readwrite("magnetic_gradient", &BVPSolverInputs::magnetic_gradient,
            "Magnetic gradient function")
        
        .def("__repr__", [](const BVPSolverInputs&) {
            return "<BVPSolverInputs>";
        });

    // ========================================================
    // IVP SOLVER OUTPUTS
    // ========================================================
    py::class_<IVPSolverOutputs>(m, "IVPSolverOutputs")
        .def(py::init<>(), "Default constructor")
        
        .def_readwrite("p", &IVPSolverOutputs::p,
            "Position trajectory (3×N matrix)")
        .def_readwrite("R", &IVPSolverOutputs::R,
            "Rotation matrices (3×3×N tensor, stored as 9×N)")
        .def_readwrite("n", &IVPSolverOutputs::n,
            "Internal force (3×N)")
        .def_readwrite("m", &IVPSolverOutputs::m,
            "Internal moment (3×N)")
        .def_readwrite("v", &IVPSolverOutputs::v,
            "Linear strain (3×N)")
        .def_readwrite("u", &IVPSolverOutputs::u,
            "Angular strain (3×N)")
        
        .def_readwrite("s_vals", &IVPSolverOutputs::s_vals,
            "Arc length values (N×1)")
        
        .def("get_tip_position", [](const IVPSolverOutputs& out) {
            return out.p.col(out.p.cols() - 1);
        }, "Get tip position (last column of p)")
        
        .def("get_tip_rotation", [](const IVPSolverOutputs& out) {
            int N = out.R.cols();
            Eigen::Matrix3d R_tip;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    R_tip(i, j) = out.R(i * 3 + j, N - 1);
                }
            }
            return R_tip;
        }, "Get tip rotation matrix")
        
        .def("__repr__", [](const IVPSolverOutputs& out) {
            return "<IVPSolverOutputs N=" + std::to_string(out.s_vals.size()) + ">";
        });

    // ========================================================
    // BVP SOLVER OUTPUTS
    // ========================================================
    py::class_<BVPSolverOutputs>(m, "BVPSolverOutputs")
        .def(py::init<>(), "Default constructor")
        
        .def_readwrite("p", &BVPSolverOutputs::p,
            "Position trajectory (3×N)")
        .def_readwrite("R", &BVPSolverOutputs::R,
            "Rotation matrices (9×N)")
        .def_readwrite("n", &BVPSolverOutputs::n,
            "Internal force (3×N)")
        .def_readwrite("m", &BVPSolverOutputs::m,
            "Internal moment (3×N)")
        .def_readwrite("v", &BVPSolverOutputs::v,
            "Linear strain (3×N)")
        .def_readwrite("u", &BVPSolverOutputs::u,
            "Angular strain (3×N)")
        
        .def_readwrite("n0_solution", &BVPSolverOutputs::n0_solution,
            "Solved initial internal force (3×1)")
        .def_readwrite("m0_solution", &BVPSolverOutputs::m0_solution,
            "Solved initial internal moment (3×1)")
        
        .def_readwrite("residual", &BVPSolverOutputs::residual,
            "BVP solver residual")
        .def_readwrite("iterations", &BVPSolverOutputs::iterations,
            "Number of iterations")
        .def_readwrite("success", &BVPSolverOutputs::success,
            "Convergence flag")
        
        .def("__repr__", [](const BVPSolverOutputs& out) {
            return "<BVPSolverOutputs success=" + 
                   std::string(out.success ? "True" : "False") + ">";
        });

    // ========================================================
    // FORWARD KINEMATICS
    // ========================================================
    m.def("solve_ivp", 
        &CRM_IVP_Solver,
        py::arg("catheter"),
        py::arg("inputs"),
        R"doc(
        Solve Initial Value Problem (IVP) for catheter shape.
        
        Given initial conditions (p0, R0, n0, m0) and external loads,
        integrate the Cosserat rod equations forward along the arc length.
        
        Parameters
        ----------
        catheter : CatheterClass
            Catheter physical parameters
        inputs : IVPSolverInputs
            Initial conditions and loads
        
        Returns
        -------
        IVPSolverOutputs
            Complete catheter configuration (position, rotation, forces)
        
        Example
        -------
        >>> catheter = CatheterClass()
        >>> catheter.L = 0.5  # 50cm catheter
        >>> catheter.N = 100
        >>> 
        >>> inputs = IVPSolverInputs()
        >>> inputs.p0 = np.zeros(3)
        >>> inputs.R0 = np.eye(3)
        >>> inputs.n0 = np.array([0, 0, 0.01])  # 10mN compression
        >>> inputs.m0 = np.zeros(3)
        >>> 
        >>> outputs = solve_ivp(catheter, inputs)
        >>> tip_pos = outputs.get_tip_position()
        )doc");

    // ========================================================
    // BOUNDARY VALUE PROBLEM
    // ========================================================
    m.def("solve_bvp",
        &CRM_BVP_Solver,
        py::arg("catheter"),
        py::arg("inputs"),
        R"doc(
        Solve Boundary Value Problem (BVP) for catheter shape.
        
        Given boundary conditions (p0, R0 at base; tip load at tip),
        find the initial internal forces/moments that satisfy equilibrium.
        
        Parameters
        ----------
        catheter : CatheterClass
            Catheter physical parameters
        inputs : BVPSolverInputs
            Boundary conditions (base pose, tip load)
        
        Returns
        -------
        BVPSolverOutputs
            Complete solution including solved n0, m0
        
        Example
        -------
        >>> catheter = CatheterClass()
        >>> inputs = BVPSolverInputs()
        >>> inputs.p0 = np.zeros(3)
        >>> inputs.R0 = np.eye(3)
        >>> inputs.tip_load = np.array([0.01, 0, 0, 0, 0, 0])  # 10mN tip force
        >>> 
        >>> outputs = solve_bvp(catheter, inputs)
        >>> if outputs.success:
        >>>     print(f"Solved n0: {outputs.n0_solution}")
        )doc");

    // ========================================================
    // JACOBIAN COMPUTATION
    // ========================================================
    m.def("compute_jacobian",
        &CRM_Compute_Jacobian,
        py::arg("catheter"),
        py::arg("inputs"),
        py::arg("outputs"),
        R"doc(
        Compute analytical Jacobian: ∂(tip_pose)/∂(n0, m0)
        
        Uses variational approach to compute sensitivity of tip position
        and orientation with respect to initial internal forces/moments.
        
        Parameters
        ----------
        catheter : CatheterClass
            Catheter parameters
        inputs : IVPSolverInputs
            Current configuration inputs
        outputs : IVPSolverOutputs
            Current configuration outputs
        
        Returns
        -------
        Jacobian : ndarray (6×6)
            [∂p_tip/∂n0, ∂p_tip/∂m0; ∂R_tip/∂n0, ∂R_tip/∂m0]
        
        Example
        -------
        >>> J = compute_jacobian(catheter, inputs, outputs)
        >>> # J[0:3, 0:3] = ∂p_tip/∂n0
        >>> # J[3:6, 3:6] = ∂R_tip/∂m0 (rotation sensitivity)
        )doc");

    // ========================================================
    // DYNAMICS STRUCTURES
    // ========================================================
    py::class_<DynamicsInputs>(m, "DynamicsInputs")
        .def(py::init<>(), "Default constructor")
        
        .def_readwrite("state", &DynamicsInputs::state,
            "Current state vector: [p; R_vec; n; m; p_dot; omega]")
        
        .def_readwrite("control", &DynamicsInputs::control,
            "Control inputs (e.g., coil currents)")
        
        .def_readwrite("t", &DynamicsInputs::t,
            "Current time (s)")
        
        .def_readwrite("external_wrench", &DynamicsInputs::external_wrench,
            "External wrench at tip (6×1)")
        
        .def("__repr__", [](const DynamicsInputs& d) {
            return "<DynamicsInputs t=" + std::to_string(d.t) + "s>";
        });

    py::class_<DynamicsOutputs>(m, "DynamicsOutputs")
        .def(py::init<>(), "Default constructor")
        
        .def_readwrite("state_dot", &DynamicsOutputs::state_dot,
            "Time derivative of state vector")
        
        .def_readwrite("p", &DynamicsOutputs::p,
            "Position trajectory (3×N)")
        .def_readwrite("R", &DynamicsOutputs::R,
            "Rotation trajectory (9×N)")
        
        .def_readwrite("kinetic_energy", &DynamicsOutputs::kinetic_energy,
            "Total kinetic energy (J)")
        .def_readwrite("potential_energy", &DynamicsOutputs::potential_energy,
            "Total potential energy (J)")
        
        .def("__repr__", [](const DynamicsOutputs&) {
            return "<DynamicsOutputs>";
        });

    // ========================================================
    // DYNAMICS SOLVER
    // ========================================================
    m.def("compute_dynamics",
        &CRMDYN_ComputeDynamics,
        py::arg("catheter"),
        py::arg("inputs"),
        R"doc(
        Compute dynamics: state_dot = f(state, control, t)
        
        Evaluates the continuous-time dynamics of the catheter system
        including inertial effects, elastic forces, and external loads.
        
        Parameters
        ----------
        catheter : CatheterClass
            Catheter parameters including mass distribution
        inputs : DynamicsInputs
            Current state, control, time
        
        Returns
        -------
        DynamicsOutputs
            Time derivatives and energetics
        
        Notes
        -----
        Use with ODE integrator (RK4, etc.) for time-stepping simulation.
        
        Example
        -------
        >>> dyn_in = DynamicsInputs()
        >>> dyn_in.state = initial_state  # [p; R; n; m; v; omega]
        >>> dyn_in.control = coil_currents
        >>> dyn_in.t = 0.0
        >>> 
        >>> dyn_out = compute_dynamics(catheter, dyn_in)
        >>> state_dot = dyn_out.state_dot
        >>> # Integrate: state_{k+1} = state_k + dt * state_dot
        )doc");

    // ========================================================
    // MATRIX OPERATIONS
    // ========================================================
    m.def("hat", 
        &CRM_hat,
        py::arg("v"),
        R"doc(
        Skew-symmetric matrix from 3D vector (hat operator).
        
        Maps v ∈ ℝ³ to [v]× ∈ so(3) such that [v]× w = v × w
        
        Parameters
        ----------
        v : ndarray (3,)
            3D vector
        
        Returns
        -------
        V : ndarray (3, 3)
            Skew-symmetric matrix
        
        Example
        -------
        >>> v = np.array([1, 2, 3])
        >>> V = hat(v)
        >>> # V @ w == np.cross(v, w) for any w
        )doc");

    m.def("vee",
        &CRM_vee,
        py::arg("V"),
        R"doc(
        Extract vector from skew-symmetric matrix (vee operator).
        
        Inverse of hat: so(3) → ℝ³
        
        Parameters
        ----------
        V : ndarray (3, 3)
            Skew-symmetric matrix
        
        Returns
        -------
        v : ndarray (3,)
            Extracted vector
        )doc");

    m.def("adjoint",
        &CRM_Adjoint,
        py::arg("g"),
        R"doc(
        Compute adjoint representation of SE(3) element.
        
        For g = (R, p) ∈ SE(3), computes Ad_g ∈ ℝ⁶ˣ⁶
        
        Parameters
        ----------
        g : tuple (R, p)
            R: rotation matrix (3×3), p: position (3×1)
        
        Returns
        -------
        Ad_g : ndarray (6, 6)
            Adjoint matrix
        )doc");

    // ========================================================
    // STATE VECTOR MANIPULATION
    // ========================================================
    m.def("pack_state",
        [](const Eigen::VectorXd& p, const Eigen::MatrixXd& R,
           const Eigen::VectorXd& n, const Eigen::VectorXd& m,
           const Eigen::VectorXd& v, const Eigen::VectorXd& omega) {
            return CRM_PackStateVector(p, R, n, m, v, omega);
        },
        py::arg("p"), py::arg("R"), py::arg("n"), 
        py::arg("m"), py::arg("v"), py::arg("omega"),
        R"doc(
        Pack kinematic/dynamic variables into state vector.
        
        State = [p; vec(R); n; m; v; omega]
        
        Parameters
        ----------
        p : ndarray (3×N)
            Position
        R : ndarray (3×3×N)
            Rotation matrices (as 9×N)
        n, m : ndarray (3×N)
            Internal forces/moments
        v, omega : ndarray (3×N)
            Linear/angular velocities
        
        Returns
        -------
        state : ndarray (dim,)
            Packed state vector
        )doc");

    m.def("unpack_state",
        &CRM_UnpackStateVector,
        py::arg("state"),
        py::arg("N"),
        R"doc(
        Unpack state vector into individual components.
        
        Returns
        -------
        tuple : (p, R, n, m, v, omega)
        )doc");

    // ========================================================
    // UTILITY FUNCTIONS
    // ========================================================
    m.def("rotation_to_axis_angle",
        [](const Eigen::Matrix3d& R) {
            // Compute axis-angle representation
            double angle = std::acos((R.trace() - 1.0) / 2.0);
            Eigen::Vector3d axis;
            if (std::abs(angle) < 1e-10) {
                axis = Eigen::Vector3d::Zero();
            } else {
                axis << R(2,1) - R(1,2), 
                        R(0,2) - R(2,0), 
                        R(1,0) - R(0,1);
                axis = axis / (2.0 * std::sin(angle));
            }
            return std::make_pair(axis, angle);
        },
        py::arg("R"),
        "Convert rotation matrix to axis-angle representation");

    m.def("axis_angle_to_rotation",
        [](const Eigen::Vector3d& axis, double angle) {
            // Rodrigues formula
            Eigen::Matrix3d K;
            K << 0, -axis(2), axis(1),
                 axis(2), 0, -axis(0),
                 -axis(1), axis(0), 0;
            
            Eigen::Matrix3d R = Eigen::Matrix3d::Identity() +
                                std::sin(angle) * K +
                                (1.0 - std::cos(angle)) * K * K;
            return R;
        },
        py::arg("axis"), py::arg("angle"),
        "Convert axis-angle to rotation matrix");

    m.def("create_default_catheter",
        []() {
            CatheterClass catheter;
            
            // Typical values for MRI-actuated catheter
            catheter.L = 0.5;      // 50 cm
            catheter.d = 0.003;    // 3 mm diameter
            catheter.E = 1e6;      // 1 MPa (soft polymer)
            catheter.G = catheter.E / (2.0 * (1.0 + 0.3)); // ν ≈ 0.3
            catheter.rho = 1200;   // kg/m³ (polymer density)
            catheter.N = 100;      // 100 segments
            
            // Compute geometric properties
            catheter.A = M_PI * catheter.d * catheter.d / 4.0;
            catheter.I = M_PI * std::pow(catheter.d, 4) / 64.0;
            catheter.J = 2.0 * catheter.I;
            
            // Compute stiffness
            catheter.K_se = Eigen::Matrix3d::Identity() * catheter.G * catheter.A;
            catheter.K_se(2, 2) = catheter.E * catheter.A; // Axial stiffness
            
            catheter.K_bt = Eigen::Matrix3d::Identity() * catheter.E * catheter.I;
            catheter.K_bt(2, 2) = catheter.G * catheter.J; // Torsional stiffness
            
            // Default magnetic properties
            catheter.m = Eigen::Vector3d(0, 0, 1e-6); // 1 μA·m²
            catheter.B0 = Eigen::Vector3d(0, 0, 1.5);  // 1.5T (typical MRI)
            
            return catheter;
        },
        "Create catheter with typical default parameters");

    // ========================================================
    // VERSION INFO
    // ========================================================
    m.attr("__version__") = "1.0.0";
    m.attr("__author__") = "CRM Catheter Development Team";
    
    // ========================================================
    // CONSTANTS
    // ========================================================
    m.attr("PI") = M_PI;
    m.attr("GRAVITY") = 9.81;  // m/s²
}
