#include <torch/extension.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <array>
#include <stdexcept>
#include <string>
#include <vector>

#include "CRM.hpp"
#include "CRM_BVPIVP_APIDeclarations.hpp"
#include "CRM_IVPJacobian.hpp"

namespace py = pybind11;

namespace CRMCatheterModel {

	struct EquilibriumConfig {
		CRMCatheterModelParams cath_params;
		CatheterConfiguration cath_config;
		ContactModeType contact_mode;
		double tip_force[3];
		double tip_constraint_point[3];
		double integration_step_size;
		double residual_threshold;

		EquilibriumConfig(const std::string& cath_params_path,
			const std::string& cath_config_path,
			double integration_step,
			const std::vector<double>& tip_force_in,
			double residual_threshold_in)
			: cath_params(Load_CRMCatheterModelParams(cath_params_path.c_str())),
			  cath_config(Load_CatheterConfiguration(cath_config_path.c_str())),
			  contact_mode(ContactModeType::FREE_TIP),
			  integration_step_size(integration_step),
			  residual_threshold(residual_threshold_in) {
			if (tip_force_in.size() != 3) {
				throw std::runtime_error("tip_force must have 3 elements.");
			}
			for (int i = 0; i < 3; i++) {
				tip_force[i] = tip_force_in[i];
				tip_constraint_point[i] = 0.0;
			}
		}
	};

	static void CheckInputs(const torch::Tensor& u, const torch::Tensor& Li) {
		TORCH_CHECK(u.device().is_cpu(), "u must be a CPU tensor.");
		TORCH_CHECK(Li.device().is_cpu(), "Li must be a CPU tensor.");
		TORCH_CHECK(u.dtype() == torch::kFloat64, "u must be float64.");
		TORCH_CHECK(Li.dtype() == torch::kFloat64, "Li must be float64.");
		TORCH_CHECK(u.dim() == 3, "u must have shape [B, N_act, 3].");
		TORCH_CHECK(u.size(2) == 3, "u must have shape [B, N_act, 3].");
		TORCH_CHECK(u.size(1) == NUM_ACT_SET, "u.size(1) must equal NUM_ACT_SET.");
		TORCH_CHECK(NUM_ACT_SET == 1, "v1.0 requires NUM_ACT_SET == 1.");
		TORCH_CHECK(Li.numel() == 1 || Li.numel() == u.size(0), "Li must be scalar or shape [B].");
		TORCH_CHECK(!Li.requires_grad(), "Li must have requires_grad=False.");
	}

	static CRMShootingMethodParams BuildShootingParams(const EquilibriumConfig& cfg,
		double Li_val,
		double actuation_currents[NUM_ACT_SET][3]) {
		double tip_constraint_point[3];
		double tip_force[3];
		for (int i = 0; i < 3; i++) {
			tip_constraint_point[i] = cfg.tip_constraint_point[i];
			tip_force[i] = cfg.tip_force[i];
		}
		return CRMConstructShootingMethodParamSet(
			cfg.cath_params,
			cfg.cath_config,
			Li_val,
			actuation_currents,
			cfg.contact_mode,
			tip_constraint_point,
			tip_force,
			cfg.integration_step_size);
	}

	std::vector<torch::Tensor> crm_equilibrium_forward(torch::Tensor u, torch::Tensor Li, std::shared_ptr<EquilibriumConfig> cfg) {
		TORCH_CHECK(cfg != nullptr, "cfg_equil is null.");
		CheckInputs(u, Li);

		u = u.contiguous();
		Li = Li.contiguous();

		const auto B = u.size(0);
		const bool Li_is_scalar = (Li.numel() == 1);
		const double Li_scalar = Li_is_scalar ? Li.item<double>() : 0.0;
		auto u_acc = u.accessor<double, 3>();
		const double* Li_ptr = Li_is_scalar ? nullptr : Li.data_ptr<double>();

		auto options = u.options();
		auto p_tip = torch::zeros({ B, 3 }, options);
		auto deltau0 = torch::zeros({ B, 3 }, options);
		auto ftip = torch::zeros({ B, 3 }, options);
		auto localmin = torch::zeros({ B }, torch::TensorOptions().dtype(torch::kInt32).device(torch::kCPU));
		auto residual_norm = torch::zeros({ B }, options);

		auto p_tip_acc = p_tip.accessor<double, 2>();
		auto deltau0_acc = deltau0.accessor<double, 2>();
		auto ftip_acc = ftip.accessor<double, 2>();
		auto localmin_acc = localmin.accessor<int32_t, 1>();
		auto residual_acc = residual_norm.accessor<double, 1>();

		for (int64_t b = 0; b < B; b++) {
			double act_currents[NUM_ACT_SET][3];
			for (int i = 0; i < NUM_ACT_SET; i++) {
				for (int j = 0; j < 3; j++) {
					act_currents[i][j] = u_acc[b][i][j];
				}
			}

			double Li_val = Li_is_scalar ? Li_scalar : Li_ptr[b];
			CRMShootingMethodParams params = BuildShootingParams(*cfg, Li_val, act_currents);

			double deltau0_guess[3] = { 0.0, 0.0, 0.0 };
			double ftip_guess[3] = { 0.0, 0.0, 0.0 };
			double deltau0_out[3];
			double ftip_out[3];
			int localmin_out = 0;

			CRMShootingMethodBVP(params, deltau0_guess, ftip_guess, deltau0_out, ftip_out, localmin_out);

			double xf[NUM_STATES];
			double moment_residual[3];
			double dummyPE = 0.0;
			std::vector<double> out_p_atLocMarkers(params.no_locmarkers * 3);
			std::vector<double> out_R_atActuators(params.no_act_set * 9);
			std::vector<double> out_p_atActuators(params.no_act_set * 3);

			CRMSolverIVP(
				params,
				deltau0_out,
				ftip_out,
				false,
				true,
				xf,
				moment_residual,
				dummyPE,
				reinterpret_cast<double(*)[3]>(out_p_atLocMarkers.data()),
				reinterpret_cast<double(*)[9]>(out_R_atActuators.data()),
				reinterpret_cast<double(*)[3]>(out_p_atActuators.data()));

			Eigen::Matrix<double, 3, 1> F_scaled = CRM_EquilibriumResidual_FScaled(params, deltau0_out);
			double res_norm = std::sqrt(F_scaled.squaredNorm());

			for (int i = 0; i < 3; i++) {
				p_tip_acc[b][i] = xf[i];
				deltau0_acc[b][i] = deltau0_out[i];
				ftip_acc[b][i] = ftip_out[i];
			}
			localmin_acc[b] = static_cast<int32_t>(localmin_out);
			residual_acc[b] = res_norm;
		}

		return { p_tip, deltau0, ftip, localmin, residual_norm };
	}

	torch::Tensor crm_equilibrium_backward(
		torch::Tensor grad_p,
		torch::Tensor u,
		torch::Tensor Li,
		torch::Tensor deltau0,
		torch::Tensor ftip,
		torch::Tensor localmin,
		torch::Tensor residual_norm,
		std::shared_ptr<EquilibriumConfig> cfg) {

		TORCH_CHECK(cfg != nullptr, "cfg_equil is null.");
		CheckInputs(u, Li);
		TORCH_CHECK(grad_p.device().is_cpu(), "grad_p must be a CPU tensor.");
		TORCH_CHECK(grad_p.dtype() == torch::kFloat64, "grad_p must be float64.");
		TORCH_CHECK(grad_p.dim() == 2 && grad_p.size(1) == 3, "grad_p must have shape [B, 3].");

		u = u.contiguous();
		Li = Li.contiguous();
		grad_p = grad_p.contiguous();
		deltau0 = deltau0.contiguous();
		ftip = ftip.contiguous();
		localmin = localmin.contiguous();
		residual_norm = residual_norm.contiguous();

		const auto B = u.size(0);
		const bool Li_is_scalar = (Li.numel() == 1);
		const double Li_scalar = Li_is_scalar ? Li.item<double>() : 0.0;
		auto u_acc = u.accessor<double, 3>();
		const double* Li_ptr = Li_is_scalar ? nullptr : Li.data_ptr<double>();
		auto grad_p_acc = grad_p.accessor<double, 2>();
		auto deltau0_acc = deltau0.accessor<double, 2>();
		auto ftip_acc = ftip.accessor<double, 2>();
		auto localmin_acc = localmin.accessor<int32_t, 1>();
		auto residual_acc = residual_norm.accessor<double, 1>();

		auto grad_u = torch::zeros_like(u);
		auto grad_u_acc = grad_u.accessor<double, 3>();

		for (int64_t b = 0; b < B; b++) {
			const int32_t localmin_b = localmin_acc[b];
			const double res_norm_b = residual_acc[b];
			if (localmin_b != 0 || res_norm_b > cfg->residual_threshold) {
				TORCH_CHECK(false,
					"[crm_equilibrium_backward] non-convergence: batch=",
					b,
					" step=-1 solver_exit_flag=",
					localmin_b,
					" residual_norm=",
					res_norm_b);
			}

			double act_currents[NUM_ACT_SET][3];
			for (int i = 0; i < NUM_ACT_SET; i++) {
				for (int j = 0; j < 3; j++) {
					act_currents[i][j] = u_acc[b][i][j];
				}
			}

			double Li_val = Li_is_scalar ? Li_scalar : Li_ptr[b];
			CRMShootingMethodParams params = BuildShootingParams(*cfg, Li_val, act_currents);

			double deltau0_b[3];
			double ftip_b[3];
			for (int i = 0; i < 3; i++) {
				deltau0_b[i] = deltau0_acc[b][i];
				ftip_b[i] = ftip_acc[b][i];
			}

			const Eigen::Matrix<double, 3, 3, Eigen::RowMajor> A = CRM_EquilibriumResidualJacobian_A(params, deltau0_b);
			const Eigen::Matrix<double, 3, CURRENT_ACT_VECTOR_DIM, Eigen::RowMajor> Bmat =
				CRM_EquilibriumResidualJacobian_B(params, deltau0_b, ftip_b);
			const CRMIVPTipJacobiansRaw Jraw = CRMSolverIVPJacobian_TipRaw(params, deltau0_b, ftip_b);

			const Eigen::Matrix<double, 3, 3, Eigen::RowMajor> y_x = Jraw.p_u0 * IVALUE_SCALE_DU;
			const Eigen::Matrix<double, 3, CURRENT_ACT_VECTOR_DIM, Eigen::RowMajor> y_u_direct = Jraw.p_zc;

			Eigen::Matrix<double, 3, 1> g_p;
			g_p << grad_p_acc[b][0], grad_p_acc[b][1], grad_p_acc[b][2];

			const Eigen::Matrix<double, 3, 1> rhs = y_x.transpose() * g_p;
			Eigen::FullPivLU<Eigen::Matrix<double, 3, 3>> lu(A.transpose());
			TORCH_CHECK(lu.isInvertible(),
				"[crm_equilibrium_backward] A^T is singular: batch=",
				b,
				" step=-1");

			const Eigen::Matrix<double, 3, 1> lambda = lu.solve(rhs);
			const Eigen::Matrix<double, CURRENT_ACT_VECTOR_DIM, 1> g_u_vec =
				y_u_direct.transpose() * g_p - Bmat.transpose() * lambda;

			for (int j = 0; j < 3; j++) {
				grad_u_acc[b][0][j] = g_u_vec(j);
			}
		}

		return grad_u;
	}

}  // namespace CRMCatheterModel

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {
	py::class_<CRMCatheterModel::EquilibriumConfig, std::shared_ptr<CRMCatheterModel::EquilibriumConfig>>(m, "EquilibriumConfig")
		.def(py::init<const std::string&, const std::string&, double, const std::vector<double>&, double>(),
			py::arg("cath_params_path"),
			py::arg("cath_config_path"),
			py::arg("integration_step"),
			py::arg("tip_force"),
			py::arg("residual_threshold") = CRMCatheterModel::TRUSTREGION_TOLERANCE);

	m.def("crm_equilibrium_forward", &CRMCatheterModel::crm_equilibrium_forward, "CRM equilibrium forward (TIP_POSITION_ONLY)");
	m.def("crm_equilibrium_backward", &CRMCatheterModel::crm_equilibrium_backward, "CRM equilibrium backward (TIP_POSITION_ONLY)");
}
