#include <torch/extension.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <algorithm>
#include <array>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "CRM.hpp"
#include "CRMDYN.hpp"

namespace py = pybind11;

namespace CRMCatheterModel {

	static bool FDReferenceEnabled() {
		const char* flag = std::getenv("CRM_DIFFSIMS_ENABLE_FD_REFERENCE");
		return flag != nullptr && std::string(flag) == "1";
	}

	struct DynamicsConfig {
		CRMCatheterModelParams cath_params;
		CatheterConfiguration cath_config;
		ContactModeType contact_mode;
		double tip_force[3];
		double tip_constraint_point[3];
		double integration_step_size;
		double act_inertia[NUM_ACT_SET][9];
		double damping[NUM_ACT_SET][6];
		double delta_t;
		double residual_threshold;
		double rot_tol;

		DynamicsConfig(const std::string& cath_params_path,
			const std::string& cath_config_path,
			double integration_step,
			const std::vector<double>& tip_force_in,
			const std::vector<double>& act_inertia_in,
			const std::vector<double>& damping_in,
			double delta_t_in,
			double residual_threshold_in,
			double rot_tol_in)
			: cath_params(Load_CRMCatheterModelParams(cath_params_path.c_str())),
			  cath_config(Load_CatheterConfiguration(cath_config_path.c_str())),
			  contact_mode(ContactModeType::FREE_TIP),
			  integration_step_size(integration_step),
			  delta_t(delta_t_in),
			  residual_threshold(residual_threshold_in),
			  rot_tol(rot_tol_in) {
			if (tip_force_in.size() != 3) {
				throw std::runtime_error("tip_force must have 3 elements.");
			}
			if (act_inertia_in.size() != NUM_ACT_SET * 9) {
				throw std::runtime_error("act_inertia must have NUM_ACT_SET*9 elements.");
			}
			if (damping_in.size() != NUM_ACT_SET * 6) {
				throw std::runtime_error("damping must have NUM_ACT_SET*6 elements.");
			}
			for (int i = 0; i < 3; i++) {
				tip_force[i] = tip_force_in[i];
				tip_constraint_point[i] = 0.0;
			}
			for (int act = 0; act < NUM_ACT_SET; act++) {
				for (int i = 0; i < 9; i++) {
					act_inertia[act][i] = act_inertia_in[act * 9 + i];
				}
				for (int i = 0; i < 6; i++) {
					damping[act][i] = damping_in[act * 6 + i];
				}
			}
		}
	};

	static void CheckStepInputs(const torch::Tensor& x_t, const torch::Tensor& u_t, const torch::Tensor& Li) {
		TORCH_CHECK(x_t.device().is_cpu(), "x_t must be a CPU tensor.");
		TORCH_CHECK(u_t.device().is_cpu(), "u_t must be a CPU tensor.");
		TORCH_CHECK(Li.device().is_cpu(), "Li must be a CPU tensor.");
		TORCH_CHECK(x_t.dtype() == torch::kFloat64, "x_t must be float64.");
		TORCH_CHECK(u_t.dtype() == torch::kFloat64, "u_t must be float64.");
		TORCH_CHECK(Li.dtype() == torch::kFloat64, "Li must be float64.");
		TORCH_CHECK(x_t.dim() == 2, "x_t must have shape [B, 24*N_act+15].");
		TORCH_CHECK(u_t.dim() == 3 && u_t.size(2) == 3, "u_t must have shape [B, N_act, 3].");
		TORCH_CHECK(u_t.size(1) == NUM_ACT_SET, "u_t.size(1) must equal NUM_ACT_SET.");
		TORCH_CHECK(x_t.size(1) == 24 * NUM_ACT_SET + 15, "x_t has incorrect packed dimension.");
		TORCH_CHECK(NUM_ACT_SET == 1, "v1.1 requires NUM_ACT_SET == 1.");
		TORCH_CHECK(Li.numel() == 1 || Li.numel() == x_t.size(0), "Li must be scalar or shape [B].");
		TORCH_CHECK(!Li.requires_grad(), "Li must have requires_grad=False.");
	}

	static bool CheckRotationMatrix(const double R[9], double tol) {
		double RtR[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		for (int r = 0; r < 3; r++) {
			for (int c = 0; c < 3; c++) {
				double sum = 0.0;
				for (int k = 0; k < 3; k++) {
					sum += R[k * 3 + r] * R[k * 3 + c];
				}
				RtR[r * 3 + c] = sum;
			}
		}
		double diff_sq = 0.0;
		for (int i = 0; i < 9; i++) {
			const double target = (i % 4 == 0) ? 1.0 : 0.0;
			const double d = RtR[i] - target;
			diff_sq += d * d;
		}
		const double frob = std::sqrt(diff_sq);
		const double det =
			R[0] * (R[4] * R[8] - R[5] * R[7]) -
			R[1] * (R[3] * R[8] - R[5] * R[6]) +
			R[2] * (R[3] * R[7] - R[4] * R[6]);
		return std::isfinite(frob) && (frob < tol) && (det > 0.0);
	}

	static void EvaluateDynResidual(CRMShootingMethodParams& params,
		const double xf_pre[NUM_STATES],
		double mL_initialguess[NUM_ACT_SET][3],
		double nL_initialguess[NUM_ACT_SET][3],
		double eta_scaled[NUM_DYN_RESIDUAL],
		double out_residual[NUM_DYN_RESIDUAL],
		double out_u0[3],
		double out_tau[NUM_ACT_SET][3]) {

		double x_0[NUM_STATES];
		for (int i = 0; i < NUM_STATES; i++) {
			if (i < 3) x_0[i] = params.p0[i];
			else if (i < 12) x_0[i] = params.R0[i - 3];
			else x_0[i] = 0.0;
		}

		DYNNLEqnParams dyn_params(params.no_flex_seg, params.no_rigid_seg, params.no_act_set, params.no_locmarkers, params.no_fcum_steps);
		bool FinalValueOnly = true;
		CRMDYNSolverIVP_Prep(params.no_flex_seg, params.no_rigid_seg, params.no_act_set, params.no_locmarkers, params.no_fcum_steps,
			x_0, params.IntegrationStepSize,
			params.Li, params.dlambdainv, params.rho, params.SegmentTypes,
			params.SegEndLambdas, params.LocMarkerLambdas,
			params.K, params.Kinv, params.ustar,
			params.MagMoment, params.fcumlambda, params.CoilAlignmentTurnAreaMatrix,
			params.B0, params.g, params.ActMass, params.actInertia, params.damping, params.DELTA_T,
			params.v_L_pre, params.w_L_pre, params.p_pre, params.R_pre,
			mL_initialguess, nL_initialguess, FinalValueOnly, dyn_params);

		dyn_params.ContactMode = params.ContactMode;
		mCopy_AB<3>(params.TipConstraintPoint, dyn_params.TipConstraintPoint);
		for (int i = 0; i < 3; i++) {
			dyn_params.TipForce[i] = params.TipForce[i];
			dyn_params.ftip_initialguess[i] = 0.0;
		}
		for (int i = 0; i < NUM_STATES; i++) {
			dyn_params.xf[i] = xf_pre[i];
		}

		double out_tau_flat[NUM_ACT_SET * 3];
		DYNNLEquation(eta_scaled, out_residual, dyn_params, out_u0, out_tau_flat);

		for (int act = 0; act < NUM_ACT_SET; act++) {
			for (int i = 0; i < 3; i++) {
				out_tau[act][i] = out_tau_flat[act * 3 + i];
			}
		}

	}

	static void EvaluateDynHOutput(CRMShootingMethodParams& params,
		const double xf_pre[NUM_STATES],
		double eta_scaled[NUM_DYN_RESIDUAL],
		double mL_initialguess[NUM_ACT_SET][3],
		double nL_initialguess[NUM_ACT_SET][3],
		double out_x_tp1[24 * NUM_ACT_SET + 15]) {

		double u0_tmp[3];
		double tau_tmp[NUM_ACT_SET][3];
		double residual_tmp[NUM_DYN_RESIDUAL];
		EvaluateDynResidual(params, xf_pre, mL_initialguess, nL_initialguess, eta_scaled, residual_tmp, u0_tmp, tau_tmp);

		double out_mL[NUM_ACT_SET][3];
		double out_nL[NUM_ACT_SET][3];
		for (int act = 0; act < NUM_ACT_SET; act++) {
			for (int i = 0; i < 3; i++) {
				out_mL[act][i] = IVALUE_SCALE_M * eta_scaled[act * 6 + i];
				out_nL[act][i] = IVALUE_SCALE_N * eta_scaled[act * 6 + 3 + i];
			}
		}

		double out_xf[NUM_STATES];
		double out_coil_state[NUM_ACT_SET][NUM_COIL_STATES];
		std::vector<double> out_p_atLocMarkers(params.no_locmarkers * 3);
		double ftip[3];
		for (int i = 0; i < 3; i++) {
			ftip[i] = params.TipForce[i];
		}

		DYNSolverIVP(params, u0_tmp, out_mL, out_nL, tau_tmp, ftip,
			true, out_xf, out_coil_state, reinterpret_cast<double(*)[3]>(out_p_atLocMarkers.data()));

		int64_t out_offset = 0;
		for (int act = 0; act < NUM_ACT_SET; act++) {
			for (int i = 0; i < 3; i++) out_x_tp1[out_offset + act * 3 + i] = out_coil_state[act][i];
		}
		out_offset += 3 * NUM_ACT_SET;
		for (int act = 0; act < NUM_ACT_SET; act++) {
			for (int i = 0; i < 3; i++) out_x_tp1[out_offset + act * 3 + i] = out_coil_state[act][i + 3];
		}
		out_offset += 3 * NUM_ACT_SET;
		for (int act = 0; act < NUM_ACT_SET; act++) {
			for (int i = 0; i < 3; i++) out_x_tp1[out_offset + act * 3 + i] = out_mL[act][i];
		}
		out_offset += 3 * NUM_ACT_SET;
		for (int act = 0; act < NUM_ACT_SET; act++) {
			for (int i = 0; i < 3; i++) out_x_tp1[out_offset + act * 3 + i] = out_nL[act][i];
		}
		out_offset += 3 * NUM_ACT_SET;
		for (int act = 0; act < NUM_ACT_SET; act++) {
			for (int i = 0; i < 3; i++) out_x_tp1[out_offset + act * 3 + i] = out_coil_state[act][i + 6];
		}
		out_offset += 3 * NUM_ACT_SET;
		for (int act = 0; act < NUM_ACT_SET; act++) {
			for (int i = 0; i < 9; i++) out_x_tp1[out_offset + act * 9 + i] = out_coil_state[act][i + 9];
		}
		out_offset += 9 * NUM_ACT_SET;
		for (int i = 0; i < NUM_STATES; i++) out_x_tp1[out_offset + i] = out_xf[i];
	}

	std::vector<torch::Tensor> crm_step_forward(torch::Tensor x_t, torch::Tensor u_t, torch::Tensor Li, std::shared_ptr<DynamicsConfig> cfg) {
		TORCH_CHECK(cfg != nullptr, "cfg_dyn is null.");
		CheckStepInputs(x_t, u_t, Li);

		x_t = x_t.contiguous();
		u_t = u_t.contiguous();
		Li = Li.contiguous();

		const auto B = x_t.size(0);
		const bool Li_is_scalar = (Li.numel() == 1);
		const double Li_scalar = Li_is_scalar ? Li.item<double>() : 0.0;
		const double* Li_ptr = Li_is_scalar ? nullptr : Li.data_ptr<double>();

		auto x_acc = x_t.accessor<double, 2>();
		auto u_acc = u_t.accessor<double, 3>();

		auto options = x_t.options();
		auto x_tp1 = torch::zeros({ B, 24 * NUM_ACT_SET + 15 }, options);
		auto eta = torch::zeros({ B, NUM_DYN_RESIDUAL }, options);
		auto solver_exit = torch::zeros({ B }, torch::TensorOptions().dtype(torch::kInt32).device(torch::kCPU));
		auto residual_norm = torch::zeros({ B }, options);
		auto out_u0 = torch::zeros({ B, 3 }, options);
		auto out_tau = torch::zeros({ B, NUM_ACT_SET, 3 }, options);

		auto x_tp1_acc = x_tp1.accessor<double, 2>();
		auto eta_acc = eta.accessor<double, 2>();
		auto solver_exit_acc = solver_exit.accessor<int32_t, 1>();
		auto residual_acc = residual_norm.accessor<double, 1>();
		auto out_u0_acc = out_u0.accessor<double, 2>();
		auto out_tau_acc = out_tau.accessor<double, 3>();

		for (int64_t b = 0; b < B; b++) {
			double v_L_pre[NUM_ACT_SET][3];
			double w_L_pre[NUM_ACT_SET][3];
			double mL_initialguess[NUM_ACT_SET][3];
			double nL_initialguess[NUM_ACT_SET][3];
			double p_pre[NUM_ACT_SET][3];
			double R_pre[NUM_ACT_SET][9];
			double xf_pre[NUM_STATES];

			int64_t offset = 0;
			for (int act = 0; act < NUM_ACT_SET; act++) {
				for (int i = 0; i < 3; i++) v_L_pre[act][i] = x_acc[b][offset + act * 3 + i];
			}
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) {
				for (int i = 0; i < 3; i++) w_L_pre[act][i] = x_acc[b][offset + act * 3 + i];
			}
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) {
				for (int i = 0; i < 3; i++) mL_initialguess[act][i] = x_acc[b][offset + act * 3 + i];
			}
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) {
				for (int i = 0; i < 3; i++) nL_initialguess[act][i] = x_acc[b][offset + act * 3 + i];
			}
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) {
				for (int i = 0; i < 3; i++) p_pre[act][i] = x_acc[b][offset + act * 3 + i];
			}
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) {
				for (int i = 0; i < 9; i++) R_pre[act][i] = x_acc[b][offset + act * 9 + i];
			}
			offset += 9 * NUM_ACT_SET;
			for (int i = 0; i < NUM_STATES; i++) xf_pre[i] = x_acc[b][offset + i];

			bool rotations_ok = true;
			for (int act = 0; act < NUM_ACT_SET; act++) {
				if (!CheckRotationMatrix(R_pre[act], cfg->rot_tol)) {
					rotations_ok = false;
					break;
				}
			}
			if (rotations_ok && !CheckRotationMatrix(&xf_pre[3], cfg->rot_tol)) {
				rotations_ok = false;
			}

			if (!rotations_ok) {
				for (int i = 0; i < x_tp1.size(1); i++) x_tp1_acc[b][i] = std::numeric_limits<double>::quiet_NaN();
				for (int i = 0; i < NUM_DYN_RESIDUAL; i++) eta_acc[b][i] = std::numeric_limits<double>::quiet_NaN();
				for (int i = 0; i < 3; i++) out_u0_acc[b][i] = std::numeric_limits<double>::quiet_NaN();
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) out_tau_acc[b][act][i] = std::numeric_limits<double>::quiet_NaN();
				solver_exit_acc[b] = -1;
				residual_acc[b] = std::numeric_limits<double>::infinity();
				continue;
			}

			double act_currents[NUM_ACT_SET][3];
			for (int act = 0; act < NUM_ACT_SET; act++) {
				for (int i = 0; i < 3; i++) act_currents[act][i] = u_acc[b][act][i];
			}

			double Li_val = Li_is_scalar ? Li_scalar : Li_ptr[b];
			CRMShootingMethodParams params = CRMDYNConstructShootingMethodParamSet(
				cfg->cath_params,
				cfg->cath_config,
				Li_val,
				act_currents,
				cfg->contact_mode,
				cfg->tip_constraint_point,
				cfg->tip_force,
				cfg->integration_step_size,
				cfg->act_inertia,
				v_L_pre,
				w_L_pre,
				p_pre,
				R_pre,
				cfg->damping,
				cfg->delta_t);

			double ftip_initialguess[3] = { 0.0, 0.0, 0.0 };
			double out_mL[NUM_ACT_SET][3];
			double out_nL[NUM_ACT_SET][3];
			double out_tau_arr[NUM_ACT_SET][3];
			double out_ftip[3];
			double out_u0_arr[3];
			int localmin = 0;

			DynamicsBVP(params, xf_pre, mL_initialguess, nL_initialguess, ftip_initialguess,
				out_u0_arr, out_mL, out_nL, out_tau_arr, out_ftip, localmin);

			double out_xf[NUM_STATES];
			double out_coil_state[NUM_ACT_SET][NUM_COIL_STATES];
			std::vector<double> out_p_atLocMarkers(params.no_locmarkers * 3);
			DYNSolverIVP(params, out_u0_arr, out_mL, out_nL, out_tau_arr, out_ftip,
				true, out_xf, out_coil_state, reinterpret_cast<double(*)[3]>(out_p_atLocMarkers.data()));

			double eta_scaled[NUM_DYN_RESIDUAL];
			for (int act = 0; act < NUM_ACT_SET; act++) {
				for (int i = 0; i < 3; i++) {
					eta_scaled[act * 6 + i] = out_mL[act][i] / IVALUE_SCALE_M;
					eta_scaled[act * 6 + 3 + i] = out_nL[act][i] / IVALUE_SCALE_N;
				}
			}

			double residual_vec[NUM_DYN_RESIDUAL];
			double u0_tmp[3];
			double tau_tmp[NUM_ACT_SET][3];
			EvaluateDynResidual(params, xf_pre, mL_initialguess, nL_initialguess,
				eta_scaled, residual_vec, u0_tmp, tau_tmp);
			double res_norm_sq = 0.0;
			for (int i = 0; i < NUM_DYN_RESIDUAL; i++) res_norm_sq += residual_vec[i] * residual_vec[i];
			double res_norm = std::sqrt(res_norm_sq);

			int64_t out_offset = 0;
			for (int act = 0; act < NUM_ACT_SET; act++) {
				for (int i = 0; i < 3; i++) x_tp1_acc[b][out_offset + act * 3 + i] = out_coil_state[act][i];
			}
			out_offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) {
				for (int i = 0; i < 3; i++) x_tp1_acc[b][out_offset + act * 3 + i] = out_coil_state[act][i + 3];
			}
			out_offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) {
				for (int i = 0; i < 3; i++) x_tp1_acc[b][out_offset + act * 3 + i] = out_mL[act][i];
			}
			out_offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) {
				for (int i = 0; i < 3; i++) x_tp1_acc[b][out_offset + act * 3 + i] = out_nL[act][i];
			}
			out_offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) {
				for (int i = 0; i < 3; i++) x_tp1_acc[b][out_offset + act * 3 + i] = out_coil_state[act][i + 6];
			}
			out_offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) {
				for (int i = 0; i < 9; i++) x_tp1_acc[b][out_offset + act * 9 + i] = out_coil_state[act][i + 9];
			}
			out_offset += 9 * NUM_ACT_SET;
			for (int i = 0; i < NUM_STATES; i++) x_tp1_acc[b][out_offset + i] = out_xf[i];

			for (int i = 0; i < NUM_DYN_RESIDUAL; i++) eta_acc[b][i] = eta_scaled[i];
			for (int i = 0; i < 3; i++) out_u0_acc[b][i] = out_u0_arr[i];
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) out_tau_acc[b][act][i] = out_tau_arr[act][i];
			solver_exit_acc[b] = static_cast<int32_t>(localmin);
			residual_acc[b] = res_norm;
		}

		return { x_tp1, eta, solver_exit, residual_norm, out_u0, out_tau };
	}

	std::vector<torch::Tensor> crm_dyn_residual(torch::Tensor eta, torch::Tensor x_t, torch::Tensor u_t, torch::Tensor Li, std::shared_ptr<DynamicsConfig> cfg) {
		TORCH_CHECK(cfg != nullptr, "cfg_dyn is null.");
		CheckStepInputs(x_t, u_t, Li);
		TORCH_CHECK(eta.device().is_cpu(), "eta must be a CPU tensor.");
		TORCH_CHECK(eta.dtype() == torch::kFloat64, "eta must be float64.");
		TORCH_CHECK(eta.dim() == 2 && eta.size(1) == NUM_DYN_RESIDUAL, "eta must have shape [B, NUM_DYN_RESIDUAL].");

		x_t = x_t.contiguous();
		u_t = u_t.contiguous();
		Li = Li.contiguous();
		eta = eta.contiguous();

		const auto B = x_t.size(0);
		const bool Li_is_scalar = (Li.numel() == 1);
		const double Li_scalar = Li_is_scalar ? Li.item<double>() : 0.0;
		const double* Li_ptr = Li_is_scalar ? nullptr : Li.data_ptr<double>();

		auto x_acc = x_t.accessor<double, 2>();
		auto u_acc = u_t.accessor<double, 3>();
		auto eta_acc = eta.accessor<double, 2>();

		auto options = x_t.options();
		auto residual = torch::zeros({ B, NUM_DYN_RESIDUAL }, options);
		auto out_u0 = torch::zeros({ B, 3 }, options);
		auto out_tau = torch::zeros({ B, NUM_ACT_SET, 3 }, options);

		auto residual_acc = residual.accessor<double, 2>();
		auto out_u0_acc = out_u0.accessor<double, 2>();
		auto out_tau_acc = out_tau.accessor<double, 3>();

		for (int64_t b = 0; b < B; b++) {
			double v_L_pre[NUM_ACT_SET][3];
			double w_L_pre[NUM_ACT_SET][3];
			double mL_initialguess[NUM_ACT_SET][3];
			double nL_initialguess[NUM_ACT_SET][3];
			double p_pre[NUM_ACT_SET][3];
			double R_pre[NUM_ACT_SET][9];
			double xf_pre[NUM_STATES];

			int64_t offset = 0;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) v_L_pre[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) w_L_pre[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) mL_initialguess[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) nL_initialguess[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) p_pre[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 9; i++) R_pre[act][i] = x_acc[b][offset + act * 9 + i];
			offset += 9 * NUM_ACT_SET;
			for (int i = 0; i < NUM_STATES; i++) xf_pre[i] = x_acc[b][offset + i];

			bool rotations_ok = true;
			for (int act = 0; act < NUM_ACT_SET; act++) {
				if (!CheckRotationMatrix(R_pre[act], cfg->rot_tol)) {
					rotations_ok = false;
					break;
				}
			}
			if (rotations_ok && !CheckRotationMatrix(&xf_pre[3], cfg->rot_tol)) {
				rotations_ok = false;
			}
			if (!rotations_ok) {
				for (int i = 0; i < NUM_DYN_RESIDUAL; i++) residual_acc[b][i] = std::numeric_limits<double>::quiet_NaN();
				for (int i = 0; i < 3; i++) out_u0_acc[b][i] = std::numeric_limits<double>::quiet_NaN();
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) out_tau_acc[b][act][i] = std::numeric_limits<double>::quiet_NaN();
				continue;
			}

			double act_currents[NUM_ACT_SET][3];
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) act_currents[act][i] = u_acc[b][act][i];

			double Li_val = Li_is_scalar ? Li_scalar : Li_ptr[b];
			CRMShootingMethodParams params = CRMDYNConstructShootingMethodParamSet(
				cfg->cath_params,
				cfg->cath_config,
				Li_val,
				act_currents,
				cfg->contact_mode,
				cfg->tip_constraint_point,
				cfg->tip_force,
				cfg->integration_step_size,
				cfg->act_inertia,
				v_L_pre,
				w_L_pre,
				p_pre,
				R_pre,
				cfg->damping,
				cfg->delta_t);

			double eta_scaled[NUM_DYN_RESIDUAL];
			for (int i = 0; i < NUM_DYN_RESIDUAL; i++) eta_scaled[i] = eta_acc[b][i];

			double residual_vec[NUM_DYN_RESIDUAL];
			double u0_tmp[3];
			double tau_tmp[NUM_ACT_SET][3];
			EvaluateDynResidual(params, xf_pre, mL_initialguess, nL_initialguess,
				eta_scaled, residual_vec, u0_tmp, tau_tmp);

			for (int i = 0; i < NUM_DYN_RESIDUAL; i++) residual_acc[b][i] = residual_vec[i];
			for (int i = 0; i < 3; i++) out_u0_acc[b][i] = u0_tmp[i];
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) out_tau_acc[b][act][i] = tau_tmp[act][i];
		}

		return { residual, out_u0, out_tau };
	}

	std::vector<torch::Tensor> crm_dyn_jacobians_fd(torch::Tensor eta, torch::Tensor x_t, torch::Tensor u_t, torch::Tensor Li, std::shared_ptr<DynamicsConfig> cfg) {
		TORCH_CHECK(FDReferenceEnabled(), "FD reference path disabled; set CRM_DIFFSIMS_ENABLE_FD_REFERENCE=1.");
		TORCH_CHECK(cfg != nullptr, "cfg_dyn is null.");
		CheckStepInputs(x_t, u_t, Li);
		TORCH_CHECK(eta.device().is_cpu(), "eta must be a CPU tensor.");
		TORCH_CHECK(eta.dtype() == torch::kFloat64, "eta must be float64.");
		TORCH_CHECK(eta.dim() == 2 && eta.size(1) == NUM_DYN_RESIDUAL, "eta must have shape [B, NUM_DYN_RESIDUAL].");

		x_t = x_t.contiguous();
		u_t = u_t.contiguous();
		Li = Li.contiguous();
		eta = eta.contiguous();

		const auto B = x_t.size(0);
		const bool Li_is_scalar = (Li.numel() == 1);
		const double Li_scalar = Li_is_scalar ? Li.item<double>() : 0.0;
		const double* Li_ptr = Li_is_scalar ? nullptr : Li.data_ptr<double>();

		auto x_acc = x_t.accessor<double, 2>();
		auto u_acc = u_t.accessor<double, 3>();
		auto eta_acc = eta.accessor<double, 2>();

		auto options = x_t.options();
		auto A_step = torch::zeros({ B, NUM_DYN_RESIDUAL, NUM_DYN_RESIDUAL }, options);
		auto B_step = torch::zeros({ B, NUM_DYN_RESIDUAL, 3 * NUM_ACT_SET }, options);
		auto C_step = torch::zeros({ B, NUM_DYN_RESIDUAL, 24 * NUM_ACT_SET + 15 }, options);

		auto A_acc = A_step.accessor<double, 3>();
		auto B_acc = B_step.accessor<double, 3>();
		auto C_acc = C_step.accessor<double, 3>();

		const double fd_eps = 1e-6;
		const int n = NUM_DYN_RESIDUAL;
		const int x_dim = 24 * NUM_ACT_SET + 15;
		const int u_dim = 3 * NUM_ACT_SET;

		for (int64_t b = 0; b < B; b++) {
			double v_L_pre[NUM_ACT_SET][3];
			double w_L_pre[NUM_ACT_SET][3];
			double mL_initialguess[NUM_ACT_SET][3];
			double nL_initialguess[NUM_ACT_SET][3];
			double p_pre[NUM_ACT_SET][3];
			double R_pre[NUM_ACT_SET][9];
			double xf_pre[NUM_STATES];

			int64_t offset = 0;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) v_L_pre[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) w_L_pre[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) mL_initialguess[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) nL_initialguess[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) p_pre[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 9; i++) R_pre[act][i] = x_acc[b][offset + act * 9 + i];
			offset += 9 * NUM_ACT_SET;
			for (int i = 0; i < NUM_STATES; i++) xf_pre[i] = x_acc[b][offset + i];

			bool rotations_ok = true;
			for (int act = 0; act < NUM_ACT_SET; act++) {
				if (!CheckRotationMatrix(R_pre[act], cfg->rot_tol)) {
					rotations_ok = false;
					break;
				}
			}
			if (rotations_ok && !CheckRotationMatrix(&xf_pre[3], cfg->rot_tol)) {
				rotations_ok = false;
			}
			if (!rotations_ok) {
				for (int i = 0; i < NUM_DYN_RESIDUAL; i++) {
					for (int j = 0; j < NUM_DYN_RESIDUAL; j++) A_acc[b][i][j] = std::numeric_limits<double>::quiet_NaN();
					for (int j = 0; j < u_dim; j++) B_acc[b][i][j] = std::numeric_limits<double>::quiet_NaN();
					for (int j = 0; j < x_dim; j++) C_acc[b][i][j] = std::numeric_limits<double>::quiet_NaN();
				}
				continue;
			}

			double act_currents[NUM_ACT_SET][3];
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) act_currents[act][i] = u_acc[b][act][i];

			double Li_val = Li_is_scalar ? Li_scalar : Li_ptr[b];
			CRMShootingMethodParams params = CRMDYNConstructShootingMethodParamSet(
				cfg->cath_params,
				cfg->cath_config,
				Li_val,
				act_currents,
				cfg->contact_mode,
				cfg->tip_constraint_point,
				cfg->tip_force,
				cfg->integration_step_size,
				cfg->act_inertia,
				v_L_pre,
				w_L_pre,
				p_pre,
				R_pre,
				cfg->damping,
				cfg->delta_t);

			double eta_scaled[NUM_DYN_RESIDUAL];
			for (int i = 0; i < NUM_DYN_RESIDUAL; i++) eta_scaled[i] = eta_acc[b][i];

			for (int j = 0; j < n; j++) {
				double eta_plus[NUM_DYN_RESIDUAL];
				double eta_minus[NUM_DYN_RESIDUAL];
				for (int i = 0; i < NUM_DYN_RESIDUAL; i++) {
					eta_plus[i] = eta_scaled[i];
					eta_minus[i] = eta_scaled[i];
				}
				eta_plus[j] += fd_eps;
				eta_minus[j] -= fd_eps;

				double residual_plus[NUM_DYN_RESIDUAL];
				double residual_minus[NUM_DYN_RESIDUAL];
				double u0_plus[3];
				double u0_minus[3];
				double tau_plus[NUM_ACT_SET][3];
				double tau_minus[NUM_ACT_SET][3];
				EvaluateDynResidual(params, xf_pre, mL_initialguess, nL_initialguess, eta_plus, residual_plus, u0_plus, tau_plus);
				EvaluateDynResidual(params, xf_pre, mL_initialguess, nL_initialguess, eta_minus, residual_minus, u0_minus, tau_minus);

				const double inv = 0.5 / fd_eps;
				for (int i = 0; i < n; i++) {
					A_acc[b][i][j] = (residual_plus[i] - residual_minus[i]) * inv;
				}
			}

			for (int j = 0; j < u_dim; j++) {
				double u_plus[NUM_ACT_SET][3];
				double u_minus[NUM_ACT_SET][3];
				for (int act = 0; act < NUM_ACT_SET; act++) {
					for (int i = 0; i < 3; i++) {
						u_plus[act][i] = act_currents[act][i];
						u_minus[act][i] = act_currents[act][i];
					}
				}
				const int act = j / 3;
				const int idx = j % 3;
				u_plus[act][idx] += fd_eps;
				u_minus[act][idx] -= fd_eps;

				CRMShootingMethodParams params_plus = CRMDYNConstructShootingMethodParamSet(
					cfg->cath_params,
					cfg->cath_config,
					Li_val,
					u_plus,
					cfg->contact_mode,
					cfg->tip_constraint_point,
					cfg->tip_force,
					cfg->integration_step_size,
					cfg->act_inertia,
					v_L_pre,
					w_L_pre,
					p_pre,
					R_pre,
					cfg->damping,
					cfg->delta_t);
				CRMShootingMethodParams params_minus = CRMDYNConstructShootingMethodParamSet(
					cfg->cath_params,
					cfg->cath_config,
					Li_val,
					u_minus,
					cfg->contact_mode,
					cfg->tip_constraint_point,
					cfg->tip_force,
					cfg->integration_step_size,
					cfg->act_inertia,
					v_L_pre,
					w_L_pre,
					p_pre,
					R_pre,
					cfg->damping,
					cfg->delta_t);

				double residual_plus[NUM_DYN_RESIDUAL];
				double residual_minus[NUM_DYN_RESIDUAL];
				double u0_plus[3];
				double u0_minus[3];
				double tau_plus[NUM_ACT_SET][3];
				double tau_minus[NUM_ACT_SET][3];
				EvaluateDynResidual(params_plus, xf_pre, mL_initialguess, nL_initialguess, eta_scaled, residual_plus, u0_plus, tau_plus);
				EvaluateDynResidual(params_minus, xf_pre, mL_initialguess, nL_initialguess, eta_scaled, residual_minus, u0_minus, tau_minus);

				const double inv = 0.5 / fd_eps;
				for (int i = 0; i < n; i++) {
					B_acc[b][i][j] = (residual_plus[i] - residual_minus[i]) * inv;
				}
			}

			for (int j = 0; j < x_dim; j++) {
				double x_plus[24 * NUM_ACT_SET + 15];
				double x_minus[24 * NUM_ACT_SET + 15];
				for (int i = 0; i < x_dim; i++) {
					x_plus[i] = x_acc[b][i];
					x_minus[i] = x_acc[b][i];
				}
				x_plus[j] += fd_eps;
				x_minus[j] -= fd_eps;

				double v_L_pre_p[NUM_ACT_SET][3];
				double w_L_pre_p[NUM_ACT_SET][3];
				double mL_guess_p[NUM_ACT_SET][3];
				double nL_guess_p[NUM_ACT_SET][3];
				double p_pre_p[NUM_ACT_SET][3];
				double R_pre_p[NUM_ACT_SET][9];
				double xf_pre_p[NUM_STATES];

				int64_t offset_p = 0;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) v_L_pre_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) w_L_pre_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) mL_guess_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) nL_guess_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) p_pre_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 9; i++) R_pre_p[act][i] = x_plus[offset_p + act * 9 + i];
				offset_p += 9 * NUM_ACT_SET;
				for (int i = 0; i < NUM_STATES; i++) xf_pre_p[i] = x_plus[offset_p + i];

				double v_L_pre_m[NUM_ACT_SET][3];
				double w_L_pre_m[NUM_ACT_SET][3];
				double mL_guess_m[NUM_ACT_SET][3];
				double nL_guess_m[NUM_ACT_SET][3];
				double p_pre_m[NUM_ACT_SET][3];
				double R_pre_m[NUM_ACT_SET][9];
				double xf_pre_m[NUM_STATES];

				int64_t offset_m = 0;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) v_L_pre_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) w_L_pre_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) mL_guess_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) nL_guess_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) p_pre_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 9; i++) R_pre_m[act][i] = x_minus[offset_m + act * 9 + i];
				offset_m += 9 * NUM_ACT_SET;
				for (int i = 0; i < NUM_STATES; i++) xf_pre_m[i] = x_minus[offset_m + i];

				CRMShootingMethodParams params_plus = CRMDYNConstructShootingMethodParamSet(
					cfg->cath_params,
					cfg->cath_config,
					Li_val,
					act_currents,
					cfg->contact_mode,
					cfg->tip_constraint_point,
					cfg->tip_force,
					cfg->integration_step_size,
					cfg->act_inertia,
					v_L_pre_p,
					w_L_pre_p,
					p_pre_p,
					R_pre_p,
					cfg->damping,
					cfg->delta_t);
				CRMShootingMethodParams params_minus = CRMDYNConstructShootingMethodParamSet(
					cfg->cath_params,
					cfg->cath_config,
					Li_val,
					act_currents,
					cfg->contact_mode,
					cfg->tip_constraint_point,
					cfg->tip_force,
					cfg->integration_step_size,
					cfg->act_inertia,
					v_L_pre_m,
					w_L_pre_m,
					p_pre_m,
					R_pre_m,
					cfg->damping,
					cfg->delta_t);

				double residual_plus[NUM_DYN_RESIDUAL];
				double residual_minus[NUM_DYN_RESIDUAL];
				double u0_plus[3];
				double u0_minus[3];
				double tau_plus[NUM_ACT_SET][3];
				double tau_minus[NUM_ACT_SET][3];
				EvaluateDynResidual(params_plus, xf_pre_p, mL_guess_p, nL_guess_p, eta_scaled, residual_plus, u0_plus, tau_plus);
				EvaluateDynResidual(params_minus, xf_pre_m, mL_guess_m, nL_guess_m, eta_scaled, residual_minus, u0_minus, tau_minus);

				const double inv = 0.5 / fd_eps;
				for (int i = 0; i < n; i++) {
					C_acc[b][i][j] = (residual_plus[i] - residual_minus[i]) * inv;
				}
			}
		}

		return { A_step, B_step, C_step };
	}

	std::vector<torch::Tensor> crm_dyn_jacobians(torch::Tensor eta, torch::Tensor x_t, torch::Tensor u_t, torch::Tensor Li, std::shared_ptr<DynamicsConfig> cfg) {
		TORCH_CHECK(cfg != nullptr, "cfg_dyn is null.");
		CheckStepInputs(x_t, u_t, Li);
		TORCH_CHECK(eta.device().is_cpu(), "eta must be a CPU tensor.");
		TORCH_CHECK(eta.dtype() == torch::kFloat64, "eta must be float64.");
		TORCH_CHECK(eta.dim() == 2 && eta.size(1) == NUM_DYN_RESIDUAL, "eta must have shape [B, NUM_DYN_RESIDUAL].");

		x_t = x_t.contiguous();
		u_t = u_t.contiguous();
		Li = Li.contiguous();
		eta = eta.contiguous();

		const auto B = x_t.size(0);
		const bool Li_is_scalar = (Li.numel() == 1);
		const double Li_scalar = Li_is_scalar ? Li.item<double>() : 0.0;
		const double* Li_ptr = Li_is_scalar ? nullptr : Li.data_ptr<double>();

		auto x_acc = x_t.accessor<double, 2>();
		auto u_acc = u_t.accessor<double, 3>();
		auto eta_acc = eta.accessor<double, 2>();

		auto options = x_t.options();
		auto A_step = torch::zeros({ B, NUM_DYN_RESIDUAL, NUM_DYN_RESIDUAL }, options);
		auto B_step = torch::zeros({ B, NUM_DYN_RESIDUAL, 3 * NUM_ACT_SET }, options);
		auto C_step = torch::zeros({ B, NUM_DYN_RESIDUAL, 24 * NUM_ACT_SET + 15 }, options);
		auto A_acc = A_step.accessor<double, 3>();
		auto B_acc = B_step.accessor<double, 3>();
		auto C_acc = C_step.accessor<double, 3>();

		const double fd_eps = 1e-6;
		const int n = NUM_DYN_RESIDUAL;
		const int u_dim = 3 * NUM_ACT_SET;
		const int x_dim = 24 * NUM_ACT_SET + 15;

		for (int64_t b = 0; b < B; b++) {
			double v_L_pre[NUM_ACT_SET][3];
			double w_L_pre[NUM_ACT_SET][3];
			double mL_initialguess[NUM_ACT_SET][3];
			double nL_initialguess[NUM_ACT_SET][3];
			double p_pre[NUM_ACT_SET][3];
			double R_pre[NUM_ACT_SET][9];
			double xf_pre[NUM_STATES];

			int64_t offset = 0;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) v_L_pre[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) w_L_pre[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) mL_initialguess[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) nL_initialguess[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) p_pre[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 9; i++) R_pre[act][i] = x_acc[b][offset + act * 9 + i];
			offset += 9 * NUM_ACT_SET;
			for (int i = 0; i < NUM_STATES; i++) xf_pre[i] = x_acc[b][offset + i];

			bool rotations_ok = true;
			for (int act = 0; act < NUM_ACT_SET; act++) {
				if (!CheckRotationMatrix(R_pre[act], cfg->rot_tol)) {
					rotations_ok = false;
					break;
				}
			}
			if (rotations_ok && !CheckRotationMatrix(&xf_pre[3], cfg->rot_tol)) {
				rotations_ok = false;
			}
			if (!rotations_ok) {
				for (int i = 0; i < NUM_DYN_RESIDUAL; i++) {
					for (int j = 0; j < NUM_DYN_RESIDUAL; j++) A_acc[b][i][j] = std::numeric_limits<double>::quiet_NaN();
				}
				continue;
			}

			double act_currents[NUM_ACT_SET][3];
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) act_currents[act][i] = u_acc[b][act][i];

			double Li_val = Li_is_scalar ? Li_scalar : Li_ptr[b];
			CRMShootingMethodParams params = CRMDYNConstructShootingMethodParamSet(
				cfg->cath_params,
				cfg->cath_config,
				Li_val,
				act_currents,
				cfg->contact_mode,
				cfg->tip_constraint_point,
				cfg->tip_force,
				cfg->integration_step_size,
				cfg->act_inertia,
				v_L_pre,
				w_L_pre,
				p_pre,
				R_pre,
				cfg->damping,
				cfg->delta_t);

			double eta_scaled[NUM_DYN_RESIDUAL];
			for (int i = 0; i < NUM_DYN_RESIDUAL; i++) eta_scaled[i] = eta_acc[b][i];

			for (int j = 0; j < n; j++) {
				double eta_plus[NUM_DYN_RESIDUAL];
				double eta_minus[NUM_DYN_RESIDUAL];
				for (int i = 0; i < NUM_DYN_RESIDUAL; i++) {
					eta_plus[i] = eta_scaled[i];
					eta_minus[i] = eta_scaled[i];
				}

				const double step = (j < 3 * NUM_ACT_SET) ? (fd_eps / IVALUE_SCALE_M) : (fd_eps / IVALUE_SCALE_N);
				eta_plus[j] += step;
				eta_minus[j] -= step;

				double residual_plus[NUM_DYN_RESIDUAL];
				double residual_minus[NUM_DYN_RESIDUAL];
				double u0_plus[3];
				double u0_minus[3];
				double tau_plus[NUM_ACT_SET][3];
				double tau_minus[NUM_ACT_SET][3];
				EvaluateDynResidual(params, xf_pre, mL_initialguess, nL_initialguess, eta_plus, residual_plus, u0_plus, tau_plus);
				EvaluateDynResidual(params, xf_pre, mL_initialguess, nL_initialguess, eta_minus, residual_minus, u0_minus, tau_minus);

				const double inv = 0.5 / step;
				for (int i = 0; i < n; i++) {
					A_acc[b][i][j] = (residual_plus[i] - residual_minus[i]) * inv;
				}
			}

			for (int j = 0; j < u_dim; j++) {
				double u_plus[NUM_ACT_SET][3];
				double u_minus[NUM_ACT_SET][3];
				for (int act = 0; act < NUM_ACT_SET; act++) {
					for (int i = 0; i < 3; i++) {
						u_plus[act][i] = act_currents[act][i];
						u_minus[act][i] = act_currents[act][i];
					}
				}
				const int act = j / 3;
				const int idx = j % 3;
				const double step = fd_eps * std::max(1.0, std::abs(act_currents[act][idx]));
				u_plus[act][idx] += step;
				u_minus[act][idx] -= step;

				CRMShootingMethodParams params_plus = CRMDYNConstructShootingMethodParamSet(
					cfg->cath_params,
					cfg->cath_config,
					Li_val,
					u_plus,
					cfg->contact_mode,
					cfg->tip_constraint_point,
					cfg->tip_force,
					cfg->integration_step_size,
					cfg->act_inertia,
					v_L_pre,
					w_L_pre,
					p_pre,
					R_pre,
					cfg->damping,
					cfg->delta_t);
				CRMShootingMethodParams params_minus = CRMDYNConstructShootingMethodParamSet(
					cfg->cath_params,
					cfg->cath_config,
					Li_val,
					u_minus,
					cfg->contact_mode,
					cfg->tip_constraint_point,
					cfg->tip_force,
					cfg->integration_step_size,
					cfg->act_inertia,
					v_L_pre,
					w_L_pre,
					p_pre,
					R_pre,
					cfg->damping,
					cfg->delta_t);

				double residual_plus[NUM_DYN_RESIDUAL];
				double residual_minus[NUM_DYN_RESIDUAL];
				double u0_plus[3];
				double u0_minus[3];
				double tau_plus[NUM_ACT_SET][3];
				double tau_minus[NUM_ACT_SET][3];
				EvaluateDynResidual(params_plus, xf_pre, mL_initialguess, nL_initialguess, eta_scaled, residual_plus, u0_plus, tau_plus);
				EvaluateDynResidual(params_minus, xf_pre, mL_initialguess, nL_initialguess, eta_scaled, residual_minus, u0_minus, tau_minus);

				const double inv = 0.5 / step;
				for (int i = 0; i < n; i++) {
					B_acc[b][i][j] = (residual_plus[i] - residual_minus[i]) * inv;
				}
			}

			for (int j = 0; j < x_dim; j++) {
				double x_plus[24 * NUM_ACT_SET + 15];
				double x_minus[24 * NUM_ACT_SET + 15];
				for (int i = 0; i < x_dim; i++) {
					x_plus[i] = x_acc[b][i];
					x_minus[i] = x_acc[b][i];
				}

				const double base = x_acc[b][j];
				const double step = fd_eps * (1.0 + 0.01 * std::abs(base));
				x_plus[j] += step;
				x_minus[j] -= step;

				double v_L_pre_p[NUM_ACT_SET][3];
				double w_L_pre_p[NUM_ACT_SET][3];
				double mL_guess_p[NUM_ACT_SET][3];
				double nL_guess_p[NUM_ACT_SET][3];
				double p_pre_p[NUM_ACT_SET][3];
				double R_pre_p[NUM_ACT_SET][9];
				double xf_pre_p[NUM_STATES];

				int64_t offset_p = 0;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) v_L_pre_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) w_L_pre_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) mL_guess_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) nL_guess_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) p_pre_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 9; i++) R_pre_p[act][i] = x_plus[offset_p + act * 9 + i];
				offset_p += 9 * NUM_ACT_SET;
				for (int i = 0; i < NUM_STATES; i++) xf_pre_p[i] = x_plus[offset_p + i];

				double v_L_pre_m[NUM_ACT_SET][3];
				double w_L_pre_m[NUM_ACT_SET][3];
				double mL_guess_m[NUM_ACT_SET][3];
				double nL_guess_m[NUM_ACT_SET][3];
				double p_pre_m[NUM_ACT_SET][3];
				double R_pre_m[NUM_ACT_SET][9];
				double xf_pre_m[NUM_STATES];

				int64_t offset_m = 0;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) v_L_pre_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) w_L_pre_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) mL_guess_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) nL_guess_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) p_pre_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 9; i++) R_pre_m[act][i] = x_minus[offset_m + act * 9 + i];
				offset_m += 9 * NUM_ACT_SET;
				for (int i = 0; i < NUM_STATES; i++) xf_pre_m[i] = x_minus[offset_m + i];

				CRMShootingMethodParams params_plus = CRMDYNConstructShootingMethodParamSet(
					cfg->cath_params,
					cfg->cath_config,
					Li_val,
					act_currents,
					cfg->contact_mode,
					cfg->tip_constraint_point,
					cfg->tip_force,
					cfg->integration_step_size,
					cfg->act_inertia,
					v_L_pre_p,
					w_L_pre_p,
					p_pre_p,
					R_pre_p,
					cfg->damping,
					cfg->delta_t);
				CRMShootingMethodParams params_minus = CRMDYNConstructShootingMethodParamSet(
					cfg->cath_params,
					cfg->cath_config,
					Li_val,
					act_currents,
					cfg->contact_mode,
					cfg->tip_constraint_point,
					cfg->tip_force,
					cfg->integration_step_size,
					cfg->act_inertia,
					v_L_pre_m,
					w_L_pre_m,
					p_pre_m,
					R_pre_m,
					cfg->damping,
					cfg->delta_t);

				double residual_plus[NUM_DYN_RESIDUAL];
				double residual_minus[NUM_DYN_RESIDUAL];
				double u0_plus[3];
				double u0_minus[3];
				double tau_plus[NUM_ACT_SET][3];
				double tau_minus[NUM_ACT_SET][3];
				EvaluateDynResidual(params_plus, xf_pre_p, mL_guess_p, nL_guess_p, eta_scaled, residual_plus, u0_plus, tau_plus);
				EvaluateDynResidual(params_minus, xf_pre_m, mL_guess_m, nL_guess_m, eta_scaled, residual_minus, u0_minus, tau_minus);

				const double inv = 0.5 / step;
				for (int i = 0; i < n; i++) {
					C_acc[b][i][j] = (residual_plus[i] - residual_minus[i]) * inv;
				}
			}
		}

		return { A_step, B_step, C_step };
	}

	std::vector<torch::Tensor> crm_dyn_vjp_fd(torch::Tensor grad_x_tp1, torch::Tensor eta, torch::Tensor x_t, torch::Tensor u_t, torch::Tensor Li, std::shared_ptr<DynamicsConfig> cfg) {
		TORCH_CHECK(FDReferenceEnabled(), "FD reference path disabled; set CRM_DIFFSIMS_ENABLE_FD_REFERENCE=1.");
		TORCH_CHECK(cfg != nullptr, "cfg_dyn is null.");
		CheckStepInputs(x_t, u_t, Li);
		TORCH_CHECK(eta.device().is_cpu(), "eta must be a CPU tensor.");
		TORCH_CHECK(eta.dtype() == torch::kFloat64, "eta must be float64.");
		TORCH_CHECK(eta.dim() == 2 && eta.size(1) == NUM_DYN_RESIDUAL, "eta must have shape [B, NUM_DYN_RESIDUAL].");
		TORCH_CHECK(grad_x_tp1.device().is_cpu(), "grad_x_tp1 must be a CPU tensor.");
		TORCH_CHECK(grad_x_tp1.dtype() == torch::kFloat64, "grad_x_tp1 must be float64.");
		TORCH_CHECK(grad_x_tp1.dim() == 2 && grad_x_tp1.size(1) == 24 * NUM_ACT_SET + 15,
			"grad_x_tp1 must have shape [B, 24*N_act+15].");

		x_t = x_t.contiguous();
		u_t = u_t.contiguous();
		Li = Li.contiguous();
		eta = eta.contiguous();
		grad_x_tp1 = grad_x_tp1.contiguous();

		const auto B = x_t.size(0);
		const bool Li_is_scalar = (Li.numel() == 1);
		const double Li_scalar = Li_is_scalar ? Li.item<double>() : 0.0;
		const double* Li_ptr = Li_is_scalar ? nullptr : Li.data_ptr<double>();

		auto x_acc = x_t.accessor<double, 2>();
		auto u_acc = u_t.accessor<double, 3>();
		auto eta_acc = eta.accessor<double, 2>();
		auto grad_acc = grad_x_tp1.accessor<double, 2>();

		auto options = x_t.options();
		auto vjp_eta = torch::zeros({ B, NUM_DYN_RESIDUAL }, options);
		auto vjp_u = torch::zeros({ B, NUM_ACT_SET, 3 }, options);
		auto vjp_x = torch::zeros({ B, 24 * NUM_ACT_SET + 15 }, options);

		auto vjp_eta_acc = vjp_eta.accessor<double, 2>();
		auto vjp_u_acc = vjp_u.accessor<double, 3>();
		auto vjp_x_acc = vjp_x.accessor<double, 2>();

		const double fd_eps = 1e-6;
		const int x_dim = 24 * NUM_ACT_SET + 15;
		const int u_dim = 3 * NUM_ACT_SET;

		for (int64_t b = 0; b < B; b++) {
			double v_L_pre[NUM_ACT_SET][3];
			double w_L_pre[NUM_ACT_SET][3];
			double mL_initialguess[NUM_ACT_SET][3];
			double nL_initialguess[NUM_ACT_SET][3];
			double p_pre[NUM_ACT_SET][3];
			double R_pre[NUM_ACT_SET][9];
			double xf_pre[NUM_STATES];

			int64_t offset = 0;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) v_L_pre[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) w_L_pre[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) mL_initialguess[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) nL_initialguess[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) p_pre[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 9; i++) R_pre[act][i] = x_acc[b][offset + act * 9 + i];
			offset += 9 * NUM_ACT_SET;
			for (int i = 0; i < NUM_STATES; i++) xf_pre[i] = x_acc[b][offset + i];

			bool rotations_ok = true;
			for (int act = 0; act < NUM_ACT_SET; act++) {
				if (!CheckRotationMatrix(R_pre[act], cfg->rot_tol)) {
					rotations_ok = false;
					break;
				}
			}
			if (rotations_ok && !CheckRotationMatrix(&xf_pre[3], cfg->rot_tol)) {
				rotations_ok = false;
			}
			if (!rotations_ok) {
				for (int i = 0; i < NUM_DYN_RESIDUAL; i++) vjp_eta_acc[b][i] = std::numeric_limits<double>::quiet_NaN();
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) vjp_u_acc[b][act][i] = std::numeric_limits<double>::quiet_NaN();
				for (int i = 0; i < x_dim; i++) vjp_x_acc[b][i] = std::numeric_limits<double>::quiet_NaN();
				continue;
			}

			double act_currents[NUM_ACT_SET][3];
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) act_currents[act][i] = u_acc[b][act][i];

			double Li_val = Li_is_scalar ? Li_scalar : Li_ptr[b];
			CRMShootingMethodParams params = CRMDYNConstructShootingMethodParamSet(
				cfg->cath_params,
				cfg->cath_config,
				Li_val,
				act_currents,
				cfg->contact_mode,
				cfg->tip_constraint_point,
				cfg->tip_force,
				cfg->integration_step_size,
				cfg->act_inertia,
				v_L_pre,
				w_L_pre,
				p_pre,
				R_pre,
				cfg->damping,
				cfg->delta_t);

			double eta_scaled[NUM_DYN_RESIDUAL];
			for (int i = 0; i < NUM_DYN_RESIDUAL; i++) eta_scaled[i] = eta_acc[b][i];

			for (int j = 0; j < NUM_DYN_RESIDUAL; j++) {
				double eta_plus[NUM_DYN_RESIDUAL];
				double eta_minus[NUM_DYN_RESIDUAL];
				for (int i = 0; i < NUM_DYN_RESIDUAL; i++) {
					eta_plus[i] = eta_scaled[i];
					eta_minus[i] = eta_scaled[i];
				}
				eta_plus[j] += fd_eps;
				eta_minus[j] -= fd_eps;

				double x_tp1_plus[24 * NUM_ACT_SET + 15];
				double x_tp1_minus[24 * NUM_ACT_SET + 15];
				EvaluateDynHOutput(params, xf_pre, eta_plus, mL_initialguess, nL_initialguess, x_tp1_plus);
				EvaluateDynHOutput(params, xf_pre, eta_minus, mL_initialguess, nL_initialguess, x_tp1_minus);

				double acc = 0.0;
				const double inv = 0.5 / fd_eps;
				for (int i = 0; i < x_dim; i++) {
					const double diff = (x_tp1_plus[i] - x_tp1_minus[i]) * inv;
					acc += diff * grad_acc[b][i];
				}
				vjp_eta_acc[b][j] = acc;
			}

			for (int j = 0; j < u_dim; j++) {
				double u_plus[NUM_ACT_SET][3];
				double u_minus[NUM_ACT_SET][3];
				for (int act = 0; act < NUM_ACT_SET; act++) {
					for (int i = 0; i < 3; i++) {
						u_plus[act][i] = act_currents[act][i];
						u_minus[act][i] = act_currents[act][i];
					}
				}
				const int act = j / 3;
				const int idx = j % 3;
				u_plus[act][idx] += fd_eps;
				u_minus[act][idx] -= fd_eps;

				CRMShootingMethodParams params_plus = CRMDYNConstructShootingMethodParamSet(
					cfg->cath_params,
					cfg->cath_config,
					Li_val,
					u_plus,
					cfg->contact_mode,
					cfg->tip_constraint_point,
					cfg->tip_force,
					cfg->integration_step_size,
					cfg->act_inertia,
					v_L_pre,
					w_L_pre,
					p_pre,
					R_pre,
					cfg->damping,
					cfg->delta_t);
				CRMShootingMethodParams params_minus = CRMDYNConstructShootingMethodParamSet(
					cfg->cath_params,
					cfg->cath_config,
					Li_val,
					u_minus,
					cfg->contact_mode,
					cfg->tip_constraint_point,
					cfg->tip_force,
					cfg->integration_step_size,
					cfg->act_inertia,
					v_L_pre,
					w_L_pre,
					p_pre,
					R_pre,
					cfg->damping,
					cfg->delta_t);

				double x_tp1_plus[24 * NUM_ACT_SET + 15];
				double x_tp1_minus[24 * NUM_ACT_SET + 15];
				EvaluateDynHOutput(params_plus, xf_pre, eta_scaled, mL_initialguess, nL_initialguess, x_tp1_plus);
				EvaluateDynHOutput(params_minus, xf_pre, eta_scaled, mL_initialguess, nL_initialguess, x_tp1_minus);

				double acc = 0.0;
				const double inv = 0.5 / fd_eps;
				for (int i = 0; i < x_dim; i++) {
					const double diff = (x_tp1_plus[i] - x_tp1_minus[i]) * inv;
					acc += diff * grad_acc[b][i];
				}
				vjp_u_acc[b][act][idx] = acc;
			}

			for (int j = 0; j < x_dim; j++) {
				double x_plus[24 * NUM_ACT_SET + 15];
				double x_minus[24 * NUM_ACT_SET + 15];
				for (int i = 0; i < x_dim; i++) {
					x_plus[i] = x_acc[b][i];
					x_minus[i] = x_acc[b][i];
				}
				x_plus[j] += fd_eps;
				x_minus[j] -= fd_eps;

				double v_L_pre_p[NUM_ACT_SET][3];
				double w_L_pre_p[NUM_ACT_SET][3];
				double mL_guess_p[NUM_ACT_SET][3];
				double nL_guess_p[NUM_ACT_SET][3];
				double p_pre_p[NUM_ACT_SET][3];
				double R_pre_p[NUM_ACT_SET][9];
				double xf_pre_p[NUM_STATES];

				int64_t offset_p = 0;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) v_L_pre_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) w_L_pre_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) mL_guess_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) nL_guess_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) p_pre_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 9; i++) R_pre_p[act][i] = x_plus[offset_p + act * 9 + i];
				offset_p += 9 * NUM_ACT_SET;
				for (int i = 0; i < NUM_STATES; i++) xf_pre_p[i] = x_plus[offset_p + i];

				double v_L_pre_m[NUM_ACT_SET][3];
				double w_L_pre_m[NUM_ACT_SET][3];
				double mL_guess_m[NUM_ACT_SET][3];
				double nL_guess_m[NUM_ACT_SET][3];
				double p_pre_m[NUM_ACT_SET][3];
				double R_pre_m[NUM_ACT_SET][9];
				double xf_pre_m[NUM_STATES];

				int64_t offset_m = 0;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) v_L_pre_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) w_L_pre_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) mL_guess_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) nL_guess_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) p_pre_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 9; i++) R_pre_m[act][i] = x_minus[offset_m + act * 9 + i];
				offset_m += 9 * NUM_ACT_SET;
				for (int i = 0; i < NUM_STATES; i++) xf_pre_m[i] = x_minus[offset_m + i];

				CRMShootingMethodParams params_plus = CRMDYNConstructShootingMethodParamSet(
					cfg->cath_params,
					cfg->cath_config,
					Li_val,
					act_currents,
					cfg->contact_mode,
					cfg->tip_constraint_point,
					cfg->tip_force,
					cfg->integration_step_size,
					cfg->act_inertia,
					v_L_pre_p,
					w_L_pre_p,
					p_pre_p,
					R_pre_p,
					cfg->damping,
					cfg->delta_t);
				CRMShootingMethodParams params_minus = CRMDYNConstructShootingMethodParamSet(
					cfg->cath_params,
					cfg->cath_config,
					Li_val,
					act_currents,
					cfg->contact_mode,
					cfg->tip_constraint_point,
					cfg->tip_force,
					cfg->integration_step_size,
					cfg->act_inertia,
					v_L_pre_m,
					w_L_pre_m,
					p_pre_m,
					R_pre_m,
					cfg->damping,
					cfg->delta_t);

				double x_tp1_plus[24 * NUM_ACT_SET + 15];
				double x_tp1_minus[24 * NUM_ACT_SET + 15];
				EvaluateDynHOutput(params_plus, xf_pre_p, eta_scaled, mL_guess_p, nL_guess_p, x_tp1_plus);
				EvaluateDynHOutput(params_minus, xf_pre_m, eta_scaled, mL_guess_m, nL_guess_m, x_tp1_minus);

				double acc = 0.0;
				const double inv = 0.5 / fd_eps;
				for (int i = 0; i < x_dim; i++) {
					const double diff = (x_tp1_plus[i] - x_tp1_minus[i]) * inv;
					acc += diff * grad_acc[b][i];
				}
				vjp_x_acc[b][j] = acc;
			}
		}

		return { vjp_eta, vjp_u, vjp_x };
	}

	std::vector<torch::Tensor> crm_dyn_vjp(torch::Tensor grad_x_tp1, torch::Tensor eta, torch::Tensor x_t, torch::Tensor u_t, torch::Tensor Li, std::shared_ptr<DynamicsConfig> cfg) {
		TORCH_CHECK(cfg != nullptr, "cfg_dyn is null.");
		CheckStepInputs(x_t, u_t, Li);
		TORCH_CHECK(eta.device().is_cpu(), "eta must be a CPU tensor.");
		TORCH_CHECK(eta.dtype() == torch::kFloat64, "eta must be float64.");
		TORCH_CHECK(eta.dim() == 2 && eta.size(1) == NUM_DYN_RESIDUAL, "eta must have shape [B, NUM_DYN_RESIDUAL].");
		TORCH_CHECK(grad_x_tp1.device().is_cpu(), "grad_x_tp1 must be a CPU tensor.");
		TORCH_CHECK(grad_x_tp1.dtype() == torch::kFloat64, "grad_x_tp1 must be float64.");
		TORCH_CHECK(grad_x_tp1.dim() == 2 && grad_x_tp1.size(1) == 24 * NUM_ACT_SET + 15,
			"grad_x_tp1 must have shape [B, 24*N_act+15].");

		x_t = x_t.contiguous();
		u_t = u_t.contiguous();
		Li = Li.contiguous();
		eta = eta.contiguous();
		grad_x_tp1 = grad_x_tp1.contiguous();

		const auto B = x_t.size(0);
		const bool Li_is_scalar = (Li.numel() == 1);
		const double Li_scalar = Li_is_scalar ? Li.item<double>() : 0.0;
		const double* Li_ptr = Li_is_scalar ? nullptr : Li.data_ptr<double>();

		auto x_acc = x_t.accessor<double, 2>();
		auto u_acc = u_t.accessor<double, 3>();
		auto eta_acc = eta.accessor<double, 2>();
		auto grad_acc = grad_x_tp1.accessor<double, 2>();

		auto options = x_t.options();
		auto vjp_eta = torch::zeros({ B, NUM_DYN_RESIDUAL }, options);
		auto vjp_u = torch::zeros({ B, NUM_ACT_SET, 3 }, options);
		auto vjp_x = torch::zeros({ B, 24 * NUM_ACT_SET + 15 }, options);

		auto vjp_eta_acc = vjp_eta.accessor<double, 2>();
		auto vjp_u_acc = vjp_u.accessor<double, 3>();
		auto vjp_x_acc = vjp_x.accessor<double, 2>();

		const double fd_eps = 1e-6;
		const int x_dim = 24 * NUM_ACT_SET + 15;
		const int u_dim = 3 * NUM_ACT_SET;

		const int offset_v = 0;
		const int offset_w = 3 * NUM_ACT_SET;
		const int offset_m = offset_w + 3 * NUM_ACT_SET;
		const int offset_n = offset_m + 3 * NUM_ACT_SET;

		for (int64_t b = 0; b < B; b++) {
			double v_L_pre[NUM_ACT_SET][3];
			double w_L_pre[NUM_ACT_SET][3];
			double mL_initialguess[NUM_ACT_SET][3];
			double nL_initialguess[NUM_ACT_SET][3];
			double p_pre[NUM_ACT_SET][3];
			double R_pre[NUM_ACT_SET][9];
			double xf_pre[NUM_STATES];

			int64_t offset = 0;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) v_L_pre[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) w_L_pre[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) mL_initialguess[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) nL_initialguess[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) p_pre[act][i] = x_acc[b][offset + act * 3 + i];
			offset += 3 * NUM_ACT_SET;
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 9; i++) R_pre[act][i] = x_acc[b][offset + act * 9 + i];
			offset += 9 * NUM_ACT_SET;
			for (int i = 0; i < NUM_STATES; i++) xf_pre[i] = x_acc[b][offset + i];

			bool rotations_ok = true;
			for (int act = 0; act < NUM_ACT_SET; act++) {
				if (!CheckRotationMatrix(R_pre[act], cfg->rot_tol)) {
					rotations_ok = false;
					break;
				}
			}
			if (rotations_ok && !CheckRotationMatrix(&xf_pre[3], cfg->rot_tol)) {
				rotations_ok = false;
			}
			if (!rotations_ok) {
				for (int i = 0; i < NUM_DYN_RESIDUAL; i++) vjp_eta_acc[b][i] = std::numeric_limits<double>::quiet_NaN();
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) vjp_u_acc[b][act][i] = std::numeric_limits<double>::quiet_NaN();
				for (int i = 0; i < x_dim; i++) vjp_x_acc[b][i] = std::numeric_limits<double>::quiet_NaN();
				continue;
			}

			double act_currents[NUM_ACT_SET][3];
			for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) act_currents[act][i] = u_acc[b][act][i];

			double Li_val = Li_is_scalar ? Li_scalar : Li_ptr[b];
			CRMShootingMethodParams params = CRMDYNConstructShootingMethodParamSet(
				cfg->cath_params,
				cfg->cath_config,
				Li_val,
				act_currents,
				cfg->contact_mode,
				cfg->tip_constraint_point,
				cfg->tip_force,
				cfg->integration_step_size,
				cfg->act_inertia,
				v_L_pre,
				w_L_pre,
				p_pre,
				R_pre,
				cfg->damping,
				cfg->delta_t);

			double eta_scaled[NUM_DYN_RESIDUAL];
			for (int i = 0; i < NUM_DYN_RESIDUAL; i++) eta_scaled[i] = eta_acc[b][i];

			double grad_fd[24 * NUM_ACT_SET + 15];
			for (int i = 0; i < x_dim; i++) grad_fd[i] = grad_acc[b][i];
			for (int act = 0; act < NUM_ACT_SET; act++) {
				for (int i = 0; i < 3; i++) {
					grad_fd[offset_m + act * 3 + i] = 0.0;
					grad_fd[offset_n + act * 3 + i] = 0.0;
				}
			}

			double vjp_eta_direct[NUM_DYN_RESIDUAL] = { 0.0 };
			for (int act = 0; act < NUM_ACT_SET; act++) {
				for (int i = 0; i < 3; i++) {
					vjp_eta_direct[act * 6 + i] += IVALUE_SCALE_M * grad_acc[b][offset_m + act * 3 + i];
					vjp_eta_direct[act * 6 + 3 + i] += IVALUE_SCALE_N * grad_acc[b][offset_n + act * 3 + i];
				}
			}

			for (int j = 0; j < NUM_DYN_RESIDUAL; j++) {
				double eta_plus[NUM_DYN_RESIDUAL];
				double eta_minus[NUM_DYN_RESIDUAL];
				for (int i = 0; i < NUM_DYN_RESIDUAL; i++) {
					eta_plus[i] = eta_scaled[i];
					eta_minus[i] = eta_scaled[i];
				}
				eta_plus[j] += fd_eps;
				eta_minus[j] -= fd_eps;

				double x_tp1_plus[24 * NUM_ACT_SET + 15];
				double x_tp1_minus[24 * NUM_ACT_SET + 15];
				EvaluateDynHOutput(params, xf_pre, eta_plus, mL_initialguess, nL_initialguess, x_tp1_plus);
				EvaluateDynHOutput(params, xf_pre, eta_minus, mL_initialguess, nL_initialguess, x_tp1_minus);

				double acc = 0.0;
				const double inv = 0.5 / fd_eps;
				for (int i = 0; i < x_dim; i++) {
					const double diff = (x_tp1_plus[i] - x_tp1_minus[i]) * inv;
					acc += diff * grad_fd[i];
				}
				vjp_eta_acc[b][j] = acc + vjp_eta_direct[j];
			}

			for (int j = 0; j < u_dim; j++) {
				double u_plus[NUM_ACT_SET][3];
				double u_minus[NUM_ACT_SET][3];
				for (int act = 0; act < NUM_ACT_SET; act++) {
					for (int i = 0; i < 3; i++) {
						u_plus[act][i] = act_currents[act][i];
						u_minus[act][i] = act_currents[act][i];
					}
				}
				const int act = j / 3;
				const int idx = j % 3;
				u_plus[act][idx] += fd_eps;
				u_minus[act][idx] -= fd_eps;

				CRMShootingMethodParams params_plus = CRMDYNConstructShootingMethodParamSet(
					cfg->cath_params,
					cfg->cath_config,
					Li_val,
					u_plus,
					cfg->contact_mode,
					cfg->tip_constraint_point,
					cfg->tip_force,
					cfg->integration_step_size,
					cfg->act_inertia,
					v_L_pre,
					w_L_pre,
					p_pre,
					R_pre,
					cfg->damping,
					cfg->delta_t);
				CRMShootingMethodParams params_minus = CRMDYNConstructShootingMethodParamSet(
					cfg->cath_params,
					cfg->cath_config,
					Li_val,
					u_minus,
					cfg->contact_mode,
					cfg->tip_constraint_point,
					cfg->tip_force,
					cfg->integration_step_size,
					cfg->act_inertia,
					v_L_pre,
					w_L_pre,
					p_pre,
					R_pre,
					cfg->damping,
					cfg->delta_t);

				double x_tp1_plus[24 * NUM_ACT_SET + 15];
				double x_tp1_minus[24 * NUM_ACT_SET + 15];
				EvaluateDynHOutput(params_plus, xf_pre, eta_scaled, mL_initialguess, nL_initialguess, x_tp1_plus);
				EvaluateDynHOutput(params_minus, xf_pre, eta_scaled, mL_initialguess, nL_initialguess, x_tp1_minus);

				double acc = 0.0;
				const double inv = 0.5 / fd_eps;
				for (int i = 0; i < x_dim; i++) {
					const double diff = (x_tp1_plus[i] - x_tp1_minus[i]) * inv;
					acc += diff * grad_acc[b][i];
				}
				vjp_u_acc[b][act][idx] = acc;
			}

			for (int j = 0; j < x_dim; j++) {
				double x_plus[24 * NUM_ACT_SET + 15];
				double x_minus[24 * NUM_ACT_SET + 15];
				for (int i = 0; i < x_dim; i++) {
					x_plus[i] = x_acc[b][i];
					x_minus[i] = x_acc[b][i];
				}
				x_plus[j] += fd_eps;
				x_minus[j] -= fd_eps;

				double v_L_pre_p[NUM_ACT_SET][3];
				double w_L_pre_p[NUM_ACT_SET][3];
				double mL_guess_p[NUM_ACT_SET][3];
				double nL_guess_p[NUM_ACT_SET][3];
				double p_pre_p[NUM_ACT_SET][3];
				double R_pre_p[NUM_ACT_SET][9];
				double xf_pre_p[NUM_STATES];

				int64_t offset_p = 0;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) v_L_pre_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) w_L_pre_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) mL_guess_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) nL_guess_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) p_pre_p[act][i] = x_plus[offset_p + act * 3 + i];
				offset_p += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 9; i++) R_pre_p[act][i] = x_plus[offset_p + act * 9 + i];
				offset_p += 9 * NUM_ACT_SET;
				for (int i = 0; i < NUM_STATES; i++) xf_pre_p[i] = x_plus[offset_p + i];

				double v_L_pre_m[NUM_ACT_SET][3];
				double w_L_pre_m[NUM_ACT_SET][3];
				double mL_guess_m[NUM_ACT_SET][3];
				double nL_guess_m[NUM_ACT_SET][3];
				double p_pre_m[NUM_ACT_SET][3];
				double R_pre_m[NUM_ACT_SET][9];
				double xf_pre_m[NUM_STATES];

				int64_t offset_m = 0;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) v_L_pre_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) w_L_pre_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) mL_guess_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) nL_guess_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 3; i++) p_pre_m[act][i] = x_minus[offset_m + act * 3 + i];
				offset_m += 3 * NUM_ACT_SET;
				for (int act = 0; act < NUM_ACT_SET; act++) for (int i = 0; i < 9; i++) R_pre_m[act][i] = x_minus[offset_m + act * 9 + i];
				offset_m += 9 * NUM_ACT_SET;
				for (int i = 0; i < NUM_STATES; i++) xf_pre_m[i] = x_minus[offset_m + i];

				CRMShootingMethodParams params_plus = CRMDYNConstructShootingMethodParamSet(
					cfg->cath_params,
					cfg->cath_config,
					Li_val,
					act_currents,
					cfg->contact_mode,
					cfg->tip_constraint_point,
					cfg->tip_force,
					cfg->integration_step_size,
					cfg->act_inertia,
					v_L_pre_p,
					w_L_pre_p,
					p_pre_p,
					R_pre_p,
					cfg->damping,
					cfg->delta_t);
				CRMShootingMethodParams params_minus = CRMDYNConstructShootingMethodParamSet(
					cfg->cath_params,
					cfg->cath_config,
					Li_val,
					act_currents,
					cfg->contact_mode,
					cfg->tip_constraint_point,
					cfg->tip_force,
					cfg->integration_step_size,
					cfg->act_inertia,
					v_L_pre_m,
					w_L_pre_m,
					p_pre_m,
					R_pre_m,
					cfg->damping,
					cfg->delta_t);

				double x_tp1_plus[24 * NUM_ACT_SET + 15];
				double x_tp1_minus[24 * NUM_ACT_SET + 15];
				EvaluateDynHOutput(params_plus, xf_pre_p, eta_scaled, mL_guess_p, nL_guess_p, x_tp1_plus);
				EvaluateDynHOutput(params_minus, xf_pre_m, eta_scaled, mL_guess_m, nL_guess_m, x_tp1_minus);

				double acc = 0.0;
				const double inv = 0.5 / fd_eps;
				for (int i = 0; i < x_dim; i++) {
					const double diff = (x_tp1_plus[i] - x_tp1_minus[i]) * inv;
					acc += diff * grad_acc[b][i];
				}
				vjp_x_acc[b][j] = acc;
			}
		}

		return { vjp_eta, vjp_u, vjp_x };
	}

}  // namespace CRMCatheterModel

// v1.3 implicit backward scaffolding (no math yet).
namespace CRMCatheterModel {

torch::Tensor crm_dyn_jacobian_layout_to_row_major(torch::Tensor j_col_major) {
	// TODO(v1.3): implement canonical column-major -> row-major conversion.
	return j_col_major.contiguous();
}

std::vector<torch::Tensor> crm_dyn_vjp_v1_3(
	torch::Tensor grad_x_tp1,
	torch::Tensor eta,
	torch::Tensor x_t,
	torch::Tensor u_t,
	torch::Tensor Li,
	std::shared_ptr<DynamicsConfig> cfg) {
	TORCH_CHECK(cfg != nullptr, "cfg_dyn is null.");
	CheckStepInputs(x_t, u_t, Li);
	TORCH_CHECK(eta.device().is_cpu(), "eta must be a CPU tensor.");
	TORCH_CHECK(eta.dtype() == torch::kFloat64, "eta must be float64.");
	TORCH_CHECK(eta.dim() == 2 && eta.size(1) == NUM_DYN_RESIDUAL, "eta must have shape [B, NUM_DYN_RESIDUAL].");
	TORCH_CHECK(grad_x_tp1.device().is_cpu(), "grad_x_tp1 must be a CPU tensor.");
	TORCH_CHECK(grad_x_tp1.dtype() == torch::kFloat64, "grad_x_tp1 must be float64.");
	TORCH_CHECK(grad_x_tp1.dim() == 2 && grad_x_tp1.size(1) == 24 * NUM_ACT_SET + 15,
		"grad_x_tp1 must have shape [B, 24*N_act+15].");

	auto direct = crm_dyn_vjp(grad_x_tp1, eta, x_t, u_t, Li, cfg);
	auto A_B_C = crm_dyn_jacobians(eta, x_t, u_t, Li, cfg);

	auto vjp_eta_direct = direct[0];
	auto vjp_u_direct = direct[1];
	auto vjp_x_direct = direct[2];

	auto A_step = crm_dyn_jacobian_layout_to_row_major(A_B_C[0]);
	auto B_step = A_B_C[1];
	auto C_step = A_B_C[2];

	const auto B = x_t.size(0);
	const int u_dim = 3 * NUM_ACT_SET;
	const int x_dim = 24 * NUM_ACT_SET + 15;

	auto options = x_t.options();
	auto vjp_u = torch::zeros({ B, NUM_ACT_SET, 3 }, options);
	auto vjp_x = torch::zeros({ B, x_dim }, options);

	for (int64_t b = 0; b < B; b++) {
		auto J_eta = A_step[b];
		auto J_u = B_step[b];
		auto J_x = C_step[b];

		auto rhs = vjp_eta_direct[b].unsqueeze(1);
		auto lambda = at::linalg_solve(J_eta.transpose(0, 1), rhs).squeeze(1);

		auto vjp_x_b = vjp_x_direct[b] - torch::matmul(J_x.transpose(0, 1), lambda);
		auto vjp_u_flat = vjp_u_direct[b].reshape({ u_dim }) - torch::matmul(J_u.transpose(0, 1), lambda);

		vjp_x[b] = vjp_x_b;
		vjp_u[b] = vjp_u_flat.reshape({ NUM_ACT_SET, 3 });
	}

	return { vjp_x, vjp_u };
}

std::vector<torch::Tensor> crm_dyn_backward_v1_3(
	torch::Tensor grad_x_tp1,
	torch::Tensor eta,
	torch::Tensor x_t,
	torch::Tensor u_t,
	torch::Tensor Li,
	std::shared_ptr<DynamicsConfig> cfg) {
	TORCH_CHECK(false, "crm_dyn_backward_v1_3 not implemented (use crm_dyn_vjp_v1_3)");
	return { grad_x_tp1, eta, x_t, u_t, Li };
}

}  // namespace CRMCatheterModel

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {
	py::class_<CRMCatheterModel::DynamicsConfig, std::shared_ptr<CRMCatheterModel::DynamicsConfig>>(m, "DynamicsConfig")
		.def(py::init<const std::string&, const std::string&, double, const std::vector<double>&,
			const std::vector<double>&, const std::vector<double>&, double, double, double>(),
			py::arg("cath_params_path"),
			py::arg("cath_config_path"),
			py::arg("integration_step"),
			py::arg("tip_force"),
			py::arg("act_inertia"),
			py::arg("damping"),
			py::arg("delta_t"),
			py::arg("residual_threshold") = 1e-5,
			py::arg("rot_tol") = 1e-5)
		.def_readonly("residual_threshold", &CRMCatheterModel::DynamicsConfig::residual_threshold)
		.def_readonly("rot_tol", &CRMCatheterModel::DynamicsConfig::rot_tol);

	m.def("crm_step_forward", &CRMCatheterModel::crm_step_forward, "CRM dynamics step forward (v1.1)");
	m.def("crm_dyn_residual", &CRMCatheterModel::crm_dyn_residual, "CRM dynamics residual G(eta, x, u, Li)");
	m.def("crm_dyn_jacobians", &CRMCatheterModel::crm_dyn_jacobians, "CRM dynamics Jacobians A_step, B_step, C_step");
	m.def("crm_dyn_jacobians_fd", &CRMCatheterModel::crm_dyn_jacobians_fd, "CRM dynamics Jacobians A_step, B_step, C_step (pure FD reference)");
	m.def("crm_dyn_vjp", &CRMCatheterModel::crm_dyn_vjp, "CRM dynamics H VJP products (eta, u, x)");
	m.def("crm_dyn_vjp_fd", &CRMCatheterModel::crm_dyn_vjp_fd, "CRM dynamics H VJP products (pure FD reference)");
	m.def("crm_dyn_vjp_v1_3", &CRMCatheterModel::crm_dyn_vjp_v1_3, "CRM dynamics v1.3 implicit VJP stub");
	m.def("crm_dyn_backward_v1_3", &CRMCatheterModel::crm_dyn_backward_v1_3, "CRM dynamics v1.3 backward stub");
	m.def("crm_dyn_jacobian_layout_to_row_major", &CRMCatheterModel::crm_dyn_jacobian_layout_to_row_major, "CRM Jacobian layout conversion helper");
}
