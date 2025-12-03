# MRI-Actuated Catheter ML/RL Framework

Physics-informed machine learning and reinforcement learning for MRI-guided magnetically actuated robotic catheters.

## Overview

This framework integrates:
- **Analytical Model**: Your C++ Cosserat rod solver (CRM + CRMDYN)
- **ML Residual Learning**: Neural networks learn corrections to analytical predictions
- **RL Control**: TD3 agent for autonomous catheter navigation

Based on:
- Kinematics: "Real-Time Forward Kinematics and Jacobian for Control of an MRI-Guided Magnetically Actuated Robotic Catheter" (TMECH 2024)
- Dynamics: "Free-Space Dynamic Modeling of an MRI-Actuated Robotic Catheter" (TMECH 2025)
- Cosserat Theory: https://www.cosseratrods.org/

## Architecture

┌─────────────────────────────────────────────┐
│ C++ Cosserat Rod Solver (CRM + CRMDYN) │
│ - Hybrid CRDM (6.2ms/step) │
│ - BVP/IVP solvers │
│ - Jacobian computation │
└──────────────┬──────────────────────────────┘
│
├──────────────┬────────────────┐
▼ ▼ ▼
┌─────────────┐ ┌─────────────┐ ┌────────────┐
│ Kinematics │ │ Dynamics │ │ RL │
│ (Hybrid) │ │ (Hybrid) │ │ (TD3 Agent)│
└─────────────┘ └─────────────┘ └────────────┘
│ │ │
▼ ▼ ▼
┌─────────────┐ ┌─────────────┐ ┌────────────┐
│ML Correction│ │ML Correction│ │Policy Net │
│ (residual) │ │ (residual) │ │(currents) │
└─────────────┘ └─────────────┘ └────────────┘

text

## Installation

### Prerequisites

- Python 3.10+
- C++17 compiler
- Your existing CRM C++ codebase
- CMake 3.15+

### Build

Install Python dependencies

pip install -r requirements.txt
Build C++ bindings

cd crm_ml_rl/bindings
mkdir build && cd build
cmake ..
make -j8
cd ../../..
Install package

pip install -e .

text

## Usage

### 1. Generate Training Data

python examples/01_generate_data.py
--n_samples 50000
--output_dir data/

text

### 2. Train Kinematics Residual Model

python examples/02_train_kinematics.py
--data_path data/kinematics_data.npz
--epochs 100
--batch_size 512

text

### 3. Train Dynamics Model

python examples/03_train_dynamics.py
--data_path data/dynamics_data.npz
--epochs 100

text

### 4. Train RL Agent

python examples/04_train_rl.py
--n_episodes 5000
--dynamics_model models/dynamics.pth

text

### 5. Evaluate

python examples/05_evaluate.py
--agent models/td3_agent.pth
--n_trials 100

text

## Catheter Specifications

From experimental prototype:
- Total length: 146 mm
- Flexible segments: 2 (104mm, 29mm)
- Rigid segment (coil): 1 (13mm)
- Coils per actuator: 3 (1 axial + 2 side)
- Material: Pebax 35D
- Young's Modulus: 31.03 MPa
- Shear Modulus: 8.11 MPa
- MRI Field: 3T (constant)

## Control

Action (6D):

    Coil currents: I ∈ [-2, 2]⁶ Amps
    [I_axial_1, I_side1_1, I_side2_1,
    I_axial_2, I_side1_2, I_side2_2]

text

## Hybrid Architecture

Kinematics prediction

p_predicted = p_analytical(q, I) + residual_NN(q, I)
Dynamics prediction

state_next = dynamics_analytical(state, I) + correction_NN(state, I)
RL uses hybrid dynamics for environment

env = MRICatheterEnv(dynamics=hybrid_dynamics)
agent = TD3Agent(state_dim=18, action_dim=6)

text

## References

1. Itsarachaiyot et al. "Real-Time Forward Kinematics and Jacobian for Control of an MRI-Guided Magnetically Actuated Robotic Catheter." IEEE TMECH, 2024.

2. Hao et al. "Free-Space Dynamic Modeling of an MRI-Actuated Robotic Catheter." IEEE/ASME TMECH, 2025.

3. Cosserat Rods Theory: https://www.cosseratrods.org/

## License

[Your License]

## Directory Structure
crm_ml_rl/
├── README.md                           # File 3
├── requirements.txt                    # File 1
├── setup.py                           # File 2
├── crm_ml_rl/
│   ├── __init__.py                    # File 4
│   ├── bindings/
│   │   ├── __init__.py                # File 5
│   │   ├── CMakeLists.txt             # File 6
│   │   └── crm_bindings.cpp           # File 7 (YOUR C++ code)
│   ├── api/
│   │   ├── __init__.py                # File 8
│   │   ├── kinematics.py              # File 9 (with Jacobian)
│   │   ├── dynamics.py                # File 10
│   │   └── hybrid.py                  # File 11
│   ├── models/
│   │   ├── __init__.py                # File 12
│   │   ├── kinematics_residual.py     # File 13
│   │   └── dynamics_model.py          # File 14
│   ├── rl/
│   │   ├── __init__.py                # File 15
│   │   ├── environment.py             # File 16
│   │   ├── td3_agent.py               # File 17
│   │   └── training.py                # File 18
│   └── utils/
│       ├── __init__.py                # File 19
│       ├── config.py                  # File 20
│       └── visualization.py           # File 21
├── examples/
│   ├── 01_generate_data.py            # File 22
│   ├── 02_train_kinematics.py         # File 23
│   ├── 03_train_dynamics.py           # File 24
│   ├── 04_train_rl.py                 # File 25
│   └── 05_evaluate.py                 # File 26
└── tests/
    ├── __init__.py                    # File 27
    ├── test_kinematics.py             # File 28
    ├── test_dynamics.py               # File 29
    └── test_rl.py                     # File 30

