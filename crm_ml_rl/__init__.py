"""
MRI-Actuated Catheter ML/RL Framework

Physics-informed machine learning and reinforcement learning for 
MRI-guided magnetically actuated robotic catheters.

Based on:
- Kinematics: Real-Time Forward Kinematics (TMECH 2024)
- Dynamics: Free-Space Dynamic Modeling (TMECH 2025)
- Cosserat Theory: https://www.cosseratrods.org/
"""

__version__ = "1.0.0"

from crm_ml_rl.api import (
    ForwardKinematicsAPI,
    InverseKinematicsAPI,
    DynamicsAPI,
    HybridKinematics,
    HybridDynamics
)

from crm_ml_rl.models import (
    KinematicsResidualNN,
    DynamicsModelNN
)

from crm_ml_rl.rl import (
    MRICatheterEnv,
    TD3Agent,
    train_td3
)

__all__ = [
    # API
    'ForwardKinematicsAPI',
    'InverseKinematicsAPI',
    'DynamicsAPI',
    'HybridKinematics',
    'HybridDynamics',
    # Models
    'KinematicsResidualNN',
    'DynamicsModelNN',
    # RL
    'MRICatheterEnv',
    'TD3Agent',
    'train_td3',
]
