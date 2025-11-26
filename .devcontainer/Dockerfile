FROM mcr.microsoft.com/devcontainers/cpp:1-ubuntu-22.04

# ============================================================
# SYSTEM DEPENDENCIES & BUILD TOOLS
# ============================================================

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    ninja-build \
    git \
    pkg-config \
    gdb \
    valgrind \
    clang \
    clang-format \
    && rm -rf /var/lib/apt/lists/*

# ============================================================
# C++ LIBRARIES (Eigen3, pybind11)
# ============================================================

RUN apt-get update && apt-get install -y --no-install-recommends \
    libeigen3-dev \
    pybind11-dev \
    && rm -rf /var/lib/apt/lists/*

# ============================================================
# PYTHON 3 & ML STACK
# ============================================================

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Create virtual environment in container
ENV VENV_PATH="/opt/venv"
RUN python3 -m venv ${VENV_PATH}
ENV PATH="${VENV_PATH}/bin:$PATH"

# Upgrade pip and install core ML/scientific packages
RUN pip install --upgrade pip setuptools wheel && \
    pip install \
    numpy \
    scipy \
    matplotlib \
    scikit-learn \
    pandas \
    jupyterlab \
    ipython \
    pybind11[global] \
    && rm -rf /tmp/pip-*

# Install PyTorch (CPU version; adjust for CUDA if needed)
RUN pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu

# ============================================================
# OPTIONAL: Medical/Robotics ML Libraries
# ============================================================

RUN pip install \
    pillow \
    opencv-python \
    scikit-image \
    && rm -rf /tmp/pip-*

# ============================================================
# DOCUMENTATION & ADDITIONAL TOOLS
# ============================================================

RUN apt-get update && apt-get install -y --no-install-recommends \
    doxygen \
    graphviz \
    && rm -rf /var/lib/apt/lists/*

# ============================================================
# CLEANUP
# ============================================================

RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
