FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# Minimal system dependencies
RUN apt-get update && apt-get install --no-install-recommends -y \
    build-essential \
    cmake \
    ninja-build \
    python3.10 \
    python3-pip \
    python3.10-dev \
    libeigen3-dev \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Python setup
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.10 1 && \
    update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1

# Minimal pip packages
RUN pip install --no-cache-dir pip setuptools wheel "pybind11[global]"

WORKDIR /workspace

ENV PYTHONPATH=/workspace:/workspace/crm_ml_rl:$PYTHONPATH
ENV PATH=/workspace/scripts:$PATH

CMD ["/bin/bash"]
