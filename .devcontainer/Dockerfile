FROM ubuntu:22.04 AS base

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake ninja-build git pkg-config \
    gdb valgrind clang clang-format clang-tidy \
    wget curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*

FROM base AS cpp-deps

RUN apt-get update && apt-get install -y --no-install-recommends \
    libeigen3-dev pybind11-dev \
    && rm -rf /var/lib/apt/lists/*

FROM cpp-deps AS python-env

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.10 python3.10-dev python3-pip python3-venv \
    && rm -rf /var/lib/apt/lists/*

ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN pip install --no-cache-dir --upgrade pip setuptools wheel

FROM python-env AS development

COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt && rm /tmp/requirements.txt

RUN apt-get update && apt-get install -y --no-install-recommends \
    doxygen graphviz vim nano htop tmux sudo \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

ARG USERNAME=developer
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

USER $USERNAME

RUN echo "source /opt/venv/bin/activate" >> /home/$USERNAME/.bashrc

CMD ["/bin/bash"]
