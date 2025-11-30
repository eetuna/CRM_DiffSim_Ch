#!/bin/bash
set -e

if [ -f /etc/os-release ]; then
    . /etc/os-release
    OS=$ID
else
    echo "Cannot detect OS"
    exit 1
fi

echo "Installing dependencies for $OS..."

case $OS in
    ubuntu|debian)
        sudo apt-get update
        sudo apt-get install -y \
            build-essential cmake ninja-build \
            libeigen3-dev pybind11-dev \
            python3-dev python3-pip python3-venv
        ;;
    fedora|rhel|centos)
        sudo dnf install -y \
            gcc-c++ cmake ninja-build \
            eigen3-devel pybind11-devel \
            python3-devel python3-pip
        ;;
    arch)
        sudo pacman -S --noconfirm \
            base-devel cmake ninja eigen pybind11 \
            python python-pip
        ;;
    *)
        echo "Unsupported OS: $OS"
        exit 1
        ;;
esac

echo "Installing Python dependencies..."
pip3 install --user -r requirements.txt

echo "Dependencies installed successfully!"
