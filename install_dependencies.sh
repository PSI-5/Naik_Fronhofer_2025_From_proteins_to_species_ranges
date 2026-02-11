#!/bin/bash
set -e

sudo apt update
sudo apt install -y \
    build-essential \
    libgsl-dev \
    openmpi-bin \
    libopenmpi-dev \
    python3 \
    python3-pip

pip3 install numpy mpi4py

