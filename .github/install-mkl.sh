#!/bin/sh
# Install Intel oneAPI MKL from Intel's signed apt repository, plus any extra
# Ubuntu packages given as arguments. Used by the CI workflows.
set -eu

wget -qO- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
    | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" \
    | sudo tee /etc/apt/sources.list.d/oneAPI.list > /dev/null
sudo apt-get update
sudo apt-get install -y --no-install-recommends intel-oneapi-mkl-devel "$@"
