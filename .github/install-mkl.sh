#!/bin/sh
# Add Intel's signed apt repository and install the packages named by the job.
set -eu

key_dir=$(mktemp -d)
trap 'rm -rf "$key_dir"' EXIT
wget -qO "$key_dir/intel.pub" \
    https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
gpg --batch --dearmor --output "$key_dir/oneapi.gpg" "$key_dir/intel.pub"
sudo install -m 0644 "$key_dir/oneapi.gpg" /usr/share/keyrings/oneapi-archive-keyring.gpg
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" \
    | sudo tee /etc/apt/sources.list.d/oneAPI.list > /dev/null
sudo apt-get update
sudo apt-get install -y --no-install-recommends "$@"
