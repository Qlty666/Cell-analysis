#!/usr/bin/env bash
set -e
cd "$(dirname "$0")"
python3 launchers/install_environment.py install full --with-ml
echo "Setup complete. Run check_new_computer.sh to verify."
