#!/usr/bin/env bash
set -e
cd "$(dirname "$0")"
python3 launchers/install_environment.py check full
echo "Environment check passed."
