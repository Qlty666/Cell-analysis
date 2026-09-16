#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PID_FILE="/tmp/liverbio-web-ui.pid"
LOG_FILE="/tmp/liverbio-web-ui.log"

cd "$ROOT_DIR"

if [ -f "$PID_FILE" ] && kill -0 "$(cat "$PID_FILE")" 2>/dev/null; then
  exit 0
fi

nohup python web/web_ui.py \
  --host 127.0.0.1 \
  --port 8000 \
  --no-browser \
  --keep-alive \
  >"$LOG_FILE" 2>&1 &

echo "$!" >"$PID_FILE"

for _ in $(seq 1 30); do
  if python -c \
    'import urllib.request; urllib.request.urlopen("http://127.0.0.1:8000/guide", timeout=1).read()' \
    >/dev/null 2>&1; then
    exit 0
  fi
  sleep 1
done

cat "$LOG_FILE" >&2 || true
exit 1
