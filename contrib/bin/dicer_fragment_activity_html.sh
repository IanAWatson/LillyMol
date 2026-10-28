#!/usr/bin/env bash
exec "${PYTHON:-python3}" "$(dirname -- "$0")/dicer_fragment_activity_html.py" "$@"
