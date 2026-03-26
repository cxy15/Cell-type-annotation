#!/usr/bin/env bash
# Convenience wrapper — implementation lives in pipeline/run_workflow.sh
exec "$(cd "$(dirname "$0")" && pwd)/pipeline/run_workflow.sh" "$@"
