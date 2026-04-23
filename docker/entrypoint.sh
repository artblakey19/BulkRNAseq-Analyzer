#!/bin/bash
# Generic UID/HOME shim used by both the pipeline and the Jupyter image.
#
# Reads host UID/GID off the bind-mounted /project (so outputs land with host
# ownership), overridable via HOST_UID / HOST_GID env vars. Points HOME at
# /tmp so conda/mamba/jupyter have a writable cache dir. Then drops
# privileges with gosu and exec's whatever inner command the image supplies.
set -euo pipefail

HOST_UID="${HOST_UID:-$(stat -c %u /project 2>/dev/null || echo 0)}"
HOST_GID="${HOST_GID:-$(stat -c %g /project 2>/dev/null || echo 0)}"

# gosu 1.17 resets HOME to the passwd entry of the target UID, falling back
# to "/" when the UID isn't in /etc/passwd — which is the common case since
# host UIDs almost never match mambauser (57439). Pass HOME=/tmp through
# `env` so it survives the privilege drop; otherwise libmamba/jupyter try to
# write caches under "/" and fail.
exec gosu "${HOST_UID}:${HOST_GID}" env HOME=/tmp PROJECT_ROOT=/project "$@"
