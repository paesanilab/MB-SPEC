#!/bin/sh
set -eu

CDPATH= cd "$(dirname "$0")"
exec autoreconf --force --install "$@"
