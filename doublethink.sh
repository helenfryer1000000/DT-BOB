#!/bin/bash
# doublethink.sh: run any Doublethink method from one config file.
#
#   ./doublethink.sh <config file>
#
# The config's 'method' line chooses the method:
#   method="mcmc"          original Doublethink MCMC (SLURM array job, code/qsub_wrapper.sh)
#   method="best_of_both"  best model with vs best model without, found by search (bob/run-bob.R)
#   method="top_m"         sum of the m best models with vs without, m = top_m <= 10 (bob/run-bob.R)
#   method="exhaustive"    exact Doublethink by enumerating all 2^nu models (nu <= 20) (bob/run-bob.R)
# A config without a 'method' line runs the MCMC, so the original configs work unchanged.
# See README.md and docs/USAGE.md.

set -e
if [ "$#" -ne 1 ]; then
	echo "Usage: $0 <config file>"
	exit 1
fi
cfg=$1
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
method="mcmc"
. "$cfg"
case "$method" in
	mcmc)
		exec "$repo/code/qsub_wrapper.sh" "$cfg" ;;
	best_of_both|top_m|exhaustive)
		exec Rscript "$repo/bob/run-bob.R" "$cfg" ;;
	*)
		echo "Unknown method '$method': use mcmc, best_of_both, top_m or exhaustive"
		exit 1 ;;
esac
