#!/bin/bash
set -u
rank=${PMIX_RANK:?}
set +e
/usr/bin/time -v -o "memory.rank-$rank.txt" ../../../build/salmon < inputfile > "stdout.rank-$rank" 2> "stderr.rank-$rank"
rc=$?
printf '%s\n' "$rc" > "exit.rank-$rank"
exit "$rc"
