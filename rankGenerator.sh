#!/bin/bash

DGE=$1
RNK=$(echo "$DGE" | sed 's/.txt/.rnk/')

sed 1d "$DGE" \
| sort -k12g \
| awk '!arr[$1]++' \
| awk '{ printf "%s\t%.6f\n", $1, $11 }' \
| sort -k2gr > "$RNK"

