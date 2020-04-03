#!/bin/bash
# Anoushka Joglekar
# April 2020

## Check if bedtools and samtools are loaded
if which bedtools >/dev/null; then
    echo bedtools found in path
else
    echo bedtools does not exist in path
    exit 127
fi

if which samtools >/dev/null; then
    echo samtools found in path
else
    echo samtools does not exist in path
    exit 127
fi
