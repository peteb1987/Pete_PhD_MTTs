#!/bin/bash
qsub -m be -l lr=1 -t 1-11 ge_script.sh