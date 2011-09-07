#!/bin/bash
qsub -m be -l lr=1 -t 1-10 ge_script.sh