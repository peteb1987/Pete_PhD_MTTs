#!/bin/bash
qsub -m be -l mem_grab=8G -t 1-8 ge_script.sh