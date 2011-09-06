#!/bin/bash
qsub -m be -l mem_grab=8G -t 1-10 ge_script.sh