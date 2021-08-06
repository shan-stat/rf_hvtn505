#!/bin/bash

for (( i = 1; i<=500; i++)); do
 sbatch  --wrap="R --no-save --no-restore < test_st.R --args 1 $i RF:tcell+glm:antibody ST" --output=out.txt
done
