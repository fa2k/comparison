
awk -v OFS='\t' '/^sampleName/ {print;next;} {$1=$1."-aa"; print;}' downsample_factors.txt > downsample_factors_fix.txt
