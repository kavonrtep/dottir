gcc -O2 -Wall -Wno-misleading-indentation \
    -o tests/golden_gen/karlin_ref tests/golden_gen/karlin_ref.c -lm
tests/golden_gen/karlin_ref > tests/golden/karlin/values.tsv
