#! /bin/bash

foamLog log >/dev/null

gnuplot -presist > /dev/null 2>&1 << EOF
        set logscale y
        set title "Residual vs. Iteration"
        set xlabel "Iteration"
        set ylabel "Residual"
        plot "logs/Ux_0" with lines,\
                "logs/p_0" with lines
EOF
