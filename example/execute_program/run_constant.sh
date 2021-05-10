lb=500
hb=1000
c=20 ## change to different value if desired (cutoff in scoring function)
s0=0.1
file=IR_constant_baseline_500_1000
for w in 12; do #change if desired
    for mu in 1.000; do # change if desired
    python Main.py -lb ${lb} -hb ${hb} -e I_long_correct/gibbs.p --freq I_long_correct/freq.p -exp ${file} -o I_${w}_${lb}_${hb}_${c}_${mu}_constant_baseline -w1 $w -w2 $w -s0 $s0 -s1 0.0192 -c ${c} -m ${mu} -d 1 --peaks peaks_500_1000
    python Main.py -lb ${lb} -hb ${hb} -e II_long_correct/gibbs.p --freq II_long_correct/freq.p -exp ${file} -o II_${w}_${lb}_${hb}_${c}_${mu}_constant_baseline -w1 $w -w2 $w -s0 $s0 -s1 0.0192 -c ${c} -m ${mu} -d 1 --peaks peaks_500_1000
done
done
lb=500 
hb=1800
c=20 ## change to different value if desired (cutoff in scoring function)
s0=0.1
file=IR_constant_baseline_500_1800
for w in 12; do ##change to different value if desired (Lorentzian Bandwidth)
    for mu in 1.000; do ##change to different value if desired (Shifting factor)
    python Main.py -lb ${lb} -hb ${hb} -e I_long_correct/gibbs.p --freq I_long_correct/freq.p -exp ${file} -o I_${w}_${lb}_${hb}_${c}_${mu}_constant_baseline -w1 $w -w2 $w -s0 $s0 -s1 0.0192 -c ${c} -m ${mu} -d 1 --peaks peaks
    python Main.py -lb ${lb} -hb ${hb} -e II_long_correct/gibbs.p --freq II_long_correct/freq.p -exp ${file} -o II_${w}_${lb}_${hb}_${c}_${mu}_constant_baseline -w1 $w -w2 $w -s0 $s0 -s1 0.0192 -c ${c} -m ${mu} -d 1 --peaks peaks
done
done
