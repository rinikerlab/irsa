## you might change
## the range lb, and hb in wavenumbers
## the cutoff used (30 cm**-1), the shift factor (mu), and the lorentzian Bandwith w
## the file accordingly (either IR_constant_baseline_500_1800.txt, or IR_constant_baseline_500_1000.txt)
## The peaks set (either peaks_500_1000.txt, or peaks.txt)
lb=500
hb=1000
c=30
s0=0.1
#file=IR_constant_baseline_500_1800.txt
file=IR_constant_baseline_500_1000.txt
for w in {12..12}; do
    for mu in 1.000; do
    python Main.py -lb ${lb} -hb ${hb} -e I_long/gibbs.p --freq I_long/freq.p -exp ${file} -o I_${w}_${lb}_${hb}_${c}_${mu}_constant_baseline -w1 $w -w2 $w -s0 $s0 -s1 0.0192 -c ${c} -m ${mu} -d 1 --peaks peaks_500_1000.txt  
    python Main.py -lb ${lb} -hb ${hb} -e II_long/gibbs.p --freq II_long/freq.p -exp ${file} -o II_${w}_${lb}_${hb}_${c}_${mu}_constant_baseline -w1 $w -w2 $w -s0 $s0 -s1 0.0192 -c ${c} -m ${mu} -d 1 --peaks peaks_500_1000.txt 
done
done
