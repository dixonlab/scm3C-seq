python allc_modifier.py ${1} &
python allc_modifier.py ${2} &
wait
python two_integrated_allc_integrate.py ${1}_CpG_mod.txt ${2}_CpG_mod.txt ${3}
