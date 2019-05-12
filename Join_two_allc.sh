python allc_modifier.py ${1} &
python allc_modifier.py ${2} &
wait
python Join_two_modified_files_to_one.py ${1}_CpG_mod.txt ${2}_CpG_mod.txt ${3}
