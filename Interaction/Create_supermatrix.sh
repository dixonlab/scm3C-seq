#Usage: sh 100kb_Super_matrix_generate.sh input reference_chrNameLength_file             
cooler makebins reference_chrNameLength_file ${2} 10000 > ${2}_10kb_bin.txt
cooler csort --nproc 2 -c1 2 -p1 3 -c2 4 -p2 5 -o ${1}.sorted ${1} ${2}_10kb_bin.txt
cooler cload pairix ${3}  ${1}.sorted ${1}_10kb_contacts.cool
cooler dump -t pixels --header --join ${1}_10kb_contacts.cool > ${1}_10kb_supermatrix.txt
