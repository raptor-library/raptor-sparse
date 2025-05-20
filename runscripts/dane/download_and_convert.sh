# I chose all the matrices that had between 20 and 30 million nonzeroes
filepath="https://suitesparse-collection-website.herokuapp.com/MM"
n_mats=9
mats=("delaunay_n23" "dielFilterV2real" "hugetrace-00020" "bone010" "spal_004" "wikipedia-20070206" "kron_g500-logn19" "ldoor" "boneS10" "Hardesty3")
fns=("DIMACS10/delaunay_n23" "Dziekonski/dielFilterV2real" "DIMACS10/hugetrace-00020" "Oberwolfach/bone010" "Mittelmann/spal_004" "Gleich/wikipedia-20070206" "DIMACS10/kron_g500-logn19" "GHS_psdef/ldoor" "Oberwolfach/boneS10" "Hardesty/Hardesty3")


for i in $(seq 0 $n_mats);
do
    mat=${mats[i]}
    echo "Testing i $i, mat ${mats[i]}"
    if test -f ${mat}.pm; then
        echo "${mat}.pm already exists"
    else
        if test -f ${mat}.mtx; then
			echo "${mat}.mtx already exists"
		else
            if test -f ${mat}; then
				echo "directory ${mat} already exists"
			else
                if test -f ${mat}.tar.gz; then
    		        echo "${mat}.tar.gz already exists"
		        else
        	        wget $filepath/${fns[i]}.tar.gz
                fi
                tar -xzvf ${mat}.tar.gz
            fi
            mv ${mat}/${mat}.mtx .
            rm -rf ${mat}/
            rm ${mat}.tar.gz
        fi
        python3 convert_petsc.py ${mat}.mtx ${mat}.pm
        #rm ${mat}.mtx
    fi 
done
