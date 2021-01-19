#!/bin/bash
function make_new_crystal_restart() {
	job_dir=( $@  ) #step_102/no_d3 step_104/no_d3 step_106/no_d3 step_107/no_d3 step_108/no_d3 ) #step_101/no_d3 step_102/no_d3 step_103/no_d3 step_104/no_d3 step_105/no_d3 step_106/no_d3 step_107/no_d3 step_108/no_d3 step_109/no_d3 step_110/no_d3 ) # step_101  step_102  step_103  step_104  step_105  step_106  step_107  step_108  step_109 ) #  step_110 ) #CoM_11NN_closespecified_spins CoM_11NN_farspecified_spins  CoM_12NNspecified_spins  CoM_13NNspecified_spins ) # Cu1NN_close  Cu1NN_far  Cu2NN  Cu3NN )
	home=$PWD
	for dir in ${job_dir[@]}; do
		cd $home/$dir
		if [ -s "job.out" ]; then # so job.out exists and is non-zero
			if grep -q "OPT END - CONVERGED" job.out; then # so job is finished
				echo "Job $dir has finished"
			else # original job is unfinished
				if [ -d "restart_1" ]; then # so there is at least 1 restart present
					res_num=$( ls -d restart_* | cut -d '_' -f2 |  sort -n | tail -n 1 )
					if grep -q "OPT END - CONVERGED" $home/$dir/restart_${res_num}/job.out; then
						echo "Job $dir/restart_${res_num} has finished"
					else # so the latest restart has not finished
						new_res_num=$((res_num+1))
						cd $home/$dir/restart_${res_num}
						optc_file=$( ls optc*[0-9] | sort | tail -n 1) # final step
						if [ -z $optc_file ]; then # so there are no optc files
			    	  echo "The job $dir has not done any steps --- ERROR"
			    	else # so there is at least 1 step but no restart
							mkdir ../restart_${new_res_num}
							crystal_to_poscar $optc_file
			    	  cp $optc_file.POSCAR.vasp ../restart_${new_res_num}
							sed -e "s/pe mpi.*/pe mpi 48/" -e "s/P Free/P Gold/" -e "s/restart_${res_num}/restart_${new_res_num}/" job_script.pbs > ../restart_${new_res_num}/job_script.pbs
			    	  cp $optc_file ../restart_${new_res_num}/fort.34
	            if [[ -s fort.9 ]]; then
	              cp fort.9 ../restart_${new_res_num}/fort.20
	            else
	              cp fort.20 ../restart_${new_res_num}/fort.20
	            fi
			    	  sed "/FMIXING/,+1 s/^[0-9].*/85/" INPUT > ../restart_${new_res_num}/INPUT # note if resing from Co then set FMIXING as high like 80 else do ~5%
							if ! grep -q GUESSP ../restart_${new_res_num}/INPUT; then
								sed -i "/ENDDFT/ a\ GUESSP" ../restart_${new_res_num}/INPUT
							fi
							echo "Job $dir/restart_${res_num} has ran steps :: run $dir/restart_${new_res_num} to continue"
						fi
					fi
				else # i.e. there is no restart present
					optc_file=$( ls optc*[0-9] | sort | tail -n 1) # final step
					if [ -z $optc_file ]; then # so there are no optc files
						echo "The job $dir has not done any steps --- ERROR"
					else # so there is at least 1 step but no restart
						res_num=1
						mkdir ./restart_${res_num}
						crystal_to_poscar $optc_file
						cp {job_script.pbs,$optc_file.POSCAR.vasp} ./restart_${res_num}
						cp $optc_file ./restart_${res_num}/fort.34
	          if [[ -s fort.9 ]]; then
	            cp fort.9 ../restart_${new_res_num}/fort.20
	          else
	            cp fort.20 ../restart_${new_res_num}/fort.20
	          fi
						title=$( head -n 1 INPUT )
						sed "/OPTGEOM/,$ !d" INPUT > ./restart_${res_num}/INPUT
						sed -i "/OPTGEOM/ i\ $title" ./restart_${res_num}/INPUT 
						sed -i "/OPTGEOM/ i\ EXTERNAL" ./restart_${res_num}/INPUT
						sed -i "/FMIXING/,+1 s/^60/75/" ./restart_${res_num}/INPUT
						sed -i "/ENDDFT/ a\ GUESSP" ./restart_${res_num}/INPUT
						echo "Job $dir has ran steps :: run $dir/restart_${res_num} to continue"
					fi 
				fi
			fi
		else
			echo "Job $dir has not started"
		fi
	done
}
