test -s ~/.alias && . ~/.alias || true
set -o vi
alias ll='ls -lrt'
alias wrk='cd /work/e05/e05/ivanmres/'
for file in $(ls ~/bash_scripts/*.sh ); do # source files to all commands
	source $file
	for funct in $( grep function $file | awk '{print $2}' | cut -d '(' -f1 ); do
		export -f $funct
	done
done
