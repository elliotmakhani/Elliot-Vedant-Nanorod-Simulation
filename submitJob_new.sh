#! /bin/bash
NAME=nanoAD

# Checks if a positional argument for ORDER was given; if not, then it prompts user for it
if [[ -n "$1" ]]; then
  ORDER=$1
else
  read -p "submit jobs? y or n"$'\n' ORDER
fi


if [ "$ORDER" == y ]; then

    # Checks if a positional argument for number of runs was given; if not, then it prompts user for it
    if [[ -n "$2" ]]; then
        RUNS=$2
    else
        read -p "how many runs? (max 1000)"$'\n' RUNS
    fi

    # Obtains full path to submission script's directory and stores in FILE_DIR
    FILE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

    # Stores the path to the compiled simulation program in SIM_EXE
    SIM_EXE="$FILE_DIR"/"$NAME"

    # Navigates to the simulation folder and compiles program, moving it to the parent directory
    cd $FILE_DIR'/nanoAD-full/build'
    cmake ..
    cmake --build .
    mv $FILE_DIR'/nanoAD-full/build/'$NAME $SIM_EXE
    cd $FILE_DIR

    # Remove all folders whose names are only digits
    for dir in $FILE_DIR/*; do
        [[ -d $dir ]] && [[ $dir =~ ^$FILE_DIR/[0-9]+$ ]] && rm -r $dir; 
    done

    # Iterates for each run, generating a folder for each, an Slurm job script and a copy of the sim executable and runs it
    for x in `seq -f "%04g" 1 $RUNS`; do
        RUN_DIR="$FILE_DIR/$x"
        SUB_FILE="$RUN_DIR/$x.sub"
        
        mkdir $RUN_DIR
        cd $RUN_DIR

        echo '#!/bin/bash' > $SUB_FILE
        echo '#SBATCH -J '_$(basename $FILE_DIR)'_'$x >> $SUB_FILE     # used to be named 'n200_21_'$x'_'$NAME'
        echo '#SBATCH -A haunj_lab' >> $SUB_FILE            # For use once lab acct is set up
        echo '#SBATCH -p standard' >> $SUB_FILE             # Set to 'free' or 'standard' depending on the partition/queue desired
        echo '#SBATCH -D '$RUN_DIR >> $SUB_FILE             # Set the run directory as the working directory
        echo '#SBATCH -t 34:00:00' >> $SUB_FILE   		    # Set the time limit to 7 days
        #echo '#SBATCH --mail-type=end' >> $SUB_FILE
        #echo '#SBATCH --mail-user=emakhani@uci.edu' >> $SUB_FILE
        echo '#SBATCH --requeue' >> $SUB_FILE               # Allows requeuing
        echo $SIM_EXE >> $SUB_FILE

        sbatch $SUB_FILE
        
		cd $FILE_DIR
    done
fi
