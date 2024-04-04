#!/bin/bash

dataset_list=("zeisel" "pbmc3k")

for dataset in ${dataset_list[@]}; do

OUTDIR="./results/simulated/"

LOGDIR="${OUTDIR}/log/${dataset}"
SHDIR="${OUTDIR}/sh/${dataset}"
BUDIR="${SHDIR}/.bu"
mkdir -p $BUDIR


files_mem_err=$(find "$LOGDIR" -maxdepth 1 -type f -name "*.log" -print0 | xargs -0 grep -rl "MemoryError")
files_time_err=$(find "$LOGDIR" -maxdepth 1 -type f -name "*.log" -print0 | xargs -0 grep -rl "ERROR_TIMEOUT")
files_alloc_err=$(find "$LOGDIR" -maxdepth 1 -type f -name "*.log" -print0 | xargs -0 grep -rl "bad_alloc")

files_err=$(find "$LOGDIR" -maxdepth 1 -type f -name "*.log" -print0 | xargs -0 grep -rl "Error")
files_err=$(comm -23 <(echo "$files_err" | sort) <(echo "$files_mem_err" | sort))
files_err=$(comm -23 <(echo "$files_err" | sort) <(echo "$files_time_err" | sort))
files_err=$(comm -23 <(echo "$files_err" | sort) <(echo "$files_alloc_err" | sort))

# MemoryError
# ERROR_TIMEOUT

TMPDIR=0G

MAXMEM=500G
MAXTIME=1440
MAXDIR=200G

IFS=$'\n'

##################
## MEMORY ERROR ##
##################

for file in $files_mem_err; do
    bn=$(basename -s ".stdout.log" $file)

    sh_file=$SHDIR/$bn.sh

    # Extract metadata

    while read -r line; do
        if [[ "$line" =~ ^#\ (threads)=(.*)$ ]]; then
            THREADS="${BASH_REMATCH[2]}"
            echo "THREATS: $THREADS"
        fi

        if [[ "$line" =~ ^#\ (memory)=(.*)$ ]]; then
            MEMORY="${BASH_REMATCH[2]}"
            echo "MEMORY: $MEMORY"
        fi

        if [[ "$line" =~ ^#\ (tmpdir)=(.*)$ ]]; then
            TMPDIR="${BASH_REMATCH[2]}"
            echo "TMPDIR: $TMPDIR"
        fi

        if [[ "$line" =~ ^#\ (t)=(.*)$ ]]; then
            MINUTES="${BASH_REMATCH[2]}"
            echo "MINUTES: $MINUTES"
        fi

    done < <(awk '/# BEGIN_MXQ/,/# END_MXQ/' "$sh_file")

    MEMORY=$MAXMEM

#if [ "$MEMORY" -lt "$MAXMEM"]; then
#	MEMORY=$((MEMORY*2))
#else
#	MEMORY=$MAXMEM



META_MXQ=$(cat << EOM
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# tmpdir=$TMPDIR
# t=$MINUTES
# END_MXQ
EOM
)

	cp $sh_file $BUDIR/

    echo "$META_MXQ" | cat - <(sed '1,/# END_MXQ/d' $sh_file) > temp.txt && mv temp.txt $sh_file
    chmod +x $sh_file

    echo "Processing file: $file"

    mxqsub --stdout=$file \
           --group-name="MEMORY_error_failed_runs_${dataset}" \
           --threads=$THREADS \
           --memory=$MEMORY \
           --tmpdir=$TMPDIR \
	   -t $MINUTES \
           bash $sh_file

    # More processing code here...
done

################
## TIME ERROR ##
################


for file in $files_time_err; do
    bn=$(basename -s ".stdout.log" $file)

    sh_file=$SHDIR/$bn.sh

    # Extract metadata

    while read -r line; do
        if [[ "$line" =~ ^#\ (threads)=(.*)$ ]]; then
            THREADS="${BASH_REMATCH[2]}"
            echo "THREATS: $THREADS"
        fi

        if [[ "$line" =~ ^#\ (memory)=(.*)$ ]]; then
            MEMORY="${BASH_REMATCH[2]}"
            echo "MEMORY: $MEMORY"
        fi

        if [[ "$line" =~ ^#\ (tmpdir)=(.*)$ ]]; then
            TMPDIR="${BASH_REMATCH[2]}"
            echo "TMPDIR: $TMPDIR"
        fi

        if [[ "$line" =~ ^#\ (t)=(.*)$ ]]; then
            MINUTES="${BASH_REMATCH[2]}"
            echo "MINUTES: $MINUTES"
        fi

    done < <(awk '/# BEGIN_MXQ/,/# END_MXQ/' "$sh_file")

    MINUTES=$MAXTIME


META_MXQ=$(cat << EOM
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# tmpdir=$TMPDIR
# t=$MINUTES
# END_MXQ
EOM
)

	cp $sh_file $BUDIR/

    echo "$META_MXQ" | cat - <(sed '1,/# END_MXQ/d' $sh_file) > temp.txt && mv temp.txt $sh_file
    chmod +x $sh_file

    echo "Processing file: $file"

    mxqsub --stdout=$file \
           --group-name="TIME_error_failed_runs_${dataset}" \
           --threads=$THREADS \
           --memory=$MEMORY \
           --tmpdir=$TMPDIR \
	   -t $MINUTES \
           bash $sh_file

done


#################
## ALLOC ERROR ##
#################


for file in $files_alloc_err; do
    bn=$(basename -s ".stdout.log" $file)

    sh_file=$SHDIR/$bn.sh

    # Extract metadata

    while read -r line; do
        if [[ "$line" =~ ^#\ (threads)=(.*)$ ]]; then
            THREADS="${BASH_REMATCH[2]}"
            echo "THREATS: $THREADS"
        fi

        if [[ "$line" =~ ^#\ (memory)=(.*)$ ]]; then
            MEMORY="${BASH_REMATCH[2]}"
            echo "MEMORY: $MEMORY"
        fi

        if [[ "$line" =~ ^#\ (tmpdir)=(.*)$ ]]; then
            TMPDIR="${BASH_REMATCH[2]}"
            echo "TMPDIR: $TMPDIR"
        fi

        if [[ "$line" =~ ^#\ (t)=(.*)$ ]]; then
            MINUTES="${BASH_REMATCH[2]}"
            echo "MINUTES: $MINUTES"
        fi

    done < <(awk '/# BEGIN_MXQ/,/# END_MXQ/' "$sh_file")

   MEMORY=$MAXMEM
   TMPDIR=200G


META_MXQ=$(cat << EOM
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# tmpdir=$TMPDIR
# t=$MINUTES
# END_MXQ
EOM
)

	cp $sh_file $BUDIR/

    echo "$META_MXQ" | cat - <(sed '1,/# END_MXQ/d' $sh_file) > temp.txt && mv temp.txt $sh_file
    chmod +x $sh_file

    echo "Processing file: $file"

    mxqsub --stdout=$file \
           --group-name="BADALLOC_failed_runs_${dataset}" \
           --threads=$THREADS \
           --memory=$MEMORY \
           --tmpdir=$TMPDIR \
	   -t $MINUTES \
           bash $sh_file

done


###################
## GENERAL ERROR ##
###################


for file in $files_err; do
    bn=$(basename -s ".stdout.log" $file)

    sh_file=$SHDIR/$bn.sh

    # Extract metadata

    while read -r line; do
        if [[ "$line" =~ ^#\ (threads)=(.*)$ ]]; then
            THREADS="${BASH_REMATCH[2]}"
            echo "THREATS: $THREADS"
        fi

        if [[ "$line" =~ ^#\ (memory)=(.*)$ ]]; then
            MEMORY="${BASH_REMATCH[2]}"
            echo "MEMORY: $MEMORY"
        fi

        if [[ "$line" =~ ^#\ (tmpdir)=(.*)$ ]]; then
            TMPDIR="${BASH_REMATCH[2]}"
            echo "TMPDIR: $TMPDIR"
        fi

        if [[ "$line" =~ ^#\ (t)=(.*)$ ]]; then
            MINUTES="${BASH_REMATCH[2]}"
            echo "MINUTES: $MINUTES"
        fi

    done < <(awk '/# BEGIN_MXQ/,/# END_MXQ/' "$sh_file")

    halfmem=$((${MAXMEM%G}/2))
    halftime=$((MAXTIME/2))
    halfdir=$((${MAXDIR%G}/2))

    if [ "${MEMORY%G}" -lt "$halfmem" ]; then
    	MEMORY=$((${MEMORY%G}*2))"G"
    else
    	MEMORY=$MAXMEM
    fi
    
    if [ "$MINUTES" -lt "$halftime" ]; then
    	MINUTES=$((MINUTES*2))
    else
    	MINUTES=$MAXTIME
    fi
    
    if [ "${TMPDIR%G}" -eq "0" ]; then 
       	TMPDIR=50G
    elif [ "${TMPDIR%G}" -lt "$halfdir" ]; then
    	TMPDIR=$((${TMPDIR%G} * 2))"G"
    else 
    	TMPDIR=$MAXDIR
    fi


META_MXQ=$(cat << EOM
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# tmpdir=$TMPDIR
# t=$MINUTES
# END_MXQ
EOM
)

	cp $sh_file $BUDIR/

    echo "$META_MXQ" | cat - <(sed '1,/# END_MXQ/d' $sh_file) > temp.txt && mv temp.txt $sh_file
    chmod +x $sh_file

    echo "Processing file: $file"

    mxqsub --stdout=$file \
           --group-name="ERROR_failed_runs_${dataset}" \
           --threads=$THREADS \
           --memory=$MEMORY \
           --tmpdir=$TMPDIR \
	   -t $MINUTES \
           bash $sh_file

done

###################
## Missing sh    ##
###################

#sh_files=$(find $SHDIR -type f -name "*.sh")
shtorun=$(find $SHDIR -maxdepth 1 -type f -name "*.sh")
all_files="$files_err\n$files_mem_err\n$files_time_err\n$files_alloc_err"

#shtorun=$(comm -23 <(echo "$sh_files" | sort) <(echo "$all_files" | sort))

for file in $shtorun; do
 
    bn=$(basename -s ".sh" $file)
    
    log_file="${LOGDIR}/${bn}.stdout.log"

	if [[ $all_files == *$log_file* ]]; then
        continue    
    fi
    
    if [[ -f "$log_file" ]]; then
        echo "${log_file} exists."
    else
        echo "ALARM: ${log_file} does not exist."
    fi

    sh_file=$SHDIR/$bn.sh

    # Extract metadata

    while read -r line; do
        if [[ "$line" =~ ^#\ (threads)=(.*)$ ]]; then
            THREADS="${BASH_REMATCH[2]}"
            echo "THREATS: $THREADS"
        fi

        if [[ "$line" =~ ^#\ (memory)=(.*)$ ]]; then
            MEMORY="${BASH_REMATCH[2]}"
            echo "MEMORY: $MEMORY"
        fi

        if [[ "$line" =~ ^#\ (tmpdir)=(.*)$ ]]; then
            TMPDIR="${BASH_REMATCH[2]}"
            echo "TMPDIR: $TMPDIR"
        fi

        if [[ "$line" =~ ^#\ (t)=(.*)$ ]]; then
            MINUTES="${BASH_REMATCH[2]}"
            echo "MINUTES: $MINUTES"
        fi

    done < <(awk '/# BEGIN_MXQ/,/# END_MXQ/' "$sh_file")

    halfmem=$((${MAXMEM%G}/2))
    halftime=$((MAXTIME/2))
    halfdir=$((${MAXDIR%G}/2))

    if [ "${MEMORY%G}" -lt "$halfmem" ]; then
        MEMORY=$((${MEMORY%G}*2))"G"
    else
        MEMORY=$MAXMEM
    fi
    
    if [ "$MINUTES" -lt "$halftime" ]; then
        MINUTES=$((MINUTES*2))
    else
        MINUTES=$MAXTIME
    fi
    
    if [ "${TMPDIR%G}" -eq "0" ]; then 
        TMPDIR=50G
    elif [ "${TMPDIR%G}" -lt "$halfdir" ]; then
        TMPDIR=$((${TMPDIR%G} * 2))"G"
    else 
        TMPDIR=$MAXDIR
    fi


META_MXQ=$(cat << EOM
#!/bin/bash

# BEGIN_MXQ
# threads=$THREADS
# memory=$MEMORY
# tmpdir=$TMPDIR
# t=$MINUTES
# END_MXQ
EOM
)

    cp $sh_file $BUDIR/

    echo "$META_MXQ" | cat - <(sed '1,/# END_MXQ/d' $sh_file) > temp.txt && mv temp.txt $sh_file
    chmod +x $sh_file

    echo "Processing file: $file"

    mxqsub --stdout=$log_file \
           --group-name="MISSING_runs_${dataset}" \
           --threads=$THREADS \
           --memory=$MEMORY \
           --tmpdir=$TMPDIR \
       -t $MINUTES \
           bash $sh_file

done


done
