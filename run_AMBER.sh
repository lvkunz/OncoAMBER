#!/bin/bash

module purge
module load python/3.8.2

set INFILE = $1
set ITER = $2
set CONFIG_NAME = $3

if ($ITER == "") then
  set ITER = 1
endif


set COUNT = 0
while ($COUNT < $ITER)
    set USER = `whoami`
    set CURRENTPATH = `pwd`
    set DATEDAY  = `date | awk '{print $3}'`
    set DATEMONTH = `date | awk '{print $2}'`
    set DATEYEAR = `date | awk '{print $6}'`
    set DATEHOUR  = `date | awk '{print $4}' | awk -F: '{print $1}'`
    set DATEMIN   = `date | awk '{print $4}' | awk -F: '{print $2}'`
    set DATE = $DATEYEAR$DATEMONTH$DATEDAY
    set UNAME = `uname`

    set DIR = $CURRENTPATH/output/$CONFIG_NAME/$INFILE-$COUNT
    if ( -d $DIR ) then
       echo Directory exists, removing and recreating $DIR
       rm -rf $DIR
    endif

    mkdir -p $DIR
    cp $INFILE $DIR
    cp ${CONFIG_NAME}.txt $DIR

    set SEED = `bash -c 'echo $RANDOM'`
    sed -i "s/seed:.*/seed: $SEED/g" $DIR/$INFILE

    set SCRIPT=$DIR/run_$INFILE-$CONFIG_NAME-$COUNT.csh

    cat - << EOF > $SCRIPT
#!/bin/bash
#BSUB -J $INFILE-$CONFIG_NAME-$COUNT
#BSUB -q normal
#BSUB -r
#BSUB -C 0
#BSUB -n 1
#BSUB -R "rusage[mem=2500]"
#BSUB -Q "140"
cd $DIR

start_time=\$(date +%s.%N)

conda run -n myenv python $INFILE $CONFIG_NAME

end_time=\$(date +%s.%N)
runtime=\$(echo "\$end_time - \$start_time" | bc)

echo "Script took \${runtime} seconds to run."

EOF

   chmod +x $SCRIPT

   bsub  -e $DIR/log.err -o $DIR/log.out < $SCRIPT

   @ COUNT = $COUNT + 1

end
