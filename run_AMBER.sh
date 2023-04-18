#!/bin/tcsh

module purge
module load python/3.7.0

set INFILE = $1
set ITER = $2
set CONFIG_NAME = $3

if ($ITER == "") then
  set ITER = 1
endif

set COUNT = `ls -d ./output/$CONFIG_NAME/$INFILE-* | wc -l`
if ($COUNT == "") then
  set COUNT = 0
endif

while ($COUNT < $ITER + $COUNT)
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
    cp *.py $DIR
    cp supportFiles/* $DIR

    set SEED = `bash -c 'echo $RANDOM'`
    echo i:Ts/Seed = $SEED >> $DIR/$INFILE

    set SCRIPT=$DIR/run_$INFILE-$CONFIG_NAME-$COUNT.csh

    cat - << EOF > $SCRIPT
#!/bin/bash
#BSUB -J $INFILE-$CONFIG_NAME-$COUNT
#BSUB -q long
#BSUB -r
#BSUB -C 0
#BSUB -n 1
#BSUB -R "rusage[mem=2500]"
#BSUB -Q "140"
cd $DIR

time python $INFILE $CONFIG_NAME

EOF

   chmod +x $SCRIPT

   bsub  -e $DIR/log.err -o $DIR/log.out < $SCRIPT

   @ COUNT = $COUNT + 1
end
