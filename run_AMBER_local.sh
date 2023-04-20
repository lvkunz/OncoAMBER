#!/bin/tcsh

set INFILE = $1
set ITER = $2
set CONFIG_NAME = $3

# List of viscosity values
set VIS_LIST = (50 100 200 500 1000)

if ($ITER == "") then
  set ITER = 1
endif

set COUNT = `ls -d ./output/$CONFIG_NAME/$INFILE-* | wc -l`
if ($COUNT == "") then
  set COUNT = 0
endif

set VIS_INDEX = 1  # Initialize viscosity index

while ($COUNT < $ITER + $COUNT)
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
       echo "Directory $DIR exists, removing and recreating..."
       rm -rf $DIR
    endif

    echo "Creating directory $DIR..."
    mkdir -p $DIR
    echo "Copying files into directory $DIR..."
    cp $INFILE $DIR
    cp *.py $DIR
    cp supportFiles/* $DIR

    # Set viscosity value
    set VIS = $VIS_LIST[$VIS_INDEX]
    echo "Setting viscosity to $VIS in $DIR/$INFILE..."
    sed -i "s/viscosity:.*/viscosity: $VIS/g" $DIR/$INFILE

    set SEED = `bash -c 'echo $RANDOM'`
    echo "Setting seed to $SEED in $DIR/$INFILE..."
    sed -i "s/seed:.*/seed: $SEED/g" $DIR/$INFILE

    cd $DIR
    echo "Running simulation in $DIR..."
    set LOGFILE = "log.out"
    set ERRFILE = "log.err"
    python $INFILE $CONFIG_NAME >& $LOGFILE || echo "Error" >& $ERRFILE

    # Increment viscosity index
    @ VIS_INDEX = ($VIS_INDEX % $#VIS_LIST) + 1

    @ COUNT = $COUNT + 1
    echo "Finished simulation $COUNT of $ITER."
end
