#!/bin/bash

module purge
module load anaconda/4.12.0

INFILE=$1
ITER=$2
CONFIG_NAME=$3

if [ -z "$ITER" ]; then
  ITER=1
fi

COUNT=0
while [ $COUNT -lt $ITER ]
do
    USER=$(whoami)
    CURRENTPATH=$(pwd)
    DATEDAY=$(date +%d)
    DATEMONTH=$(date +%m)
    DATEYEAR=$(date +%Y)
    DATEHOUR=$(date +%H)
    DATEMIN=$(date +%M)
    DATE=${DATEYEAR}${DATEMONTH}${DATEDAY}
    UNAME=$(uname)

    DIR=${CURRENTPATH}/output/${CONFIG_NAME}/${INFILE}-${COUNT}
    if [ -d $DIR ]; then
       echo "Directory exists, removing and recreating $DIR"
       rm -rf $DIR
    fi

    mkdir -p $DIR
    cp $INFILE $DIR
    cp ${CONFIG_NAME}.txt $DIR

    SCRIPT=$DIR/run_${INFILE}-${CONFIG_NAME}-${COUNT}.csh

    cat << EOF > $SCRIPT
#!/bin/bash
#BSUB -J ${INFILE}-${CONFIG_NAME}-${COUNT}
#BSUB -q normal
#BSUB -r
#BSUB -C 0
#BSUB -n 8
#BSUB -R "rusage[mem=2500]"
#BSUB -Q "140"
cd $DIR

start_time=\$(date +%s.%N)

conda activate amberenv
python $INFILE $CONFIG_NAME

end_time=\$(date +%s.%N)
runtime=\$(echo "\$end_time - \$start_time" | bc)

echo "Script took \${runtime} seconds to run."

EOF

   chmod +x $SCRIPT

   bsub  -e $DIR/log.err -o $DIR/log.out < $SCRIPT

   COUNT=$(expr $COUNT + 1)

done