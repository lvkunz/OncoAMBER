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

USER=$(whoami)
CURRENTPATH=$(pwd)
DATEDAY=$(date +%d)
DATEMONTH=$(date +%m)
DATEYEAR=$(date +%Y)
DATEHOUR=$(date +%H)
DATEMIN=$(date +%M)
DATE=${DATEYEAR}${DATEMONTH}${DATEDAY}
UNAME=$(uname)

# Create a new directory for the output
OUTPUT_DIR=${CURRENTPATH}/output/${CONFIG_NAME}_${INFILE}_${DATE}_${USER}_${UNAME}
if [ -d $OUTPUT_DIR ]; then
   echo "Output directory exists, removing and recreating $OUTPUT_DIR"
   rm -rf $OUTPUT_DIR
fi
mkdir -p $OUTPUT_DIR

# Copy input files to the output directory
cp $INFILE $OUTPUT_DIR
cp ${CONFIG_NAME}.txt $OUTPUT_DIR

while [ $COUNT -lt $ITER ]
do
    # Create a new directory for the current iteration
    DIR=${OUTPUT_DIR}/${INFILE}_${CONFIG_NAME}_dt${DT}_iter${COUNT}
    if [ -d $DIR ]; then
       echo "Iteration directory exists, removing and recreating $DIR"
       rm -rf $DIR
    fi
    mkdir -p $DIR

    # Copy input files to the iteration directory
    cp $INFILE $DIR
    cp ${CONFIG_NAME}.txt $DIR

    # Update the configuration file with the DT value
    read -p "Enter the DT value: " DT
    sed -i "s/dt: .*/dt: ${DT}/" ${DIR}/${CONFIG_NAME}.txt

    # Ask user if they want to change another parameter
    read -p "Do you want to change another parameter? (y/n): " ANS
    if [ $ANS == "n" ]; then
        pass
    elif [ $ANS == "y" ]; then
        # Ask user for parameter name and value
        read -p "Enter the parameter name: " PARAM_NAME
        read -p "Enter the parameter value: " PARAM_VALUE
        sed -i "s/${PARAM_NAME}: .*/${PARAM_NAME}: ${PARAM_VALUE}/" ${DIR}/${CONFIG_NAME}.txt

        # Modify the directory name to include the parameter name
        DIR=${OUTPUT_DIR}/${INFILE}_${CONFIG_NAME}_dt${DT}_${PARAM_NAME}_${PARAM_VALUE}_iter${COUNT}
        if [ -d $DIR ]; then
           echo "Iteration directory exists, removing and recreating $DIR"
           rm -rf $DIR
        fi
        mkdir -p $DIR

        # Copy input files to the new iteration directory
        cp $INFILE $DIR
        cp ${CONFIG_NAME}.txt $DIR
    fi

    # Generate the script for the current iteration
    SCRIPT=${DIR}/run_${INFILE}-${CONFIG_NAME}-${COUNT}.csh

    cat << EOF > $SCRIPT
#!/bin/bash
#BSUB -J ${INFILE}-${CONFIG_NAME}-${COUNT}
#BSUB -q bigmem
#BSUB -r
#BSUB -C 0
#BSUB -n 1
#BSUB -R "rusage[mem=64000]"
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
