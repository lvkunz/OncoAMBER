#!/bin/bash -e

module purge
module load anaconda/4.12.0
module load topas/3.8.0
#show the version of oncoamber installed on amberenv

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
DATESEC=$(date +%S)
DATE=${DATEYEAR}${DATEMONTH}${DATEDAY}
UNAME=$(uname)

# Create a new directory for the output
DAY_DIR=${CURRENTPATH}/output/${DATE}_${USER}_${UNAME}
OUTPUT_DIR=${DAY_DIR}/${CONFIG_NAME}_${INFILE}_${DATEHOUR}${DATEMIN}${DATESEC}

if [ -d $OUTPUT_DIR ]; then
   echo "Output directory exists, removing and recreating $OUTPUT_DIR"
   rm -rf $OUTPUT_DIR
fi
mkdir -p $OUTPUT_DIR

# Copy input files to the output directory
cp $INFILE $OUTPUT_DIR
cp ${CONFIG_NAME}.txt $OUTPUT_DIR



param_names=()

# Loop until the user enters "none"
while true; do
  read -p "Enter a parameter name to change, or type 'none' to continue: " param_name

  # Check if the user entered "none", and exit the loop if so
  if [ "$param_name" == "none" ]; then
    break
  fi

  # Add the parameter name to the array
  param_names+=("$param_name")
done

# Print the list of parameter names
echo "Parameter names to be changed: ${param_names[@]}"

for (( COUNT=0; COUNT<$ITER; COUNT++ )); do

  # Create a new directory for the current iteration
  DIR=${OUTPUT_DIR}/iter${COUNT}
  if [ -d $DIR ]; then
    echo "Iteration directory exists, removing and recreating $DIR"
    rm -rf $DIR
  fi
  mkdir -p $DIR

  # Copy input files to the iteration directory
  cp $INFILE $DIR
  cp ${CONFIG_NAME}.txt $DIR

  cp TOPAS* $DIR
  cp FRAC* $DIR

  # Loop through each parameter name and prompt the user for a value
  declare -A param_values=()
  for param_name in "${param_names[@]}"; do
    read -p "Enter the value for iteration $COUNT for the parameter $param_name: " param_value
    param_values["$param_name"]="$param_value"
  done

  # Create the CSV file with the parameter values for this iteration
  if [ $COUNT -eq 0 ]; then
    # Create the header row for the CSV file
    echo "Iteration ${param_names[*]}" > ${OUTPUT_DIR}/${CONFIG_NAME}_${INFILE}_param_values.csv

  fi
  # Append the parameter values for this iteration to the CSV file
  echo -n "$COUNT" >> ${OUTPUT_DIR}/${CONFIG_NAME}_${INFILE}_param_values.csv
  for param_name in "${param_names[@]}"; do
    echo -n " ${param_values[$param_name]}" >> ${OUTPUT_DIR}/${CONFIG_NAME}_${INFILE}_param_values.csv
  done
  echo "" >> ${OUTPUT_DIR}/${CONFIG_NAME}_${INFILE}_param_values.csv

  # Modify the configuration file with the parameter values for this iteration
  for param_name in "${!param_values[@]}"; do
    # Backup the original file
    cp ${DIR}/${CONFIG_NAME}.txt ${DIR}/${CONFIG_NAME}.txt.bak

    # Attempt to replace the parameter
    sed -i "s/${param_name}: .*/${param_name}: ${param_values[$param_name]}/" ${DIR}/${CONFIG_NAME}.txt

    # Check if the file was changed
    if cmp -s "${DIR}/${CONFIG_NAME}.txt" "${DIR}/${CONFIG_NAME}.txt.bak"; then
      echo "Error: Parameter $param_name was not found in ${DIR}/${CONFIG_NAME}.txt"
      exit 1
    fi

    # If the file was changed, remove the backup
    rm ${DIR}/${CONFIG_NAME}.txt.bak
  done

  #change the parameter running_on_cluster to True
  sed -i "s/running_on_cluster: .*/running_on_cluster: True/" "${DIR}/${CONFIG_NAME}.txt"

  # Generate the script for the current iteration
  SCRIPT=${DIR}/run_${INFILE}-${CONFIG_NAME}-${COUNT}.csh

  cat << EOF > $SCRIPT
#!/bin/bash
#BSUB -J ${INFILE}-${CONFIG_NAME}-${COUNT}
#BSUB -q normal
#BSUB -r
#BSUB -C 0
#BSUB -n 1
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
  bsub -e $DIR/log.err -o $DIR/log.out < $SCRIPT

done