#!/bin/bash -e

INFILE=$1

# Create a new directory for the output
OUTPUT_DIR=output/$(date +"%Y%m%d_%H%M%S")

mkdir -p $OUTPUT_DIR

# Copy input files to the output directory
cp $INFILE $OUTPUT_DIR

# Generate the script
SCRIPT=${OUTPUT_DIR}/run_${INFILE}.sh

cat << EOF > $SCRIPT
#!/bin/bash
#BSUB -J ${INFILE}
#BSUB -q normal
#BSUB -r
#BSUB -C 0
#BSUB -n 1
#BSUB -R "rusage[mem=2500]"
#BSUB -Q "140"
cd $OUTPUT_DIR

start_time=\$(date +%s.%N)

conda activate amberenv
python $INFILE

end_time=\$(date +%s.%N)
runtime=\$(echo "\$end_time - \$start_time" | bc)

echo "Script took \${runtime} seconds to run."

EOF

chmod +x $SCRIPT
bsub -e $OUTPUT_DIR/log.err -o $OUTPUT_DIR/log.out < $SCRIPT
