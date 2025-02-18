#!/bin/sh


seed=$7
i=1
t=0
while [ $i -le $6 ]; do
  j=`/home/jplank/cs302/Labs/Lab5/bin/sb-player $1 $2 $3 $4 $5 n n $seed | awk '{ print $NF }'`
  if [ a"$j" = a ]; then
    j=0
  fi
  t=`echo $t $j | awk '{ print $1+$2 }'`
  echo $i $j $t | awk '{ printf "Run %3d - Score: %6d  - Average %10.3lf\n", $1, $2, ($3+0.0)/($1+0.0) }'
  i=`echo $i | awk '{ print $1 + 1 }'`
  seed=`echo $seed | awk '{ print $1 + 1 }'`
done

#!/bin/bash

# Directory to store game outputs
OUTPUT_DIR="./output"
mkdir -p "$OUTPUT_DIR"

# Command to run
COMMAND="/home/jplank/cs302/Labs/Lab5/bin/sb-player 8 10 5 pbyrg /home/cmurra30/cs302/Superball/proj3/bin/sb-play n y 40"

# Run the command 50 times and save output
for i in $(seq 1 50); do
    OUTPUT_FILE="$OUTPUT_DIR/game$i.txt"
    echo "Running game $i..."
    $COMMAND > "$OUTPUT_FILE"
done

echo "All games completed. Outputs saved in $OUTPUT_DIR."
