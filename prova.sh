# get first number of the file 'src/simple_2DNS/status.prm'
#!/bin/bash

# Get the first line of the file
first_line=$(head -n 1 src/simple_2DNS/status.prm)

# Extract the first string
first_string=$(echo "$first_line" | awk '{print $1}')

# # Convert the string to an integer
# n=$(printf "%.0f" "$first_string")

echo "Integer value: $n"

echo $n
if [ $n -eq 0 ]; then
  echo "n is 0"
else
  echo "n is not 0"
fi