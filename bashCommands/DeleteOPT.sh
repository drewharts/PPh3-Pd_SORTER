#!/bin/bash

# Prompt the user to input the directory
read -p "Enter the directory path: " directory

# Check if the directory exists
if [ ! -d "$directory" ]; then
  echo "Directory not found!"
  exit 1
fi

# Navigate to the target directory
cd "$directory"

# Find and delete files with the ".opt.xyz" extension
find . -type f -name "*.opt.xyz" -delete

# Print completion message
echo "Deletion complete."
