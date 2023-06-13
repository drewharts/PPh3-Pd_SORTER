#!/bin/bash

# Prompt user for source directory
read -p "Enter the source directory path: " source_dir

# Check if source directory exists
if [ ! -d "$source_dir" ]; then
  echo "Source directory does not exist."
  exit 1
fi

# Prompt user for destination directory
read -p "Enter the destination directory path: " destination_dir

# Check if destination directory exists
if [ ! -d "$destination_dir" ]; then
  echo "Destination directory does not exist."
  exit 1
fi

# Move files with ".opt.xyz" extension from source to destination directory
find "$source_dir" -type f -name "*.opt.xyz" -exec mv {} "$destination_dir" \;

echo "Files moved successfully."