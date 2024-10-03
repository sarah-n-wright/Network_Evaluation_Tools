#!/bin/bash

# Function to process the list of files
process_files() {
    input_file="$1"

    # Loop through each file listed in the input file
    while IFS= read -r file; do
        if [ -f "$file" ]; then
            # Read the first line of the file
            first_line=$(head -n 1 "$file")
            
            # Extract the chromosome from the first line (assuming it's the first column)
            chromosome=$(echo "$first_line" | awk '{print $1}')
            
            # Write the file name to the corresponding chromosome file
            echo "$file" >> "${chromosome}_files.txt"
        else
            echo "File $file not found."
        fi
    done < "$input_file"
}

# Example usage
process_files $1

