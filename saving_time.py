import argparse
import sys
import numpy as np

if __name__ == '__main__':
    if len(sys.argv) == 3:
        FILENAME = sys.argv[1]
        SAVEFILENAME = sys.argv[2]
        
        # Initialize an empty list to store the extracted values
        time_values = []

        # Open the file and read its contents
        with open(FILENAME, 'r') as file:
            for line in file:
                # Find the index of "Time: " in the line
                time_index = line.find("Time: ")
                if time_index != -1:
                    # Extract the value after "Time: " (including " s")
                    time_str = line[time_index + len("Time: "):].strip()
                    try:
                        # Remove the " s" and convert the value to a float
                        time_value = float(time_str[:-2])
                        # Append the converted value to the list
                        time_values.append(time_value)
                    except ValueError:
                        # Handle the case where the value after "Time: " is not a valid float
                        print(f"Invalid time value: {time_str}")

        # Check if any time values were found
        if len(time_values) > 0:
            # Convert the list of values to a numpy array
            time_array = np.array(time_values)

            # Calculate mean, max, and min times
            mean_time = np.mean(time_array)
            max_time = np.max(time_array)
            min_time = np.min(time_array)

            # Print the results
            print(f"Mean Time: {mean_time}")
            print(f"Max Time: {max_time}")
            print(f"Min Time: {min_time}")

            # Save the time values to a TSV file
            np.savetxt(SAVEFILENAME, time_array, delimiter='\t', header='Time (s)', comments='', fmt='%.6f')
            print(f"Time array saved to '{SAVEFILENAME}'")
        else:
            print("No valid time values found in the file.")
    else:
        print("Usage: python script_name.py <input_file.txt> <saving_filename>")

