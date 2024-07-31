import os
import re
import csv
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import sys

# Directory containing the zip files and output directories
def run_profiling(zip_dir, runs_dir):
    ZIP_DIR = zip_dir

    # Output CSV file
    OUTPUT_CSV = "workflow_summary.csv"

    task_cpu_allocation = {
        "bwaMem": 8,
        "bamStat": 1,
        "flowcell_metrics": 1,
        "dimer_metrics": 1,
        "combine_dimer_metrics": 1,
        "runKraken": 14,
        "consolidate_flowcell_metrics": 1,
        "combine_images": 1,
        "consolidate_samples": 1,
        "fastqc": 3,
        "multiqc": 1,
        "fastaIndex": 1,
        # Add more tasks and their CPU allocations as needed
    }

    # Regex patterns for parsing the workflow.log
    start_pattern = re.compile(r'(\d+-\d+-\d+ \d+:\d+:\d+\.\d+) .*? ready :: job: "(.*?)"')
    finish_pattern = re.compile(r'(\d+-\d+-\d+ \d+:\d+:\d+\.\d+) .*? finish :: job: "(.*?)"')
    datetime_pattern = re.compile(r'(\d+-\d+-\d+ \d+:\d+:\d+\.\d+)')
    # Function to parse datetime from the log
    def parse_datetime(dt_str):
        return datetime.strptime(dt_str, '%Y-%m-%d %H:%M:%S.%f')

    # Collect data for CSV
    data = []

    # Iterate over each output directory
    for item in os.listdir(runs_dir):
        item_path = os.path.join(runs_dir, item)
        if os.path.isdir(item_path):
            workflow_log_path = os.path.join(item_path, '_LAST', "workflow.log")
            if os.path.isfile(workflow_log_path):
                print(f"Processing {item_path}")

                # Extract the base name of the directory as the run_name
                run_name = item

                # Get the corresponding zip file size
                zip_file_path = os.path.join(ZIP_DIR, f"{run_name}.zip")
                if os.path.isfile(zip_file_path):
                    zip_file_size = os.path.getsize(zip_file_path)
                else:
                    zip_file_size = 0

                # Initialize dictionaries to track start and end times
                start_times = {}
                end_times = {}
                timestamps = []

                # Parse the workflow.log to extract task runtimes
                with open(workflow_log_path, 'r') as log_file:
                    for line in log_file:
                        start_match = start_pattern.search(line)
                        if start_match:
                            dt_str, task_name = start_match.groups()
                            start_times[task_name] = parse_datetime(dt_str)

                        finish_match = finish_pattern.search(line)
                        if finish_match:
                            dt_str, task_name = finish_match.groups()
                            end_times[task_name] = parse_datetime(dt_str)
                        datetime_match = datetime_pattern.search(line)
                        if datetime_match:
                            dt_str = datetime_match.group(1)
                            timestamps.append(parse_datetime(dt_str))

                # Calculate total runtime for the workflow
                total_cpu_seconds = 0
                for task_name in start_times:
                    if task_name in end_times:
                        runtime_seconds = (end_times[task_name] - start_times[task_name]).total_seconds()

                        # Get CPU allocation for the task
                        cpus = task_cpu_allocation.get(task_name, 1)  # Default to 1 if not specified
                        cpu_seconds = runtime_seconds * cpus
                        total_cpu_seconds += cpu_seconds
                if timestamps:
                    start_time = min(timestamps)
                    end_time = max(timestamps)
                    total_runtime_seconds = (end_time - start_time).total_seconds()
                else:
                    total_runtime_seconds = 0
                # Append the information to the data list
                data.append([run_name, zip_file_size, total_runtime_seconds, total_cpu_seconds])

                # Append the information to the data list
            else:
                print(f"No workflow.log found in {item_path}")

    # Write data to CSV
    with open(OUTPUT_CSV, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(["run_name", "zip_file_size", "total_runtime_seconds", "total_cpu_seconds"])
        csv_writer.writerows(data)

    print(f"Summary written to {OUTPUT_CSV}")

    # Load the data from CSV
    df = pd.read_csv(OUTPUT_CSV)

    # Plot zip file size against total runtime
    plt.scatter(df['zip_file_size'], df['total_cpu_seconds'])
    plt.xlabel('Zip File Size (bytes)')
    plt.ylabel('Total Runtime (cpu seconds)')
    plt.title('Total Runtime vs. Zip File Size')
    plt.grid(True)
    plt.savefig('runtime_vs_zip_size.png')
    plt.show()


if __name__ == "__main__":
    run_profiling(sys.argv[1], sys.argv[2])