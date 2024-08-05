import os
import re
import csv
from datetime import datetime
import pandas as pd
import seaborn as sns
from scipy.stats import linregress
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
    task_durations = []

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
                        task_name = task_name.split("-")[1]

                        # Get CPU allocation for the task
                        cpus = task_cpu_allocation.get(task_name, 1)  # Default to 1 if not specified
                        cpu_seconds = runtime_seconds * cpus
                        total_cpu_seconds += cpu_seconds

                        # Store task durations for additional analysis
                        task_durations.append((run_name, task_name, runtime_seconds, cpu_seconds))

                if timestamps:
                    start_time = min(timestamps)
                    end_time = max(timestamps)
                    total_runtime_seconds = (end_time - start_time).total_seconds()
                else:
                    total_runtime_seconds = 0
                # Append the information to the data list
                data.append([run_name, zip_file_size, total_runtime_seconds, total_cpu_seconds])

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
    df['zip_file_size_kb'] = df['zip_file_size'] / 1024 / 1024  # Convert bytes to MB

    # Initialize the plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Plot for total CPU seconds
    sns.regplot(x='zip_file_size_kb', y='total_cpu_seconds', data=df, ax=axes[0])
    axes[0].set_xlabel('Zip File Size (MB)')
    axes[0].set_ylabel('Total CPU Runtime (cpu seconds)')
    axes[0].set_title('Total CPU Runtime vs. Zip File Size')
    slope_cpu, intercept_cpu, r_value_cpu, p_value_cpu, std_err_cpu = linregress(df['zip_file_size_kb'], df['total_cpu_seconds'])
    axes[0].text(0.05, 0.95, f'Slope: {slope_cpu:.2f} S/MB', transform=axes[0].transAxes, fontsize=12, verticalalignment='top')

    # Plot for total runtime seconds
    sns.regplot(x='zip_file_size_kb', y='total_runtime_seconds', data=df, ax=axes[1])
    axes[1].set_xlabel('Zip File Size (MB)')
    axes[1].set_ylabel('Total Runtime (seconds)')
    axes[1].set_title('Total Runtime vs. Zip File Size')
    slope_runtime, intercept_runtime, r_value_runtime, p_value_runtime, std_err_runtime = linregress(df['zip_file_size_kb'], df['total_runtime_seconds'])
    axes[1].text(0.05, 0.95, f'Slope: {slope_runtime:.2f} S/MB', transform=axes[1].transAxes, fontsize=12, verticalalignment='top')

    plt.tight_layout()
    plt.savefig('runtime_vs_zip_size_subplot.png')

    # Plot total runtime vs total CPU runtime
    fig, ax = plt.subplots()
    sns.regplot(x='total_runtime_seconds', y='total_cpu_seconds', data=df, ax=ax)
    ax.set_xlabel('Total Runtime (seconds)')
    ax.set_ylabel('Total CPU Runtime (cpu seconds)')
    ax.set_title('Total CPU Runtime vs. Total Runtime')
    plt.savefig('total_cpu_vs_total_runtime.png')
    plt.show()

    # Analyze task durations
    task_df = pd.DataFrame(task_durations, columns=["run_name", "task_name", "runtime_seconds", "cpu_seconds"])

    # Plot distribution of task runtimes
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='task_name', y='runtime_seconds', data=task_df)
    plt.xticks(rotation=90)
    plt.xlabel('Task Name')
    plt.ylabel('Runtime (seconds)')
    plt.title('Distribution of Task Runtimes')
    plt.tight_layout()
    plt.savefig('task_runtime_distribution.png')

    # Plot CPU seconds per task
    plt.figure(figsize=(10, 6))
    sns.boxplot(y='task_name', x='cpu_seconds', data=task_df)
    plt.xticks(rotation=90)
    plt.ylabel('Task Name')
    plt.xlabel('CPU Time (seconds)')
    plt.title('Distribution of CPU Time per Task')
    plt.tight_layout()
    plt.savefig('task_cpu_time_distribution.png')
    plt.show()

    # # Identify bottlenecks: tasks with high average runtime
    # avg_runtime_df = task_df.groupby('task_name')['runtime_seconds'].mean().reset_index()
    # avg_runtime_df = avg_runtime_df.sort_values(by='runtime_seconds', ascending=False)
    #
    # plt.figure(figsize=(10, 6))
    # sns.barplot(x='runtime_seconds', y='task_name', data=avg_runtime_df)
    # plt.xlabel('Average Runtime (seconds)')
    # plt.ylabel('Task Name')
    # plt.title('Average Runtime per Task')
    # plt.tight_layout()
    # plt.savefig('average_runtime_per_task.png')
    # plt.show()

if __name__ == "__main__":
    run_profiling(sys.argv[1], sys.argv[2])

#%%
