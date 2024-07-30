version development

task flowcell_metrics {
    input {
        File fq_1
        File fq_2
        String tag
    }
    command <<<
        source activate base
        conda activate myenv
        python <<CODE
import os
import json
import gzip
from Bio import SeqIO
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
import numpy as np
def extract_xy(file):
    x_coords = []
    y_coords = []
    # Handle compressed fastq files
    open_func = gzip.open if file.endswith('.gz') else open
    with open_func(file, 'rt') as handle:
        for record in SeqIO.parse(handle, "fastq"):
            name = record.id
            name = name.split(' ')[0]
            x, y = map(int, name.split(":")[-2:])
            x_coords.append(x)
            y_coords.append(y)
    return x_coords, y_coords

def plot_xy(x_coords_1, y_coords_1, tag):
    plt.figure(figsize=(10, 6))
    plt.scatter(x_coords_1, y_coords_1, s=.01, alpha=0.5, label='Read 1')
    # save the plot
    plt.title(tag + " - Strand 1 Distribution")
    plt.xlabel("X coordinate")
    plt.ylabel("Y coordinate")
    plt.legend()
    plt.savefig(tag + "_flowcell_distribution.png")
    plt.close()


def cluster_and_plot(x_coords, y_coords, tag):
    coords = np.column_stack((x_coords, y_coords))
    db = DBSCAN(eps=3, min_samples=10).fit(coords)
    labels = db.labels_

    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]

    plt.figure(figsize=(10, 6))
    for k, col in zip(unique_labels, colors):
        if k == -1:
            col = [0, 0, 0, 1]
        class_member_mask = (labels == k)
        xy = coords[class_member_mask]
        plt.scatter(xy[:, 0], xy[:, 1], s=.001, alpha=0.5, color=tuple(col))

    plt.title(tag + " - Clustered Distribution")
    plt.xlabel("X coordinate")
    plt.ylabel("Y coordinate")
    plt.savefig(tag + "_clusters.png")
    plt.close()

    # Metrics: number of clusters, noise points, average cluster size
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = list(labels).count(-1)
    cluster_sizes = [list(labels).count(label) for label in unique_labels if label != -1]
    avg_cluster_size = np.mean(cluster_sizes) if cluster_sizes else 0

    metrics = dict(
        n_clusters = n_clusters,
        n_noise = n_noise,
        avg_cluster_size = avg_cluster_size,
        cluster_sizes = cluster_sizes
        )

    with open(tag + ".flowcell_metrics.json", "w") as f:
        json.dump(metrics, f)

# Extracting coordinates and plotting
x_coords_1, y_coords_1 = extract_xy("~{fq_1}")
x_coords_2, y_coords_2 = extract_xy("~{fq_2}")

# Plot combined XY coordinates
plot_xy(x_coords_1, y_coords_1, "~{tag}")

# Clustering and plotting
all_x_coords = x_coords_1
all_y_coords = y_coords_1
cluster_and_plot(all_x_coords, all_y_coords, "~{tag}")
CODE
        >>>
output {
        File metrics = "~{tag}.flowcell_metrics.json"
        File overall_png = "~{tag}_flowcell_distribution.png"
        File clusters = "~{tag}_clusters.png"
    }
    runtime {
        docker: "docker.io/oblivious1/drap:flowcell"
    }
}

task consolidate_flowcell_metrics {
    input {
        Array[File] flowcell_metrics_files
    }

    command <<<

        python <<CODE
import json
import os

# Initialize data structures
consolidated_data = {
}

# Process each flowcell metrics file
for metrics_file in "~{sep=' ' flowcell_metrics_files}".split(' '):
    with open(metrics_file) as f:
        metrics = json.load(f)
    sample_name = os.path.basename(metrics_file).split('.')[0]
    consolidated_data[sample_name] = {}
    consolidated_data[sample_name]["n_clusters"] = metrics["n_clusters"]
    consolidated_data[sample_name]["n_noise"]= metrics["n_noise"]
    consolidated_data[sample_name]["avg_cluster_size"] = metrics["avg_cluster_size"]
    # Prepare the MultiQC-compatible JSON format
    multiqc_data = {
        "id": "flowcell_metrics",
        "section_name": "Flowcell Metrics",
        "description": "Flowcell metrics from multiple samples.",
        "plot_type": "table",
        "pconfig": {
            "id": "flowcell_metrics_table",
            "title": "Flowcell Metrics",
            'col1_header': "Read Group",
        },
        "data": consolidated_data
        }

        # Save the consolidated JSON
with open('multiqc_flowcell_metrics_mqc.json', 'w') as f:
    json.dump(multiqc_data, f, indent=4)
CODE
    >>>

    output {
        File multiqc_flowcell_metrics_json = "multiqc_flowcell_metrics_mqc.json"
    }

    runtime {
        memory: "4 GB"
        cpu: 2
        docker: "python:3.8-slim"
    }
}


task combine_images {
    input {
    Array[File] images
    String output_image_name
}

command <<<
    source activate base
    conda activate myenv
    python <<CODE
import math
import matplotlib.pyplot as plt
from PIL import Image

def optimal_grid(n):
    ncols = int(math.sqrt(n))
    nrows = math.ceil(n / ncols)
    return nrows, ncols

def combine_images(image_files, output_image):
    # Load images
    images = [Image.open(img) for img in image_files]

    # Determine grid size
    n = len(images)
    nrows, ncols = optimal_grid(n)

    # Determine size of the combined image
    widths, heights = zip(*(i.size for i in images))
    max_width = max(widths)
    max_height = max(heights)

    # Create a blank canvas for the combined image
    combined_image = Image.new('RGB', (ncols * max_width, nrows * max_height))

    # Paste each image into the combined image
    for idx, img in enumerate(images):
        row = idx // ncols
        col = idx % ncols
        combined_image.paste(img, (col * max_width, row * max_height))

    # Save the combined image
    combined_image.save(output_image)

    # Input images and output image name
image_files = "~{sep=' ' images}".split(' ')
output_image = "~{output_image_name}"

combine_images(image_files, output_image)
CODE
>>>

output {
        File combined_image = "~{output_image_name}"
        }

runtime {
  docker: "docker.io/oblivious1/drap:flowcell"
  memory: "4 GB"}
}