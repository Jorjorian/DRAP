version development

task flowcell_metrics {
    input {
        File fq_1
        File fq_2
        String tag
    }
    command <<<
        pip install sklearn
        python <<CODE
import os
import json
from Bio import SeqIO
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
import numpy as np

def extract_xy(file):
    x_coords = []
    y_coords = []
    for record in SeqIO.parse(file, "fastq"):
        # Assuming the x, y coordinates are stored in the read name as 'x_y'
        name = record.id
        x, y = map(int, name.split(":")[-2:])
        x_coords.append(x)
        y_coords.append(y)
return x_coords, y_coords

def plot_xy(x_coords_1, y_coords_1, x_coords_2, y_coords_2, tag):
    plt.figure(figsize=(10, 6))
    plt.scatter(x_coords_1, y_coords_1, s=1, alpha=0.5, label='Read 1')
    plt.scatter(x_coords_2, y_coords_2, s=1, alpha=0.5, label='Read 2')
    plt.title(tag + " - Combined Strand Distribution")
    plt.xlabel("X coordinate")
    plt.ylabel("Y coordinate")
    plt.legend()
    plt.savefig(tag + "_combined_flowcell_distribution.png")
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
        plt.scatter(xy[:, 0], xy[:, 1], s=1, alpha=0.5, color=tuple(col))

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
plot_xy(x_coords_1, y_coords_1, x_coords_2, y_coords_2, "~{tag}")

# Clustering and plotting
all_x_coords = x_coords_1 + x_coords_2
all_y_coords = y_coords_1 + y_coords_2
cluster_and_plot(all_x_coords, all_y_coords, "~{tag}")
CODE
        >>>
output {
        File metrics = "~{tag}.flowcell_metrics.json"
        File overall_png = glob("~{tag}_flowcell_distribution.png")
        File clusters = glob("~{tag}_clusters.png")
    }
    runtime {
        docker: "biopython/biopython"
    }
}