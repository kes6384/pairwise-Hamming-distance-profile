import os
import re
import numpy as np
import matplotlib.pyplot as plt
import math

# ============================================================
# User-configurable paths
# ============================================================

# Path to the ground-truth file
GROUND_TRUTH_FILE = "/Users/kesil/Desktop/Undergraduate_Thesis/pairwise-Hamming-distance-profile/test/results_fullProfile_10000_32.txt"

# Directory containing rep_1.txt, rep_2.txt, ..., rep_1000.txt
REPLICATE_DIR = "/Users/kesil/Desktop/Undergraduate_Thesis/pairwise-Hamming-distance-profile/test/replications_"

# Pattern of replicate file names
REPLICATE_PATTERN = "rep_{}.txt"

# Number of replicate experiments
NUM_REPLICATES = 1000

# Data features
SAMPLE_RATES = [0.999 , 0.5 , 0.25 , 0.1 , 0.01]
HD_MAX = 32
COLORS = "bwr"

# Output directory for figures
OUTPUT_DIR = "/Users/kesil/Desktop/Undergraduate_Thesis/pairwise-Hamming-distance-profile/test/errorFigures"

def read_hd_counts(file_path):
    """
    Read Hamming-distance counts from a text file.

    Expected format:
        Hamming Distance : Number of Pairs
        0 : 0
        1 : 10
        2 : 100
        ...

    Returns
    -------
    dict
        Dictionary mapping Hamming distance to number of pairs.
    """
    counts = {}

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()

            # Match lines such as:
            # 23 : 6981408
            match = re.match(r"^(\d+)\s*:\s*([0-9]+.?[0-9]*e?\+?[0-9]*)", line)

            if match:
                hd = int(match.group(1))
                count = float(match.group(2))
                counts[hd] = count
                #if(hd == 24):
                #    print(hd , count)

    return counts

# ============================================================
# Load ground truth
# ============================================================

ground_truth = read_hd_counts(GROUND_TRUTH_FILE)

if len(ground_truth) == 0:
    raise ValueError(
        f"No Hamming-distance counts were found in {GROUND_TRUTH_FILE}"
    )

hamming_distances = sorted(ground_truth.keys())

#hamming_distances = [hd for hd in hamming_distances if hd > 11]

print(f"Loaded ground truth for {len(hamming_distances)} Hamming distances.")
print(
    f"Hamming-distance range: "
    f"{min(hamming_distances)} - {max(hamming_distances)}"
)

# ============================================================
# Load replicate results
# ============================================================

# Must load results for each sampling rate
all_est_counts = {}
for rate in SAMPLE_RATES:
    # Shape:
    #   rows    = replicate experiments
    #   columns = Hamming distances
    estimated_counts = np.full(
        (NUM_REPLICATES, len(hamming_distances)),
        np.nan,
        dtype=float
    )

    for rep in range(0, NUM_REPLICATES):

        file_name = REPLICATE_PATTERN.format(rep)
        file_path = os.path.join((REPLICATE_DIR+str(rate)), file_name)

        if not os.path.exists(file_path):
            print(f"Warning: file not found: {file_path}")
            continue

        counts = read_hd_counts(file_path)

        for j, hd in enumerate(hamming_distances):
            if hd in counts:
                estimated_counts[rep - 1, j] = counts[hd]

    all_est_counts[rate] = estimated_counts

    print(
        f"Successfully loaded "
        f"{np.sum(~np.all(np.isnan(estimated_counts), axis=1))} "
        f"replicate files for rate = {rate}."
    )

# ============================================================
# Calculate absolute errors
# ============================================================

all_errors = {}
all_rel_errors = {}
for rate in all_est_counts:
    true_counts = np.array(
        [ground_truth[hd] for hd in hamming_distances],
        dtype=float
    )

    # Error:
    #
    #     error = p_hat - p
    #
    errors = abs(all_est_counts[rate] - true_counts[np.newaxis, :])
    all_errors[rate] = np.copy(errors)

    # Relative error:
    #
    #     relative_error = (p_hat - p) / p
    #
    # Relative error is undefined when p = 0.
    relative_errors = np.full_like(errors, np.nan)

    nonzero_mask = true_counts != 0

    relative_errors[:, nonzero_mask] = (
        errors[:, nonzero_mask]
        / true_counts[np.newaxis, nonzero_mask]
    )
    all_rel_errors[rate] = np.copy(relative_errors)

#print(all_rel_errors[0.999][-1])
#print(all_rel_errors[0.5][-1])

# ============================================================
# Calculate average absolute relative errors
# ============================================================

all_avgs = {}
for rate in all_rel_errors:
    rel_errors = all_rel_errors[rate]

    # Want average of all relative errors for the same HD
    sums = [0 for i in range(HD_MAX)]
    avgs = [np.nan for i in range(HD_MAX)]
    divs = [0 for i in range(HD_MAX)]
    #rows = reps
    #cols = HDs
    for col in range(HD_MAX):
        for row in range(NUM_REPLICATES):
            if not(math.isnan(rel_errors[row][col])):
                sums[col] += rel_errors[row][col]
                divs[col] += 1
    for i in range(HD_MAX):
        if(divs[i] != 0):
            avgs[i] = sums[i]/divs[i]
    all_avgs[rate] = avgs

# ============================================================
# Create output directory
# ============================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# Heat map of average relative error
# ============================================================

avgs = np.array([all_avgs[avg][7:] for avg in all_avgs] , dtype=float)

fig, ax = plt.subplots()
im = ax.imshow(avgs , cmap=COLORS)

#X = HD
ax.set_xticks(range(HD_MAX-7) , labels=[str(i+1) for i in range(7,HD_MAX)] , rotation=45)
#Y = Sampling rate
ax.set_yticks(range(len(SAMPLE_RATES)) , labels=[(str(i) + "*" + str(i)) for i in SAMPLE_RATES])

#for i in range(len(SAMPLE_RATES)):
#    for j in range(HD_MAX):
#        text = ax.text(j , i , avgs[i , j] , ha="center" , va="center" , color="w")

ax.set_title("Average Relative Error Based on Sampling Rate")
ax.set_xlabel("Hamming Distance")
ax.set_ylabel("Sampling Rate")
fig.tight_layout()

im.set_clim(0.0 , 0.2)
cbar = ax.figure.colorbar(im , ax=ax)
cbar.ax.set_ylabel("" , rotation=-90 , va="bottom")
#ax.spines[:].set_visible(False)

output = os.path.join(
    OUTPUT_DIR,
    "heatmap.png"
)

fig.savefig(output, dpi=300)
plt.close()

print(f"Saved error plot to: {output}")
