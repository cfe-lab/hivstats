import matplotlib.pyplot as plt
import numpy as np
import csv
from dataclasses import dataclass
from itertools import zip_longest
import statistics
import Levenshtein
from Bio import Align
from scipy.stats import gaussian_kde

from mynotebook_data import get_joined


interactive_mode = None


def show_graphics():
    if interactive_mode is True:
        plt.show()


def compute_kde_bandwidth(scores):
    """Compute bandwidth using Silverman's rule of thumb.

    h ≈ 0.9 * min(σ, IQR/1.34) * n^(-1/5)

    Args:
        scores: List of numeric values

    Returns:
        Bandwidth value for KDE
    """
    if len(scores) < 2:
        return 1.0

    scores_array = np.array(scores)
    sigma = np.std(scores_array, ddof=1)

    # Handle zero or invalid standard deviation
    if sigma == 0 or np.isnan(sigma) or np.isinf(sigma):
        return 1.0

    q75, q25 = np.percentile(scores_array, [75, 25])
    iqr = q75 - q25

    # Handle zero IQR (all values in middle 50% are identical)
    if iqr == 0:
        scale = sigma
    else:
        scale = min(sigma, iqr / 1.34)

    # Silverman's rule of thumb
    h = 0.9 * scale * (len(scores) ** (-1 / 5))

    # Ensure bandwidth is reasonable (not too small, not invalid)
    if h <= 0 or np.isnan(h) or np.isinf(h):
        return max(sigma * 0.1, 0.01)

    return max(h, 0.01)


def compute_kde(scores, num_points=1000, bw_method=None):
    """Compute KDE for a set of scores.

    Args:
        scores: List of numeric values
        num_points: Number of points to evaluate KDE at
        bw_method: Bandwidth method ('silverman', 'scott', or numeric value)
                   If None, uses Silverman's rule

    Returns:
        (x_values, density_values) tuple for plotting
    """
    if len(scores) < 2:
        return None, None

    scores_array = np.array(scores)

    # Check for zero variance (all values identical or nearly identical)
    std = scores_array.std(ddof=1)
    if std == 0 or np.isnan(std) or np.isinf(std):
        return None, None

    # Check if range is too small for meaningful KDE
    x_min, x_max = scores_array.min(), scores_array.max()
    x_range = x_max - x_min
    if x_range < 1e-10:  # Values are essentially identical
        return None, None

    # Determine bandwidth
    if bw_method is None:
        bw = compute_kde_bandwidth(scores)
        # Avoid division by zero or near-zero std
        if std > 1e-10:
            bw_method = bw / std
        else:
            return None, None

    # Create KDE
    try:
        kde = gaussian_kde(scores_array, bw_method=bw_method)
    except Exception:
        return None, None

    # Evaluate KDE
    try:
        x_vals = np.linspace(x_min - 0.1 * x_range, x_max + 0.1 * x_range, num_points)
        density = kde(x_vals)

        # Check for invalid values
        if np.any(np.isnan(density)) or np.any(np.isinf(density)):
            return None, None

    except Exception:
        return None, None

    return x_vals, density


def scale_kde_to_histogram(density, x_vals, hist_counts, bin_width):
    """Scale KDE density to match histogram count scale.

    Args:
        density: KDE density values
        x_vals: x-values where KDE was evaluated
        hist_counts: Maximum count in histogram
        bin_width: Width of histogram bins

    Returns:
        Scaled density values
    """
    if density is None or len(density) == 0:
        return density

    # Scale so KDE peak roughly matches histogram scale
    # Histogram shows counts, KDE shows density
    # To match: multiply density by (n * bin_width)
    kde_max = np.max(density)
    if kde_max > 0:
        scale_factor = hist_counts * bin_width / kde_max
        return density * scale_factor
    return density


def print_statistics(name, scores):
    # Calculate statistics
    mean = statistics.mean(scores) if scores else 0
    median = statistics.median(scores) if scores else 0
    try:
        mode = statistics.mode(scores)
    except BaseException:
        mode = None
    stdev = statistics.stdev(scores) if len(scores) > 1 else 0
    # variance = statistics.variance(scores) if len(scores) > 1 else 0
    min_score = min(scores) if scores else float("inf")
    max_score = max(scores) if scores else 0

    # Print statistics
    print(f"Name: {name}")
    print(f"Count: {len(scores)}")
    print(f"Mean: {round(mean, 2)}")
    print(f"Median: {round(median, 2)}")
    print(f"Mode: {round(mode, 2) if mode is not None else 'undefined'}")
    print(f"Standard Deviation: {round(stdev, 2)}")
    print(f"Minimum: {min_score}")
    print(f"Maximum: {max_score}")


def levenshtein_distance(s1, s2):
    return Levenshtein.distance(s1, s2)


aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -1.5
aligner.extend_gap_score = -0.2


def aligner_distance(query, reference):
    if len(query) == 0:
        return float("inf")

    orf_alignments = aligner.align(query, reference)
    if not orf_alignments:
        return float("inf")

    orf_alignment = orf_alignments[0]
    return (-1 * orf_alignment.score) / len(query) + aligner.match_score


ORFs = [
    "gag",
    "pol",
    "env",
    "vif",
    "vpr",
    "tat_exon1",
    "rev_exon1",
    "vpu",
    "tat_exon2",
    "rev_exon2",
    "nef",
]


@dataclass
class Data:
    name: str
    scores: list
    start: float
    end: float


def unranged(name, generator):
    return Data(name, list(generator), None, None)


def ranged(name, start, end, generator):
    return Data(name, list(generator), start, end)


def get_distance_scores_nw(orf, joined):
    for val in joined:
        if val["region"] == orf:
            yield val["distance"]


def get_size_scores_nw(orf, joined):
    for val in joined:
        if val["region"] == orf:
            aminos = val["protein"]
            yield len(aminos) * 3


def get_nostopcodon_size_scores_nw(orf, joined):
    for val in joined:
        if val["region"] == orf:
            aminos = val["aminoacids"]
            yield len(aminos)


def get_asize_scores_nw(orf, joined):
    for val in joined:
        if val["region"] == orf:
            astart = val["start"]
            aend = val["end"]
            yield aend - astart + 1


def get_indel_scores_nw(orf, joined):
    for val in joined:
        if val["region"] == orf:
            yield float(val["indel_impact"])


def get_scores_all(orf, metric, outliers, joined):
    # print(f"joined: {joined[:10]}")

    if metric == "distance":
        return ranged("Distance", 0, 2, get_distance_scores_nw(orf, joined))
    elif metric == "size (protein)":
        return unranged("Size", get_size_scores_nw(orf, joined))
    elif metric == "size":
        return unranged("Size", get_asize_scores_nw(orf, joined))
    elif metric == "indel impact":
        return unranged("indel_impact", get_indel_scores_nw(orf, joined))
    else:
        raise ValueError(f"Invalid choice of metric: {metric}")


def get_scores(orf, metric, outliers, joined):
    data = get_scores_all(orf, metric, outliers, joined)

    # outliers = 0.01 # percentage of outliers
    lower_threshold = np.quantile(data.scores, outliers) if data.scores else 0

    upper_threshold = np.quantile(data.scores, 1 - outliers) if data.scores else 0
    # lower_threshold = np.percentile(data.scores, outliers * 100)
    # upper_threshold = np.percentile(data.scores, 100 - (outliers * 100))

    # print(f'orf  : {orf}')
    # print(f'lower: {lower_threshold}')
    # print(f'upper: {upper_threshold}')

    data.scores = [
        x for x in data.scores if x >= lower_threshold and x <= upper_threshold
    ]
    return data


def filter_based_on_intactness(goodq, joined, metric):
    """Filter sequences based on appropriate intactness criteria for the metric.

    Uses progressive intactness filtering to avoid circular dependencies:
    - For 'size' metric: use structural intactness only
    - For 'distance' metric: use structural + size intactness
    - For 'indel impact' metric: use structural + size + distance intactness

    Args:
        goodq: If True, return intact sequences; if False, return defective
        joined: List of sequence dictionaries
        metric: Metric being analyzed ('size', 'distance', 'indel impact', etc.)

    Yields:
        Sequence dictionaries matching the intactness criteria
    """
    # Select appropriate intactness field based on metric
    if metric == "size" or metric == "size (protein)":
        intact_field = "size_structural_intact"
    elif metric == "distance":
        intact_field = "distance_intact"
    elif metric == "indel impact":
        intact_field = "indel_intact"
    else:
        # Default to most basic check for unknown metrics
        intact_field = "size_structural_intact"

    for val in joined:
        if val[intact_field] == goodq:
            yield val


def show_it(orf, data):
    scores = data.scores
    numrange = [data.start, data.end] if data.start is not None else None

    # Create histogram
    counts, bins, patches = plt.hist(
        scores, bins=30, edgecolor="black", range=numrange, align="mid", alpha=0.7
    )

    # Compute and plot KDE
    x_vals, density = compute_kde(scores)
    if x_vals is not None and density is not None:
        # Scale KDE to match histogram
        bin_width = bins[1] - bins[0] if len(bins) > 1 else 1
        scaled_density = scale_kde_to_histogram(
            density, x_vals, np.max(counts), bin_width
        )
        plt.plot(x_vals, scaled_density, "r-", linewidth=2, label="KDE", alpha=0.8)
        plt.legend()

    plt.xlabel(data.name)
    plt.ylabel("Count")
    plt.title(f"Distribution of {orf}")
    show_graphics()


def show_two(orf, data_good, data_bad):
    scores_good = data_good.scores
    scores_bad = data_bad.scores
    numrange = [data_good.start, data_good.end] if data_good.start is not None else None

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax1.set_ylabel("Intact count")
    counts_good, bins_good, _ = ax1.hist(
        scores_good,
        bins=30,
        alpha=0.5,
        label="Intact",
        edgecolor="black",
        color="black",
        range=numrange,
        align="mid",
    )

    # Compute and plot KDE for intact
    x_vals_good, density_good = compute_kde(scores_good)
    if x_vals_good is not None and density_good is not None:
        bin_width = bins_good[1] - bins_good[0] if len(bins_good) > 1 else 1
        scaled_density_good = scale_kde_to_histogram(
            density_good, x_vals_good, np.max(counts_good), bin_width
        )
        ax1.plot(
            x_vals_good,
            scaled_density_good,
            "darkblue",
            linewidth=2,
            label="Intact KDE",
            alpha=0.8,
            linestyle="--",
        )

    ax2.set_ylabel("Defective count")
    counts_bad, bins_bad, _ = ax2.hist(
        scores_bad,
        bins=30,
        alpha=0.5,
        label="Defective",
        edgecolor="black",
        color="red",
        range=numrange,
        align="mid",
    )

    # Compute and plot KDE for defective
    x_vals_bad, density_bad = compute_kde(scores_bad)
    if x_vals_bad is not None and density_bad is not None:
        bin_width = bins_bad[1] - bins_bad[0] if len(bins_bad) > 1 else 1
        scaled_density_bad = scale_kde_to_histogram(
            density_bad, x_vals_bad, np.max(counts_bad), bin_width
        )
        ax2.plot(
            x_vals_bad,
            scaled_density_bad,
            "darkred",
            linewidth=2,
            label="Defective KDE",
            alpha=0.8,
            linestyle="--",
        )

    plt.xlabel(data_good.name)
    plt.title(f"Distribution of {orf}")

    # Combine legends from both axes
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="best")

    fig.tight_layout()
    show_graphics()


def process_orf(joined, orf, metric, outliers):
    scores = get_scores(orf, metric, outliers, joined)

    # print(f"scores: {scores.scores[:10]}")

    show_it(orf, scores)
    # print(pd.DataFrame(np.array(scores.scores, dtype=float)).describe())
    print_statistics(orf, scores.scores)
    print("------------------------------------------")


def process_two_orfs(joined, orf, metric, outliers):
    good = filter_based_on_intactness(True, joined, metric)
    bad = filter_based_on_intactness(False, joined, metric)
    scores_good = get_scores(orf, metric, outliers, good)
    scores_bad = get_scores(orf, metric, outliers, bad)
    show_two(orf, scores_good, scores_bad)

    print("Intact:")
    print_statistics(orf, scores_good.scores)
    print("")
    print("Nonintact:")
    print_statistics(orf, scores_bad.scores)
    print("------------------------------------------")


def dump_data(joined):
    with open("output/data.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerow(ORFs)
        columns = [get_scores(orf, joined).scores for orf in ORFs]
        for row in zip_longest(*columns, fillvalue=""):
            writer.writerow(row)


def show_size_examples():
    if interactive_mode is not True:
        return

    import ipywidgets as widgets
    from ipywidgets import interactive

    def fx(reference, query):
        print(f"distance(a, b) = {aligner_distance(reference, query)}")
        print(f"levenshtein(a, b) = {levenshtein_distance(reference, query)}")

    acc = []

    acc.append(widgets.Label(value="Example 1:"))

    acc.append(interactive(fx, reference="AAAAAAA", query="BBBBBBB"))

    acc.append(widgets.Label(value="Example 2:"))

    acc.append(interactive(fx, reference="A", query="BBBBBBB"))

    acc.append(widgets.Label(value="Example 3:"))

    acc.append(interactive(fx, reference="VITALIY", query="VITALIK"))

    header = widgets.VBox(
        [
            widgets.Label("Distance is calculated between sequences of aminoacids."),
            widgets.Label(
                "It is based on how well those sequences align.\n\n"
                "In below examples levenshtein distance is used for comparison, it's not what we use in CFEIntact becase of its performance."
            ),
        ]
    )

    mainpart = widgets.HBox(acc)

    w = widgets.Accordion([widgets.VBox([header, mainpart])])
    w.set_title(0, "What is distance?")

    # display(w)


def show_all_orfs(joined, extractedby, select, metric, outliers):
    if metric == "distance":
        show_size_examples()

    for orf in ORFs:
        if select == "together":
            process_orf(joined, orf, metric, outliers)
        elif select == "intact":
            process_orf(
                filter_based_on_intactness(True, joined, metric), orf, metric, outliers
            )
        elif select == "nonintact":
            process_orf(
                filter_based_on_intactness(False, joined, metric), orf, metric, outliers
            )
        elif select == "separately":
            process_two_orfs(joined, orf, metric, outliers)
        else:
            raise ValueError(f"Invalid choice of select: {select}")


def jupyter_main():
    import ipywidgets as widgets
    from ipywidgets import interact

    global interactive_mode
    interactive_mode = True

    def interactable(extractedby, metric, outliers):
        if extractedby == "Los Alamos/Plasma":
            joined = list(get_joined(source="los-alamos/plasma"))
        elif extractedby == "CFEIntact/All":
            joined = list(get_joined(source="cfeintact/all"))
        elif extractedby == "CFEIntact/Plasma":
            joined = list(get_joined(source="cfeintact/plasma"))
        else:
            raise ValueError(f"Unexpected extractedby: {extractedby!r}.")

        def cont(select):
            show_all_orfs(joined, extractedby, select, metric, outliers)

        interact(
            cont,
            select=widgets.ToggleButtons(
                options=["intact", "nonintact", "together", "separately"],
                description="Intact vs Defective (CFEIntact opinion):"
                if "CFEIntact" in extractedby
                else "Intact vs Defective (Based on stop codon presense)",
            ),
        )

    slider_outliers = widgets.FloatSlider(
        value=0.01,
        min=0,
        max=0.1,
        step=0.01,
        description="Outliers:",
        disabled=False,
        continuous_update=False,
        readout_format="2.0%",
        orientation="horizontal",
        readout=True,
    )

    interact(
        interactable,
        extractedby=widgets.ToggleButtons(
            options=["CFEIntact/Plasma", "Los Alamos/Plasma", "CFEIntact/All"],
            description="Extracted by:",
        ),
        metric=widgets.ToggleButtons(
            options=["size", "size (protein)", "distance", "indel impact"],
            description="Metric:",
        ),
        outliers=slider_outliers,
    )
