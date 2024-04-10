
import matplotlib.pyplot as plt
import numpy as np
import csv
from dataclasses import dataclass
from itertools import zip_longest
import statistics
import Levenshtein
from Bio import Align

from mynotebook_data import get_joined


interactive_mode = None


def show_graphics():
    if interactive_mode is True:
        plt.show()


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
aligner.mode = 'global'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -1.5
aligner.extend_gap_score = -0.2

def aligner_distance(query, reference):
   if len(query) == 0:
       return float('inf')

   orf_alignments = aligner.align(query, reference)
   if not orf_alignments:
       return float('inf')

   orf_alignment = orf_alignments[0]
   return (-1 * orf_alignment.score) / len(query) + aligner.match_score


ORFs = ["gag", "pol", "env", "vif", "vpr", "tat_exon1", "rev_exon1", "vpu", "tat_exon2", "rev_exon2", "nef"]


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

    data.scores = [x for x in data.scores if x >= lower_threshold and x <= upper_threshold]
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

    plt.hist(scores, bins=30, edgecolor='black', range=numrange, align='mid')  # Adjust the number of bins as per your preference
    # Draw a horizontal line for the range
    # plt.axhline(y=400, xmin=0.1, xmax=0.5, linestyle='--')

    plt.xlabel(data.name)
    plt.ylabel('Count')
    plt.title(f'Distribution of {orf}')
    # plt.legend()
    show_graphics()


def show_two(orf, data_good, data_bad):
    scores_good = data_good.scores
    scores_bad = data_bad.scores
    numrange = [data_good.start, data_good.end] if data_good.start is not None else None

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax1.set_ylabel('Intact count')
    ax1.hist(scores_good, bins=30, alpha=0.5, label='Intact', edgecolor='black', color='black', range=numrange, align='mid')

    ax2.set_ylabel('Defective count')
    ax2.hist(scores_bad, bins=30, alpha=0.5, label='Defective', edgecolor='black', color='red', range=numrange, align='mid')

    plt.xlabel(data_good.name)
    plt.title(f'Distribution of {orf}')
    plt.legend(loc = 0)
    # plt.legend()
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
        for row in zip_longest(*columns, fillvalue=''):
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

    header = widgets.VBox([widgets.Label("Distance is calculated between sequences of aminoacids."),
                           widgets.Label("It is based on how well those sequences align.\n\n"
                                         "In below examples levenshtein distance is used for comparison, it's not what we use in CFEIntact becase of its performance.")])

    mainpart = widgets.HBox(acc)

    w = widgets.Accordion([widgets.VBox([header, mainpart])])
    w.set_title(0, "What is distance?")

    # display(w)


def show_all_orfs(joined, extractedby, select, metric, outliers):
    if metric == "distance":
        show_size_examples()

    for orf in ORFs:
        if select == 'together':
            process_orf(joined, orf, metric, outliers)
        elif select == 'intact':
            process_orf(filter_based_on_intactness(True, joined, metric), orf, metric, outliers)
        elif select == 'nonintact':
            process_orf(filter_based_on_intactness(False, joined, metric), orf, metric, outliers)
        elif select == 'separately':
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
            joined = list(get_joined(source='los-alamos/plasma'))
        elif extractedby == "CFEIntact/All":
            joined = list(get_joined(source='cfeintact/all'))
        elif extractedby == "CFEIntact/Plasma":
            joined = list(get_joined(source='cfeintact/plasma'))
        else:
            raise ValueError(f"Unexpected extractedby: {extractedby!r}.")

        def cont(select):
            show_all_orfs(joined, extractedby, select, metric, outliers)

        interact(cont,
                 select=widgets.ToggleButtons(options=["intact", "nonintact", "together", "separately"],
                                              description='Intact vs Defective (CFEIntact opinion):' if "CFEIntact" in extractedby \
                                              else "Intact vs Defective (Based on stop codon presense)"))

    slider_outliers=widgets.FloatSlider(
        value=0.01,
        min=0,
        max=0.1,
        step=0.01,
        description='Outliers:',
        disabled=False,
        continuous_update=False,
        readout_format='2.0%',
        orientation='horizontal',
        readout=True,
    )

    interact(interactable,
             extractedby=widgets.ToggleButtons(options=["CFEIntact/Plasma", "Los Alamos/Plasma", "CFEIntact/All"],
                                               description='Extracted by:'),
             metric=widgets.ToggleButtons(options=["size", "size (protein)", "distance", "indel impact"],
                                          description='Metric:'),
             outliers=slider_outliers)

