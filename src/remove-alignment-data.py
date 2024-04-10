#! /usr/bin/env python3

import argparse
import sys


def remove_dashes(input_seq):
    return input_seq.replace("-", "")


def format_fasta_output(header, sequence, line_width):
    formatted_sequence = "\n".join([sequence[i:i+line_width] for i in range(0, len(sequence), line_width)])
    return f">{header}\n{formatted_sequence}"


def infer_line_width(first_line):
    # Infer the line width from the length of the first line in the sequence
    return len(first_line)


def process_fasta_file(input_file, output_file, line_width):
    with open(input_file, "r") as f:
        lines = f.readlines()

    sequences = {}
    line_widths = {}
    current_header = ""
    current_sequence = ""

    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if current_header:
                sequences[current_header] = current_sequence
                current_sequence = ""
            current_header = line[1:]
        else:
            if current_header not in line_widths:
                line_widths[current_header] = len(line)
            current_sequence += line

    if current_header and current_sequence:
        sequences[current_header] = current_sequence

    with open(output_file, "w") as f:
        for header, sequence in sequences.items():
            formatted_sequence = remove_dashes(sequence)
            if line_width is None:
                line_width_inferred = line_widths[header]
            else:
                line_width_inferred = line_width
            formatted_output = format_fasta_output(header, formatted_sequence, line_width_inferred)
            f.write(formatted_output + "\n")


def main(argv) -> int:
    parser = argparse.ArgumentParser(description="Process a FASTA file by removing '-' characters from sequences while maintaining line widths.")
    parser.add_argument("input_file", help="Path to the input FASTA file")
    parser.add_argument("output_file", help="Path to the output FASTA file")
    parser.add_argument("--line_width", type=int, default=None, help="Desired line width for the output sequences")

    args = parser.parse_args(argv)
    process_fasta_file(args.input_file, args.output_file, args.line_width)
    print("FASTA processing complete.")
    return 0


if __name__ == "__main__":
    exit(main(sys.argv[1:]))
