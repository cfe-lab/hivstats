
all: output/results.txt

TMP_RESULTS = ./output/temporary_results.txt

output/results.txt: output/fullgenomes-all/regions.csv output/fullgenomes-plasma/regions.csv output/individual-plasma/joined.csv src/print_results.py src/mynotebook.py src/print_results.py src/mynotebook_data.py
	uv run -- python src/print_results.py 1>$(TMP_RESULTS)
	mv -- $(TMP_RESULTS) "$@"

serve:
	uv run -- jupyter notebook src/main.ipynb

reanalyze:
	rm -rf output/
	$(MAKE) output/results.txt
	cat output/results.txt

output/fullgenomes-plasma/regions.csv: src/run-cfeintact output/fullgenomes-plasma.fasta
	$^

output/fullgenomes-plasma.fasta: input/fullgenomes-plasma/los-alamos-plasma-sequences.fasta input/fullgenomes-plasma/other-intact-sequences.fasta
	@ mkdir -p output
	cat $^ > "$@"

output/fullgenomes-all/regions.csv: src/run-cfeintact output/fullgenomes-all.fasta
	$^

output/fullgenomes-all.fasta: input/fullgenomes-all/los-alamos-all-sequences.fasta
	@ mkdir -p output
	cat $^ > "$@"


ifeq ($(wildcard output/individual-plasma/seq.makefile),)
output/individual-plasma/joined.csv: output/individual-plasma/seq.makefile
	$(MAKE) -f output/individual-plasma/seq.makefile
else
include output/individual-plasma/seq.makefile
endif

output/individual-plasma/seq.makefile: src/generate-individual-plasma-makefile
	@ mkdir -p output/individual-plasma/
	$^ > "$@"

clean:
	rm -rf output

.PHONY: all csvs serve clean
.SECONDARY:
