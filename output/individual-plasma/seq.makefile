output/individual-plasma/joined.csv: src/join-csv-files output/individual-plasma/seq/env.csv output/individual-plasma/seq/gag.csv output/individual-plasma/seq/nef.csv output/individual-plasma/seq/pol.csv output/individual-plasma/seq/rev_exon1.csv output/individual-plasma/seq/rev_exon2.csv output/individual-plasma/seq/tat_exon1.csv output/individual-plasma/seq/tat_exon2.csv output/individual-plasma/seq/vif.csv output/individual-plasma/seq/vpr.csv output/individual-plasma/seq/vpu.csv
	@ mkdir -p output/individual-plasma
	uv run python -- $^ $@

output/individual-plasma/seq/env.csv: src/make-individual-plasma-csv input/individual-plasma/seq/env.fasta
	@ mkdir -p output/individual-plasma/seq/
	uv run python -- $^ $@

output/individual-plasma/seq/gag.csv: src/make-individual-plasma-csv input/individual-plasma/seq/gag.fasta
	@ mkdir -p output/individual-plasma/seq/
	uv run python -- $^ $@

output/individual-plasma/seq/nef.csv: src/make-individual-plasma-csv input/individual-plasma/seq/nef.fasta
	@ mkdir -p output/individual-plasma/seq/
	uv run python -- $^ $@

output/individual-plasma/seq/pol.csv: src/make-individual-plasma-csv input/individual-plasma/seq/pol.fasta
	@ mkdir -p output/individual-plasma/seq/
	uv run python -- $^ $@

output/individual-plasma/seq/rev_exon1.csv: src/make-individual-plasma-csv input/individual-plasma/seq/rev_exon1.fasta
	@ mkdir -p output/individual-plasma/seq/
	uv run python -- $^ $@

output/individual-plasma/seq/rev_exon2.csv: src/make-individual-plasma-csv input/individual-plasma/seq/rev_exon2.fasta
	@ mkdir -p output/individual-plasma/seq/
	uv run python -- $^ $@

output/individual-plasma/seq/tat_exon1.csv: src/make-individual-plasma-csv input/individual-plasma/seq/tat_exon1.fasta
	@ mkdir -p output/individual-plasma/seq/
	uv run python -- $^ $@

output/individual-plasma/seq/tat_exon2.csv: src/make-individual-plasma-csv input/individual-plasma/seq/tat_exon2.fasta
	@ mkdir -p output/individual-plasma/seq/
	uv run python -- $^ $@

output/individual-plasma/seq/vif.csv: src/make-individual-plasma-csv input/individual-plasma/seq/vif.fasta
	@ mkdir -p output/individual-plasma/seq/
	uv run python -- $^ $@

output/individual-plasma/seq/vpr.csv: src/make-individual-plasma-csv input/individual-plasma/seq/vpr.fasta
	@ mkdir -p output/individual-plasma/seq/
	uv run python -- $^ $@

output/individual-plasma/seq/vpu.csv: src/make-individual-plasma-csv input/individual-plasma/seq/vpu.fasta
	@ mkdir -p output/individual-plasma/seq/
	uv run python -- $^ $@

