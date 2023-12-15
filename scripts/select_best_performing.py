import json
import pathlib

# Load performance scores of all the previous rounds
scores = {}
for num, score in enumerate(snakemake.input.scores):
    with open(score) as scorefile:
        scores[num] = json.load(scorefile)['overall']

sorted_scores = sorted(scores.keys(), key = lambda num: scores[num])
best_scorer = pathlib.Path(snakemake.input.scores[sorted_scores[0]]).parent
print(f"Found the best performing iteration: {best_scorer}")


with open(snakemake.output.config, "wt") as f:
    f.write((best_scorer / "beers.config.yaml").read_text())

out_dir = pathlib.Path(snakemake.output.out_dir)
out_dir.mkdir(exist_ok=True)

for filename in ["coverage.txt", "coverage_summary.txt", "gc_content.txt", "seq_frequencies.json", "frag_sizes.txt"]:
    with open(best_scorer / "sample1" / filename, "rt") as f:
        contents = f.read()
    with open(out_dir / filename, "wt") as f:
        f.write(contents)
