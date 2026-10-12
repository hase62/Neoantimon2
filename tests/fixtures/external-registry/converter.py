"""Authored synthetic test converter, equivalent to the accepted R example."""
import argparse
import csv

parser = argparse.ArgumentParser()
for arg in ("--input", "--output", "--bundle"):
    parser.add_argument(arg, required=True)
args = parser.parse_args()
fields = ["schema_version", "peptide_sequence", "scope", "mhc_key", "task", "provider",
          "predictor_version", "model_version", "run_id", "metric", "metric_kind",
          "metric_semantics", "score", "score_direction", "unit", "status", "reason", "source_row"]
with open(args.input, newline="") as source, open(args.output, "w", newline="") as destination:
    writer = csv.DictWriter(destination, fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for i, row in enumerate(csv.DictReader(source, delimiter="\t"), 1):
        writer.writerow(dict(schema_version="1.0.0", peptide_sequence=row["peptide"], scope="peptide_hla",
                             mhc_key=row["hla"], task="external:synthetic_score", provider="synthetic-model",
                             predictor_version="1", model_version="demo-1", run_id="demo-run-1", metric="score",
                             metric_kind="raw_score", metric_semantics="raw", score=row["value"],
                             score_direction="higher_is_better", unit="unitless", status="evaluated", reason="",
                             source_row=f"X:{i}"))
