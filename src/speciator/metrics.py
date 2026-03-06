import csv
import dataclasses
from collections import defaultdict
from pathlib import Path
from typing import Iterable


@dataclasses.dataclass
class Metric:
    species: str
    name: str
    lower_bounds: float
    upper_bounds: float


def read_metrics(file: Path) -> dict[str, dict[str, Metric]]:
    all_metrics: dict[str, dict[str, Metric]] = defaultdict(dict)
    with open(file, "r") as metrics_fh:
        reader: Iterable[dict[str, str]] = csv.DictReader(metrics_fh, delimiter=",")
        for line in reader:
            all_metrics[line["species"]][line["metric"]] = Metric(
                line["species"],
                line["metric"],
                float(line["lower_bounds"] if line["lower_bounds"] else 0.0),
                float(line["upper_bounds"] if line["upper_bounds"] else float("inf")),
            )
    return all_metrics
