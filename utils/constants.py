from enum import Enum


class Layout(Enum):
    PE = "paired-end"
    SE = "single-end"
    SLC = "single-cell"


class PrepMethod(Enum):
    total = "total"
    mrna = "mrna"
    polya = "mrna"  # Write mrna if polya OR mrna
