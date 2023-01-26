from enum import Enum

class EndType(Enum):
    PE = "paired_end"
    SE = "single_end"
    SLC = "single_cell"
    
class PrepMethod(Enum):
    total = "total"
    mrna = "mrna"
    polya = "mrna"  # Write mrna if polya OR mrna
