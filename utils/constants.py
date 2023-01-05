from enum import Enum

class EndType(Enum):
    paired_end = "PE"
    single_end = "SE"
    single_cell = "SLC"
    
class PrepMethod(Enum):
    total = "total"
    mrna = "mrna"
