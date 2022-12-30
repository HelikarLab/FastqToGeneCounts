"""
This file is responsible for validating that each CELL-TYPE_S##R## has at least two runs per study
"""
import csv
from pathlib import Path
import re

def validate(config_file: str | Path) -> bool:
    schema = _ingest(config_file)
    return _is_valid(schema)

def _ingest(control: str | Path) -> dict[str, dict[str, list[str]]]:
    """
    This function is responsible for taking data from the control file and creating a usable data format
    :param control: The control file to read from
    :return: A dictionary containing each cell type, their studies, and each studies' replicates
    """
    schema: dict[str, dict[str, list[str]]] = {}
    
    with open(control, "r") as i_stream:
        reader = csv.reader(i_stream)
        for line in reader:
            sample = line[1]                                # "naiveB_S1R1"
            cell_type = sample.split("_")[0]                # naiveB
            tag = sample.split("_")[1]                      # S1R1
            study = re.search(r"S\d+", tag).group()         # S1
            run = re.search(r"R\d+(r\d+)?", tag).group()    # R1 (with optional replicate, i.e., R1r1)
            
            if cell_type not in schema.keys():
                schema[cell_type] = {}
            
            cell_studies = schema[cell_type]
            if study not in cell_studies.keys():
                cell_studies[study] = []
            
            study_list: list[str] = cell_studies[study]
            if run not in study_list:
                study_list.append(run)
            
    return schema
    
def _is_valid(schema: dict[str, dict[str, list[str]]]) -> bool | None:
    """
    This function is responsible for validating that each study in the schema contains at least two replicates
    :param schema: The schema generated from the ingest() function
    :return:
    """
    for cell_type in schema.keys():  # Iterate through cell types (naiveB, naivecd8)
        
        studies: dict[str, list[str]] = schema[cell_type]
        for study in studies.keys():  # Iterate through studies (S1, S2)
        
            replicates: list[str] = studies[study]  # Get replicates in study (R1, R2)
            if len(replicates) < 2:
                print(f"\n\nYou have included a study that contains a single replicate.\n"
                      f"If you are going to continue using this data with MADRID, it will result in errors during the first stage of the analysis process.\n"
                      f"Please exclude this study, or include additional replicates\n\n")
                
                raise ValueError(f"Error in `{cell_type}_{study}{replicates[0]}`. Only 1 replicate found")
            
    return True
