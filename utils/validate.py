"""
This file is responsible for validating that each CELL-TYPE_S##R## has at least two runs per study
"""
import csv
from pathlib import Path
import re
from typing import Union


class _GenomeData:
    @classmethod
    def validate(cls, config_filename, config):
        # Validate the genome input paths
        reference_flat_file: Path = Path(config["REF_FLAT_FILE"])
        rRna_interval_list: Path = Path(config["RRNA_INTERVAL_LIST"])
        bed_file: Path = Path(config["BED_FILE"])
        genome_fasta_file: Path = Path(config["GENERATE_GENOME"]["GENOME_FASTA_FILE"])
        gtf_file: Path = Path(config["GENERATE_GENOME"]["GTF_FILE"])
        
        should_perform_rnaseq_metrics = config["PERFORM_GET_RNASEQ_METRICS"]
        # should_perform_insert_size = config["PERFORM_GET_INSERT_SIZE"]
        should_perform_get_fragment_size = config["PERFORM_GET_FRAGMENT_SIZE"]

        genome_valid = True
        if not reference_flat_file.exists() and should_perform_rnaseq_metrics:
            genome_valid = False
            print("The REF_FLAT_FILE file was not found")
        
        if not bed_file.exists() and should_perform_get_fragment_size:
            genome_valid = False
            print("The BED_FILE file was not found")
        
        if not rRna_interval_list.exists():
            genome_valid = False
            print("The RRNA_INTERVAL_LIST file was not found")
        
        if not genome_fasta_file.exists():
            genome_valid = False
            print("The GENOME_FASTA_FILE file was not found")
        
        if not gtf_file.exists():
            genome_valid = False
            print("The GTF_FILE file was not found")
            
        if not genome_valid:
            print(f"Searching config file: {config_filename}")
            raise ValueError("Unable to find one or more genome-related files")
    
        return genome_valid

class _ControlData:
    @classmethod
    def ingest(cls, control: Union[str, Path]) -> dict[str, dict[str, list[str]]]:
        """
        This function is responsible for taking data from the control file and creating a usable data format
        :param control: The control file to read from
        :return: A dictionary containing each cell type, their studies, and each studies' replicates
        """
        schema: dict[str, dict[str, list[str]]] = {}
    
        with open(control, "r") as i_stream:
            dialect = csv.Sniffer().sniff(i_stream.read(1024))
            i_stream.seek(0)
            reader = csv.reader(i_stream, delimiter=str(dialect.delimiter))
            for line in reader:
                sample = line[1]  # "naiveB_S1R1"
                cell_type = sample.split("_")[0]  # naiveB
                tag = sample.split("_")[1]  # S1R1
                study = re.search(r"S\d+", tag).group()  # S1
                run = re.search(r"R\d+(r\d+)?", tag).group()  # R1 (with optional replicate, i.e., R1r1)
            
                if cell_type not in schema.keys():
                    schema[cell_type] = {}
            
                cell_studies = schema[cell_type]
                if study not in cell_studies.keys():
                    cell_studies[study] = []
            
                study_list: list[str] = cell_studies[study]
                if run not in study_list:
                    study_list.append(run)
    
        return schema

    @classmethod
    def is_valid(cls, schema: dict[str, dict[str, list[str]]]) -> Union[bool, None]:
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


def validate(config: dict) -> bool:
    config_file: str = config["MASTER_CONTROL"]

    # These are the "master" controls
    # If any of these are invalid, this function will return False
    control_valid: bool = True
    genome_valid: bool = True
    
    if config["BYPASS_REPLICATE_VALIDATION"]:
        control_valid = True
    else:
        # Validate the control file (file containing SRR, sample, etc.)
        schema = _ControlData.ingest(config_file)
        if not _ControlData.is_valid(schema):
            control_valid = False
    
    if config["BYPASS_GENOME_VALIDATION"]:
        genome_valid = True
    else:
        genome = _GenomeData.validate(config_file, config)
        if not genome:
            genome_valid = False
    
    return control_valid and genome_valid
