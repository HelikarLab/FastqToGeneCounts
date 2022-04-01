"""
This function will be responsible for implementing the simple True/False function checkers
    originally located at the top of the master "Snakefile"

This is required after splitting the master Snakefile into multiple sub-workflows, as these new sub-workflows
    do not have access to the functions implemented by the master Snakefile (such as def perform_prefetch)

As a result of this, we must implement these functions into a file that can be accessible by any Snakefile
"""
from ConfigParser import ConfigParser

config = ConfigParser("config.yaml")


def perform_prefetch() -> bool:
    if str(config["PERFORM_PREFETCH"]).lower() == "true":
        return True
    else:
        return False


def perform_trim():  # QC
    if str(config["PERFORM_TRIM"]).lower() == "true":
        return True
    else:
        return False


def perform_screen():  # QC
    if str(config["PERFORM_SCREEN"]).lower() == "true":
        return True
    else:
        return False


def perform_get_insert_size():
    if str(config["PERFORM_GET_INSERT_SIZE"]).lower() == "true":
        return True
    else:
        return False


def perform_get_fragment_size():  # for zFPKM QC
    print(str(config["PERFORM_GET_FRAGMENT_SIZE"]).lower())
    if str(config["PERFORM_GET_FRAGMENT_SIZE"]).lower() == "true":
        return True
    else:
        return False


def perform_get_rnaseq_metrics():  # QC
    if str(config["PERFORM_GET_RNASEQ_METRICS"]).lower() == "true":
        return True
    else:
        return False
