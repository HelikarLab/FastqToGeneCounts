def trim(config: dict):  # QC
    return str(config["PERFORM_TRIM"]).lower() == "true"

def screen(config: dict):  # QC
    return str(config["PERFORM_SCREEN"]).lower() == "true"


def prefetch(config: dict):
    return str(config["PERFORM_PREFETCH"]).lower() == "true"


def get_insert_size(config: dict):
    return str(config["PERFORM_GET_INSERT_SIZE"]).lower() == "true"


def get_fragment_size(config: dict):  # for zFPKM QC
    return str(config["PERFORM_GET_FRAGMENT_SIZE"]).lower() == "true"


def get_rnaseq_metrics(config: dict):  # QC
    return str(config["PERFORM_GET_RNASEQ_METRICS"]).lower() == "true"
