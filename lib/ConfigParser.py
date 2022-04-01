"""
This file will be responsible for parsing an input configuration file
This is needed after splitting the master Snakefile into various sub-snakefiles
These new sub-snakefiles do not have access to the checker-functions (like perform_prefetch) that are located in the master snakefile
As a result, we must parse the config file from snakemake, and
"""

import yaml


class ConfigParser:
    # Use __new__ here so we can return the _yaml_config immediately after instantiating the class
    def __new__(cls, config_path, *args, **kwargs):
        _i_stream = open(config_path, "r")
        _yaml_config = yaml.safe_load(_i_stream)
        return _yaml_config


if __name__ == '__main__':

    config = ConfigParser("/Users/joshl/PycharmProjects/TestSnakemake/config.yaml")
    print(config)
    print(config["ROOT_DIR"])
    print(config["TOP"]["LOWER"])
