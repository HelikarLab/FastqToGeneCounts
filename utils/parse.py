# ruff: noqa: D102 T201
import csv
import re
import shutil
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

import pandas as pd

from utils.download_genome import get_latest_release, species_from_taxon


@dataclass(frozen=True, slots=True)
class Perform:
    dump_fastq: bool
    trim: bool
    contaminant_screen: bool
    rnaseq_metrics: bool
    insert_size: bool
    fragment_size: bool


@dataclass(frozen=True, slots=True)
class Validation:
    bypass_replicate_validation: bool
    bypass_genome_validation: bool


def print_key_value_table(title: str, rows: Iterable[tuple[str, object]], *, center_block: bool = True) -> None:
    rows = [(str(k), "" if v is None else str(v)) for k, v in rows]
    outer_pad = 5
    left_width = max((len(k) for k, _ in rows), default=0)
    right_width = max((len(v) for _, v in rows), default=0)

    body_lines: list[str] = []
    for k, v in rows:
        new_k = f"|{' ' * outer_pad}{k:<{left_width}}"
        new_v = f"{v:<{right_width}}{' ' * outer_pad} |"
        body_lines.append(f"{new_k}: {new_v}")

    # body_lines = [f"| {k:<{left_width}}: {v:<{right_width}} |" for k, v in rows]
    interior_width = 1 + left_width + right_width + 1 + (outer_pad * 2) - 1  # " " + left + ": " + right + " "
    title = f"| {title:^{interior_width}} |"

    width = max([len(title), *(len(s) for s in body_lines)], default=len(title))
    top = "=" * width
    mid = f"| {'-' * interior_width} |"

    lines = [top, title, mid, *body_lines, top]
    prefix = ""
    if center_block:
        cols = shutil.get_terminal_size(fallback=(80, 24)).columns
        pad = max((cols - width) // 2, 0)
        prefix = " " * pad
    lines = [prefix + line for line in lines]
    print()
    print("\n".join(lines))
    print()


@dataclass(frozen=True, slots=True)
class Genome:
    species_dir: Path
    contaminants_dir: Path
    taxon_id: int
    version: str
    show_progress: bool
    type: Literal["primary_assembly", "toplevel"]
    _ensembl_release: str = field(init=False)

    @property
    def ensembl_release(self):
        return self._ensembl_release

    def __post_init__(self):
        """Post initialization to set the ensembl_release attribute."""
        if self.version == "latest":
            val = f"release-{get_latest_release()}"
        elif self.version.startswith("release"):
            val = self.version
        elif self.version.isdigit():
            val = f"release-{self.version}"
        else:
            raise ValueError(
                "Invalid GENOME VERSION in config.yaml file. "
                "Valid options are: 'latest', 'release-###' (i.e., 'release-112'), or an integer (i.e., 112)."
            )

        if self.type not in ["primary_assembly", "toplevel"]:
            raise ValueError(
                f"Invalid GENOME FASTA_TYPE in config.yaml file. Valid options are: 'primary_assembly' or 'toplevel'; got: '{self.type}'."
            )

        object.__setattr__(self, "_ensembl_release", val)
        # Append `-{release-value}` to the species directory
        object.__setattr__(self, "species_dir", Path(self.species_dir.as_posix() + f"_{self.ensembl_release}"))


@dataclass(frozen=True, slots=True)
class Config:
    sample_filepath: Path
    root: Path
    data_root: Path
    temp_root: Path
    como_root: Path
    logs_root: Path
    experiment_name: str
    local_fastq_filepath: Path | None
    perform: Perform
    validation: Validation
    benchmark_count: int
    genome: Genome
    species_name: str
    benchmark_dir: Path
    _fastq_files: list[Path] = field(default_factory=list, init=False)

    @staticmethod
    def _validate_config(config: dict[str, Any]):  # noqa: C901
        if config["MASTER_CONTROL"] == "":
            raise ValueError("MASTER_CONTROL cannot be an empty string.")
        if not Path(config["MASTER_CONTROL"]).exists():
            raise FileNotFoundError(f"MASTER_CONTROL path does not exist: {config['MASTER_CONTROL']}")
        
        if not str(config["BENCHMARK_TIMES"]).isdigit() or int(config["BENCHMARK_TIMES"]) < 0:
            raise ValueError("BENCHMARK_TIMES must be a non-negative integer.")

        if config["LOCAL_FASTQ_FILES"] != "" and not Path(config["LOCAL_FASTQ_FILES"]).exists():
            raise FileNotFoundError(f"LOCAL_FASTQ_FILES path does not exist: {config['LOCAL_FASTQ_FILES']}")

        if not str(config["GENOME"]["TAXONOMY_ID"]).isdigit():
            raise ValueError("GENOME TAXONOMY_ID must be an integer.")

        str_true_false = ["true", "false"]
        if str(config["PERFORM_DUMP_FASTQ"]).lower() not in str_true_false:
            raise ValueError("PERFORM_DUMP_FASTQ must be 'true' or 'false'.")
        if str(config["PERFORM_TRIM"]).lower() not in str_true_false:
            raise ValueError("PERFORM_TRIM must be 'true' or 'false'.")
        if str(config["PERFORM_SCREEN"]).lower() not in str_true_false:
            raise ValueError("PERFORM_SCREEN must be 'true' or 'false'.")
        if str(config["PERFORM_GET_RNASEQ_METRICS"]).lower() not in str_true_false:
            raise ValueError("PERFORM_GET_RNASEQ_METRICS must be 'true' or 'false'.")
        if str(config["PERFORM_GET_INSERT_SIZE"]).lower() not in str_true_false:
            raise ValueError("PERFORM_GET_INSERT_SIZE must be 'true' or 'false'.")
        if str(config["PERFORM_GET_FRAGMENT_SIZE"]).lower() not in str_true_false:
            raise ValueError("PERFORM_GET_FRAGMENT_SIZE must be 'true' or 'false'.")
        if str(config["BYPASS_REPLICATE_VALIDATION"]).lower() not in str_true_false:
            raise ValueError("BYPASS_REPLICATE_VALIDATION must be 'true' or 'false'.")
        if str(config["BYPASS_GENOME_VALIDATION"]).lower() not in str_true_false:
            raise ValueError("BYPASS_GENOME_VALIDATION must be 'true' or 'false'.")
        if str(config["GENOME"]["SHOW_PROGRESS"]).lower() not in str_true_false:
            raise ValueError("GENOME SHOW_PROGRESS must be 'true' or 'false'.")

    @classmethod
    def create(cls, config: dict[str, Any]) -> "Config":
        """Create a configuration object."""
        Config._validate_config(config)

        root = Path(config["ROOTDIR"])
        experiment_name = config["EXPERIMENT_NAME"]
        fastq_files = config["LOCAL_FASTQ_FILES"]
        taxon_id: int = int(config["GENOME"]["TAXONOMY_ID"])
        species_name: str = species_from_taxon(taxon_id=taxon_id)

        sample_filepath = Path(config["MASTER_CONTROL"])
        data_root = Path(root, experiment_name, "data")
        temp_root = Path(root, experiment_name, "temp")
        como_root = Path("COMO_input", experiment_name)
        logs_root = Path(config["LOG_DIR"], experiment_name)
        species_dir = Path(config["GENOME"]["SAVE_DIR"], species_name)

        return cls(
            sample_filepath=Path(config["MASTER_CONTROL"]),
            root=root,
            data_root=data_root,
            temp_root=temp_root,
            como_root=como_root,
            logs_root=logs_root,
            experiment_name=experiment_name,
            local_fastq_filepath=Path(fastq_files) if fastq_files else None,
            species_name=species_name,
            benchmark_count=int(config["BENCHMARK_TIMES"]),
            benchmark_dir=Path(config["BENCHMARK_DIR"], experiment_name),
            perform=Perform(
                dump_fastq=config["PERFORM_DUMP_FASTQ"],
                trim=config["PERFORM_TRIM"],
                contaminant_screen=config["PERFORM_SCREEN"],
                rnaseq_metrics=config["PERFORM_GET_RNASEQ_METRICS"],
                insert_size=config["PERFORM_GET_INSERT_SIZE"],
                fragment_size=config["PERFORM_GET_FRAGMENT_SIZE"],
            ),
            validation=Validation(
                bypass_replicate_validation=config["BYPASS_REPLICATE_VALIDATION"],
                bypass_genome_validation=config["BYPASS_GENOME_VALIDATION"],
            ),
            genome=Genome(
                species_dir=species_dir,
                contaminants_dir=Path(config["ROOTDIR"], "FastQ_Screen_Genomes"),
                taxon_id=taxon_id,
                version=config["GENOME"]["VERSION"],
                type=config["GENOME"]["FASTA_TYPE"],
                show_progress=config["GENOME"]["SHOW_PROGRESS"],
            ),
        )

    def __post_init__(self):
        """Post initialization checks."""
        if self.local_fastq_filepath is not None:
            object.__setattr__(
                self,
                "_fastq_files",
                list(self.local_fastq_filepath.rglob("*.fastq")) + list(self.local_fastq_filepath.rglob("*.fastq.gz")),
            )
        self.root.mkdir(parents=True, exist_ok=True)

    def fastq_files(self, filter_by: str = "") -> list[Path]:
        """Return a list of filepaths matching an optional filter value.

        If `filter` is an empty string, all files will be returned.
        The filter is a simple `if filter in str`, and does not support regex, glob, etc.
        """
        if not filter_by:
            return self._fastq_files
        return [i for i in self._fastq_files if filter_by in i.as_posix()]


class SampleData:
    def __init__(self, sample_filepath: Path):
        """Parse the provided filepath into a SampleData object."""
        with sample_filepath.open() as i_stream:
            # Get the delimiter from the master control file (in case it is not a comma)
            # from: https://stackoverflow.com/questions/16312104
            delimiter: str = csv.Sniffer().sniff(i_stream.readline().rstrip("\n")).delimiter
            i_stream.seek(0)  # reset i_stream read buffer
            self._sample_df: pd.DataFrame = pd.read_csv(sample_filepath, header=0, delimiter=delimiter)

        self._pairs: pd.DataFrame = self._sample_df["sample"].astype(str).str.extract(r"^(?P<tissue>.+)_(?P<tag>S\d+R\d+(?:r\d+)?)$")
        if self._pairs.isna().any().any():
            na_value = self._sample_df[self._pairs.isna().any(axis=1)]
            print(na_value)
            raise ValueError(
                "Some sample names in the MASTER_CONTROL file do not follow the expected format '<tissue>_<tag>' (e.g., effectorcd8_S1R1). "
                "The items have been printed on the previous line."
            )

        self._sample_names: list[str] = self._sample_df["sample"].to_list()
        self._tissues: list[str] = self._pairs["tissue"].tolist()
        self._tags: list[str] = self._pairs["tag"].tolist()
        self._ends: list[str] = self._sample_df["endtype"].astype(str).str.upper().tolist()
        self._studies: list[str] = self._pairs["tag"].str.extract(r"^(S\d+)")[0].tolist()
        self._tissues_paired: list[str] = []
        self._tags_paired: list[str] = []
        self._studies_paired: list[str] = []
        self._ends_paired: list[str] = []
        for row in self._sample_df.itertuples(index=False):
            if row.endtype not in ["PE", "SE"]:
                raise ValueError(f"Unexpected 'endtype': '{row.endtype}'. Expected 'PE' or 'SE'.")

            _tissue, _tag = row.sample.split("_")
            _study: re.Match[str] | None = re.match(r"S\d+", _tag)
            if not _study:
                raise ValueError(f"Unexpected tag format: '{_tag}'. Expected format 'S#R#' (e.g., 'S1R1').")
            _study_match = _study.group()

            self._tissues_paired += [_tissue, _tissue] if row.endtype == "PE" else [_tissue]
            self._tags_paired += (
                [_tag, _tag]
                if row.endtype == "PE"
                else [
                    _tag,
                ]
            )
            self._studies_paired += [_study_match, _study_match] if row.endtype == "PE" else [_study_match]
            self._ends_paired += ["1", "2"] if row.endtype == "PE" else ["S"]

    @property
    def samples(self) -> pd.DataFrame:
        return self._sample_df

    @property
    def pairs(self) -> pd.DataFrame:
        return self._pairs

    @property
    def sample_names(self) -> list[str]:
        return self._sample_names

    @property
    def tissues(self) -> list[str]:
        return self._tissues

    @property
    def tags(self) -> list[str]:
        return self._tags

    @property
    def ends(self) -> list[str]:
        return self._ends

    @property
    def studies(self) -> list[str]:
        return self._studies

    @property
    def tissues_paired(self) -> list[str]:
        return self._tissues_paired

    @property
    def tags_paired(self) -> list[str]:
        return self._tags_paired

    @property
    def studies_paired(self) -> list[str]:
        return self._studies_paired

    @property
    def ends_paired(self) -> list[str]:
        return self._ends_paired
