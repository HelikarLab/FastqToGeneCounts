# ruff: noqa: D102

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

from utils.download_genome import Utilities


@dataclass(frozen=True, slots=True)
class Genome:
    species_dir: Path
    contaminants_dir: Path
    taxon_id: int
    version: str
    show_progress: bool
    _ensembl_release: str = field(init=False)

    @property
    def ensembl_release(self):
        return self._ensembl_release

    def __post_init__(self):
        """Post initialization to set the ensembl_release attribute."""
        if self.version == "latest":
            val = f"release-{Utilities.get_latest_release()}"
        elif self.version.startswith("release"):
            val = self.version
        elif self.version.isdigit():
            val = f"release-{self.version}"
        else:
            raise ValueError(
                "Invalid GENOME VERSION in config.yaml file. "
                "Valid options are: 'latest', 'release-###' (i.e., 'release-112'), or an integer (i.e., 112)."
            )

        object.__setattr__(self, "_ensembl_release", val)


@dataclass(frozen=True, slots=True)
class Perform:
    prefetch: bool
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

    @classmethod
    def create(cls, config: dict[str, Any]) -> "Config":
        """Create a configuration object."""
        root = Path(config["ROOTDIR"])
        experiment_name = config["EXPERIMENT_NAME"]
        fastq_files = config["LOCAL_FASTQ_FILES"]
        taxon_id: int = int(config["GENOME"]["TAXONOMY_ID"])
        species_name: str = Utilities.get_species_from_taxon(taxon_id=taxon_id)

        return cls(
            sample_filepath=Path(config["MASTER_CONTROL"]),
            root=root,
            data_root=Path(root, experiment_name, "data"),
            temp_root=Path(root, experiment_name, "temp"),
            como_root=Path("COMO_input", experiment_name),
            logs_root=Path(config["LOG_DIR"], experiment_name),
            experiment_name=experiment_name,
            local_fastq_filepath=Path(fastq_files) if fastq_files else None,
            species_name=species_name,
            benchmark_count=config["BENCHMARK_TIMES"],
            benchmark_dir=Path(config["BENCHMARK_DIR"]),
            perform=Perform(
                prefetch=config["PERFORM_PREFETCH"],
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
                species_dir=Path(config["GENOME"]["SAVE_DIR"], species_name),
                contaminants_dir=Path(config["ROOTDIR"], "FastQ_Screen_Genomes"),
                taxon_id=taxon_id,
                version=config["GENOME"]["VERSION"],
                show_progress=config["GENOME"]["SHOW_PROGRESS"],
            ),
        )

    def __post_init__(self):
        """Post initialization checks."""
        if not self.sample_filepath.exists():
            raise FileNotFoundError(f"Control file path does not exist: {self.sample_filepath}")
        if self.benchmark_count < 0:
            raise ValueError("Benchmark count must be non-negative.")
        if self.local_fastq_filepath and not self.local_fastq_filepath.exists():
            raise FileNotFoundError(f"Local FASTQ file path does not exist: {self.local_fastq_filepath}")

        self.root.mkdir(parents=True, exist_ok=True)


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
            raise ValueError(
                "Some sample names in the MASTER_CONTROL file do not follow the expected format '<tissue>_<tag>' (e.g., effectorcd8_S1R1)."
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


if __name__ == "__main__":
    SampleData(Path("/Users/joshl/Projects/FastqToGeneCounts/controls/application_note.csv"))
