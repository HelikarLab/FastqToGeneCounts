from dataclasses import dataclass
from pathlib import Path
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
from tqdm import tqdm


@dataclass(frozen=True, slots=True, weakref_slot=False)
class Pearson:
    stat: float
    pval: float


@dataclass(frozen=True, slots=True, weakref_slot=False)
class Spearman:
    stat: float
    pval: float


@dataclass(frozen=True, slots=True, weakref_slot=False)
class KSTest:
    stat: float
    pval: float
    location: float
    sign: int | np.uint8

    def __post_init__(self):
        object.__setattr__(self, "stat", float(self.stat))
        object.__setattr__(self, "pval", float(self.pval))
        object.__setattr__(self, "location", float(self.location))
        object.__setattr__(self, "sign", int(self.sign))


class Statistics:
    def __init__(self, raw: pd.DataFrame, original_suffix: str = "original", reproduction_suffix: str = "reproduction"):
        self._raw: pd.DataFrame = raw
        self.original_suf: str = original_suffix.removeprefix("_")
        self.reproduction_suf: str = reproduction_suffix.removeprefix("_")
        self._names: list[str] = ["_".join(i.split("_")[:2]) for i in self._raw.columns]

    @staticmethod
    def log1p(df: pd.DataFrame) -> pd.DataFrame:
        return pd.DataFrame(data=np.log1p(df.values.copy()), index=df.index.copy(), columns=df.columns.copy())

    def pearson(self) -> tuple[dict[str, Pearson], dict[str, Pearson]]:
        """Perform Pearson statistic calculations on original vs reproduction data.

        :returns: A tuple of dictionaries:
            1) Pearson calculations on raw data
            2) Pearson calculations on log1p-transformed data
        """

        def _do_pearson(_df: pd.DataFrame) -> dict[str, Pearson]:
            correlations: dict[str, Pearson] = {}
            for name in self._names:
                result = scipy.stats.pearsonr(self._raw[f"{name}_{self.original_suf}"], self._raw[f"{name}_{self.reproduction_suf}"])
                correlations[name] = Pearson(stat=result.statistic, pval=result.pvalue)
            return correlations

        df: pd.DataFrame = self._raw.copy()
        return _do_pearson(df), _do_pearson(self.log1p(df))

    def spearman(self) -> tuple[dict[str, Spearman], dict[str, Spearman]]:
        """Perform Spearman statistic calculations on original vs reproduction data.

        :returns: A tuple of dictionaries:
            1) Spearman calculations on raw data
            2) Spearman calculations on log1p-transformed data
        """

        def _do_spearman(_df: pd.DataFrame) -> dict[str, Spearman]:
            names: list[str] = ["_".join(i.split("_")[:2]) for i in _df.columns]
            correlations: dict[str, Spearman] = {}
            for name in names:
                result = scipy.stats.spearmanr(_df[f"{name}_{self.original_suf}"], _df[f"{name}_{self.reproduction_suf}"], nan_policy="omit")
                correlations[name] = Spearman(stat=result.statistic, pval=result.pvalue)
            return correlations

        df = self._raw.copy()
        return _do_spearman(df), _do_spearman(self.log1p(df))

    def kolmogorov_smirnov(self) -> tuple[dict[str, KSTest], dict[str, KSTest]]:
        """Perform Kolmogorov-Smirnov statistic calculations on original vs reproduction data.

        :returns: A tuple of dictionaries:
            1) Kolmogorov-Smirnov calculations on raw data
            2) Kolmogorov-Smirnov calculations on log1p-transformed data
        """

        def _do_ks(_df: pd.DataFrame) -> dict[str, KSTest]:
            names: list[str] = ["_".join(i.split("_")[:2]) for i in _df.columns]

            correlations: dict[str, KSTest] = {}
            for name in self._names:
                result = scipy.stats.kstest(
                    rvs=_df[f"{name}_{self.reproduction_suf}"],
                    cdf=_df[f"{name}_{self.original_suf}"],
                    alternative="two-sided",
                )
                correlations[name] = KSTest(stat=result.statistic, pval=result.pvalue, location=result.statistic_location, sign=result.statistic_sign)
            return correlations

        df = self._raw.copy()
        return _do_ks(df), _do_ks(self.log1p(df))

    def plot_distribution(self, sample_name: str | None = None, *, log1p: bool = False):
        sample_name = sample_name or self._raw.columns[0].removesuffix(f"_{self.original_suf}")
        original_name = f"{sample_name}_{self.original_suf}"
        reproduction_name = f"{sample_name}_{self.reproduction_suf}"

        orig_data, repro_data = self._raw[original_name].to_numpy(copy=True), self._raw[reproduction_name].to_numpy(copy=True)
        if log1p:
            orig_data, repro_data = np.log1p(orig_data), np.log1p(repro_data)

        fig, ax = plt.subplots(nrows=1, ncols=1)
        limits = [
            np.min([orig_data, repro_data]),
            np.max([orig_data, repro_data]),
        ]

        ax.plot(limits, limits, color="black", linestyle="-", linewidth=0.75, zorder=0)
        plt.scatter(x=orig_data, y=repro_data, s=10, zorder=1)
        ax.set_title(f"Log1p Gene Counts for Sample: '{sample_name}'" if log1p else f"Raw Gene Counts for Sample: '{sample_name}'")
        ax.set_xlabel(f"Gene Counts: '{original_name}'")
        ax.set_ylabel(f"Gene Counts: '{reproduction_name}'")
        fig.show()
        plt.close(fig)


def read_excel(filepath: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read the excel datafile.

    The file should have sheets with names "Original" and "Reproduction"

    :param: The file to read from
    :returns: A tuple of pd.DataFrames, where:
        1) The original data
        2) The reproduction data
    """
    with pd.ExcelFile(filepath) as xlsx:
        if "Original" not in xlsx.sheet_names:
            raise ValueError("Excel file must contain a sheet named 'Original'")
        if "Reproduction" not in xlsx.sheet_names:
            raise ValueError("Excel file must contain a sheet named 'Reproduction'")
        original: pd.DataFrame = pd.read_excel(xlsx, "Original", skiprows=3, header=0)
        reproduction: pd.DataFrame = pd.read_excel(xlsx, "Reproduction", skiprows=3, header=0)

    original: pd.DataFrame = original.drop(columns=["Unnamed: 0", "Sample Names"])
    original: pd.DataFrame = original.set_index("Gene ID")

    reproduction: pd.DataFrame = reproduction.drop(columns=["Unnamed: 0", "Sample Names"])
    reproduction: pd.DataFrame = reproduction.set_index("Gene ID")
    return original, reproduction


def df_from_quantsf_files(path: Path | str, suffix: str) -> pd.DataFrame:
    suffix = suffix.removeprefix("_")
    files = list(Path(path).iterdir())
    sample_names = [f"{i.name.removesuffix('_quant.sf')}_{suffix}" for i in files]

    dfs: list[pd.DataFrame] = []
    for file, name in tqdm(zip(files, sample_names, strict=True), total=len(files), desc=f"Reading '{suffix}' files"):
        df = pd.read_csv(file, sep="\t")
        df = df.drop(columns=["Length", "EffectiveLength", "TPM"])
        df.columns = ["name", name]
        df = df.set_index("name", drop=True)
        dfs.append(df)
    return pd.concat(dfs, axis=1)


def main():
    nfcore_dir = Path(__file__).parent / "salmon_quant" / "nfcore"
    ftgc_dir = Path(__file__).parent / "salmon_quant" / "ftgc"
    
    if not nfcore_dir.exists():
        raise FileNotFoundError(f"Directory not found: {nfcore_dir}")
    if not ftgc_dir.exists():
        raise FileNotFoundError(f"Directory not found: {ftgc_dir}")
    
    nfcore = df_from_quantsf_files(nfcore_dir, suffix="nfcore")
    reproduction = df_from_quantsf_files(ftgc_dir, suffix="ftgc")
    if len(nfcore.columns) != len(reproduction.columns):
        raise ValueError("Mismatched number of columns between original and reproduction dataframes.")
    
    merged = nfcore.merge(reproduction, left_index=True, right_index=True)

    stat = Statistics(raw=merged, original_suffix="nfcore", reproduction_suffix="reproduction")
    raw_pearson, log1p_pearson = stat.pearson()
    print("\n\nRaw Pearson")
    pprint(raw_pearson)
    
    print("\nLog1p Pearson")
    pprint(log1p_pearson)

    raw_spearman, log1p_spearman = stat.spearman()
    print("\n\nRaw Spearman")
    print(raw_spearman)
    print("\nLog1p Spearman")
    pprint(log1p_spearman)

    raw_ks, log1p_ks = stat.kolmogorov_smirnov()
    print("\n\nRaw K-S")
    print(raw_ks)
    print("\nLog1p K-S")
    pprint(log1p_ks)


if __name__ == "__main__":
    main()
