import plotly.express as px
from plotly.graph_objs import Figure
from pathlib import Path
import typing
import pandas as pd

def read_file():
    """
    This is a temporary function to read benchmark files into a pandas dataframe
    :return: 
    """
    file: Path = Path("/Users/joshl/PycharmProjects/FastqToGeneCounts/benchmarking/preB/contaminant_screen/preB_S1R1_1.benchmark")
    df = pd.read_csv(file, delimiter="\t")
    
    print(df)
    df = df.rename(columns={"s": "seconds"})
    df = df["seconds"]
    print(df)


def create_dataframe(files: dict[str, list[Path]]) -> pd.DataFrame:
    """
    This function is responsible for aggregating all the input benchmark files into a pandas dataframe
    :param files: A list of benchmark files as Path objects
    :return: A pandas dataframe
    """

    dataframes: pd.DataFrame = pd.DataFrame()
    column_names: list[str] = []
    for rule_name, files in files.items():
        rule_df: pd.DataFrame = pd.DataFrame()
        for file in files:
            
            # Only get the time in seconds
            df = pd.read_csv(file, delimiter="\t")["s"]
            rule_df = pd.concat([rule_df, df], axis=0, ignore_index=True)
        
        # For each rule dataframe created, concatentate it to the main dataframe
        dataframes = pd.concat([dataframes, rule_df], axis=1, ignore_index=True)
        
        # Crate a title-cased rule name by splitting the rule name on underscores and joining the words with spaces, then uppercasing letters after a space
        new_column_name: str = " ".join(rule_name.split("_")).title()
        column_names.append(new_column_name)
    
    # Reset the column names
    dataframes.columns = column_names
    
    # Calculate the average of each rule and palce it at the bottom of the dataframe
    dataframes.loc["Average"] = dataframes.mean()
    
    return dataframes
    
def plot_data(df: pd.DataFrame):
    """
    This function is responsible for plotting the benchmark data using Plotly
    :param df: A pandas dataframe containing the benchmark data from create_dataframe
    :return:
    """
    # Create a bar graph of the "Average" row if the value is greater than 10
    df: pd.DataFrame = df.loc["Average"]
    df: pd.DataFrame = df[df > 10]
    df: pd.DataFrame = df.sort_values(ascending=True)
    fig: px.bar = px.bar(
        df, x=df.index, y=df.values,
        labels={"x": "Rule Name", "y": "Time (seconds)"},
        title="Benchmarking Rules Requiring Greater than 10 seconds to Run"
    )
    fig.show()
        

def main():
    # Get all benchmark files in the current directory
    benchmark_files: typing.Generator[Path, None, None] = (file for file in Path().rglob("*.benchmark"))
    
    # Create a dictionary of rule names and their corresponding benchmark files
    rule_benchmark_files: dict[str, list[Path]] = {}
    for file in benchmark_files:
        rule_name = file.parent.stem
        if rule_name not in rule_benchmark_files:
            rule_benchmark_files[rule_name] = []
        rule_benchmark_files[rule_name].append(file)
        
    benchmark_df: pd.DataFrame = create_dataframe(rule_benchmark_files)
    plot_data(benchmark_df)


if __name__ == '__main__':
    main()
