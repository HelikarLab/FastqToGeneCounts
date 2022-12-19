import plotly.express as px
import plotly.graph_objects as go
from plotly.graph_objs import Figure
from pathlib import Path
import typing
import pandas as pd
from _datetime import datetime

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
            
            # Only get the "time" column
            try:
                df = pd.read_csv(file, delimiter="\t")["s"]
            except pd.errors.EmptyDataError:
                print(f"No columns found in file {file}")
                exit(1)
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
    
def plot_data(df: pd.DataFrame) -> Figure:
    """
    This function is responsible for plotting the benchmark data using Plotly
    :param df: A pandas dataframe containing the benchmark data from create_dataframe
    :return:
    """
    # Create a bar graph of the "Average" row if the value is greater than 10
    df: pd.DataFrame = df.loc["Average"]
    # df: pd.DataFrame = df[df > 10]
    df: pd.DataFrame = df.sort_values(ascending=True)
    fig: px.bar = px.bar(
        df, x=df.index, y=df.values,
        labels={"x": "Rule Name", "y": "Time (seconds)"},
        title="Average Runtime for a Single File"
    )
    
    # Add a bar that is the total of the Averages row
    fig.add_trace(
        go.Bar(
            x=["Total"],
            y=[df.sum()],
            name="Total"
        )
    )
    
    return fig
    
def save_plot(fig: Figure) -> None:
    """
    This function is responsible for saving the plotly figure to a file
    :param fig: A plotly figure
    :return: None
    """
    # Write the plotly figure to an html file with the following name: benchmarking_{date}_{time}.html
    fig.write_html(f"benchmarking_{datetime.now().strftime('%Y-%m-%d')}.html")

def main():
    # Get all benchmark files in the current directory
    benchmark_files: list[Path] = [file for file in Path().rglob("*.benchmark")]

    # Create a dictionary of rule names and their corresponding benchmark files
    rule_benchmark_files: dict[str, list[Path]] = {}
    for file in benchmark_files:
        rule_name = file.parent.stem
        if rule_name not in rule_benchmark_files:
            rule_benchmark_files[rule_name] = []
        rule_benchmark_files[rule_name].append(file)
        
    benchmark_df: pd.DataFrame = create_dataframe(rule_benchmark_files)
    plotly: Figure = plot_data(benchmark_df)
    save_plot(plotly)
    
    
if __name__ == '__main__':
    main()
