import pyodbc
import pandas as pd
import geopandas as gpd
import ee
import os
import warnings

# Suppress the PyODBC UserWarning
warnings.filterwarnings("ignore", message="pandas only supports SQLAlchemy connectable")

# Suppress the DtypeWarning (Mixed Data Types in Pandas)
warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)

###### LOAD DATA ####################################################################################################
def load_data(file_path, sheet_name=None, table_dict=None, encoding="utf-8", delimiter=","):
    """
    Loads data from various file formats including CSV, Excel (XLSX), Access (MDB), Shapefile (SHP), and Geopackage (GPKG).

    Parameters:
        file_path (str): The path to the file.
        sheet_name (str, optional): The sheet name for Excel files (XLSX).
        table_dict (list of str, optional): The table name(s) for Access (MDB) and Geopackage (GPKG) databases.
        encoding (str, optional): Encoding format for CSV files (default: "utf-8").
        delimiter (str, optional): Delimiter for CSV files (default: ",").

    Returns:
        pandas.DataFrame, dictionary of pandas.DataFrame or geopandas.GeoDataFrame: The loaded data.
    """
    # Check if the path is a directory and ends with "shp" (case insensitive)
    if os.path.isdir(file_path) and file_path.lower().endswith("shp"):
        return gpd.read_file(file_path)  # Pass entire directory to geopandas

    else:
        file_extension = file_path.split('.')[-1].lower()

        if file_extension == 'csv':
            return pd.read_csv(file_path, encoding=encoding, delimiter=delimiter)

        elif file_extension in ['xls', 'xlsx']:
            if sheet_name is None:
                return pd.read_excel(file_path)
            else:
                return pd.read_excel(file_path, sheet_name=sheet_name)

        elif file_extension == 'mdb':
            if table_dict is None:
                raise ValueError("Please specify a table_name to load data from an MDB file.")

            # Define connection string
            conn_str = (
                r"DRIVER={Microsoft Access Driver (*.mdb, *.accdb)};"
                f"DBQ={file_path};"
            )

            # Connect to the database
            conn = pyodbc.connect(conn_str)

            # Function to read a table into a DataFrame
            def read_table_to_df(table_name, columns):
                """
                Read a table (e.g mdb) into a DataFrame from the database.
                Parameters:
                table_name (str): The name of the table to read.
                Returns:
                pd.DataFrame: The table as a DataFrame.
                """
                query = f"SELECT {columns} FROM {table_name}"
                return pd.read_sql(query, conn)

                # Read tables into DataFrames

            dfs = {}
            for table_name, columns in table_dict.items():
                col_str = ', '.join(columns)  # Format column names into a SQL query
                dfs[table_name] = read_table_to_df(table_name, col_str)

            # Close the connection
            conn.close()
            return dfs

        elif file_extension == 'gpkg':
            return gpd.read_file(file_path)

        else:
            raise ValueError(f"Unsupported file format: {file_extension}")

