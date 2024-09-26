from os import environ, getcwd, chdir
from os.path import join
from subprocess import call
# import rpy2's package module
import rpy2.robjects.packages as rpackages
import rpy2.robjects as robjects
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.packages import importr

import pandas as pd

from mmonitor.config import images_path, r_path

# static paths
# check config.py for file locations
_CSV_FILE = join(r_path, 'abundances.csv')
# _IMAGE_FILE = join(images_path, 'horizon.png')
_HORIZON_R_SCRIPT = join(r_path, 'horizon.r')


def generate_image(df: pd.DataFrame, width: int, height: int) -> None:
    """
    Generate a horizon plot of taxonomy abundances.
    Check static paths for file locations.

    WARNING: - Local files are being created to pass data
             - Files that should be static are created dynamically
             - Scripts are being called in subprocesses
             - Environment variables are being passed between processes

    Future considerations:
             - Find a way to integrate the R script into Python
             - Find a Python-native method to generate horizon plots
             - Find a way to generate horizon plots in Plotly
    """

    horizon_scale = _write_csv(df)
    _export_env_vars(width, height, horizon_scale)
    _call_r_script()


def _write_csv(df: pd.DataFrame) -> int:
    """
    Transpose mmonitor db to a Dataframe with only taxonomies and their abundances.
    """

    # calculate the width of bands for the horizon plot assuming 6 bands
    horizon_scale = (df['abundance'].max() - df['abundance'].min()) // 6

    # scale the abundances so that 0 is the origin for a better y axis label
    def shift(x):
        return x - horizon_scale * 3

    df['abundance'] = df['abundance'].apply(shift)

    # target df
    csv_df = pd.DataFrame()

    # create column in target df for each unique taxonomy
    taxonomies = df['taxonomy'].unique()

    # extract a taxonomy's abundances as series and
    # append it to the target df as new column
    for t in taxonomies:
        f = df['taxonomy'] == t
        csv_df[t] = pd.Series(df[f]['abundance']).reset_index(drop=True)

    # write file
    csv_df.to_csv(_CSV_FILE, index=False)

    return horizon_scale


def _export_env_vars(width: int, height: int, horizon_scale: int) -> None:
    """
    Export necessary variables to environment
    so that the R script can read them from the environment.
    """

    environ['HORI_WIDTH'] = str(width)
    environ['HORI_HEIGHT'] = str(height)
    environ['HORI_BANDWIDTH'] = str(horizon_scale)
    environ['HORI_CSV'] = _CSV_FILE
    environ['HORI_IMAGE'] = _IMAGE_FILE


def _call_r_script() -> None:
    """
    Execute the R script that generates a horizon plot.
    """

    cwd = getcwd()
    chdir(r_path)
    # utils = rpackages.importr('utils')
    # make sure lattice is installed

    #utils.install_packages('LatticeExtra', repos="https://cloud.r-project.org")

    call(['Rscript', '--no-save', _HORIZON_R_SCRIPT])
    chdir(cwd)
