import sqlite3
from typing import List, Tuple, Any

import pandas as pd
from json import loads


def _parse_dict(x):
    return pd.Series(loads(x))


# https://www.skytowner.com/explore/splitting_dictionary_into_separate_columns_in_pandas_dataframe
def _explode_metadata(df):
    return pd.concat([df, df['data'].apply(_parse_dict)], axis=1).drop(columns='data')


class MMonitorDBInterface:
    """
    Interface to an sqlite MMonitor database.

    Future:
        - Consider switching from raw queries to an ORM like SQAlchemy.
    """

    def __init__(self, db_path: str):
        self._db_path = db_path

    def query_to_dataframe(self, query: str) -> pd.DataFrame:
        """
        Query the database and receive the result as pd Dataframe
        """
        con = sqlite3.connect(self._db_path)
        df = pd.read_sql_query(query, con)
        con.close()
        return df

    def query_to_list(self, query: str) -> List[Tuple[Any]]:
        """
        Query the database and receive the result as list
        """
        con = sqlite3.connect(self._db_path)
        ls = list(con.execute(query))
        con.close()
        return ls

    def get_abundance_meta_by_taxonomy(self, taxonomy: str) -> pd.DataFrame:
        """
        Fetch a join of taxonomy and metadata by sample id as pd Dataframe
        """
        q = "SELECT mmonitor.sample_id, mmonitor.abundance, metadata.* " \
            "FROM mmonitor " \
            "INNER JOIN metadata " \
            "WHERE mmonitor.sample_id = metadata.sample_id " \
            f"AND mmonitor.taxonomy = '{taxonomy}' " \
            "ORDER BY mmonitor.sample_id"
        return _explode_metadata(self.query_to_dataframe(q))

    def get_abundance_by_taxonomy(self, taxonomy: str) -> pd.DataFrame:
        """
        Fetch sample ids and abundances of a taxonomy as pd Dataframe
        """
        q = f"SELECT sample_id, abundance FROM mmonitor WHERE taxonomy = '{taxonomy}' ORDER BY sample_id"
        return self.query_to_dataframe(q)

    def get_all_meta(self) -> pd.DataFrame:
        """
        Fetch all metadata of a sample id as pd Dataframe
        """
        q = "SELECT * FROM metadata ORDER BY sample_id"
        return _explode_metadata(self.query_to_dataframe(q))

    def get_unique_taxonomies(self) -> List[str]:
        """
        Fetch all unique taxonomies as list
        """
        q = "SELECT DISTINCT taxonomy FROM mmonitor"
        return [t[0] for t in self.query_to_list(q)]

    def get_unique_samples(self) -> List[str]:
        """
        Fetch all unique sample ids as list
        """
        q = "SELECT DISTINCT sample_id FROM mmonitor"
        return [t[0] for t in self.query_to_list(q)]
