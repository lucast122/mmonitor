import sqlite3
from json import loads, dumps
from typing import List, Tuple, Any

import pandas as pd


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

    def create_db(self, db_name):
        """
        Create a new sqlite MMonitor database on the local disk
        """
        con = sqlite3.connect(db_name)
        cursor = con.cursor()

        # mmonitor taxonomies and abundances
        create_command = f"""CREATE TABLE IF NOT EXISTS mmonitor (
            read_id INTEGER PRIMARY KEY,
            taxonomy TEXT,
            abundance INTEGER,
            sample_id INTEGER,
            project_id INTEGER,
            sample_date TEXT
        )"""

        # sample metadata
        cursor.execute(create_command)
        create_command = f"""CREATE TABLE IF NOT EXISTS metadata (
            "sample_id"	INTEGER PRIMARY KEY,
            "data" TEXT
        )"""

        cursor.execute(create_command)
        con.commit()
        con.close()
        self._db_path = db_name

    def update_table_with_emu_out(self, emu_out_path: str, tax_rank: str, sample_name: str,
                                  project_name: str, sample_date):
        """
        Update MMonitor data from a file containing emu.py output
        """
        con = sqlite3.connect(self._db_path)
        cursor = con.cursor()

        df = pd.read_csv(
            f"{emu_out_path}/{sample_name}_rel-abundance.tsv",
            sep='\t',
            header=None,
            usecols=[0, 1, 2, 3],
            names=['Taxid', 'Abundance', 'Species','Genus']
        )
        df = df.sort_values('Abundance', ascending=False)
        #drop first row
        df = df.iloc[1:]
        # format name

        # df['Name'] = df['Name'].apply(lambda s: s.strip())
        # add sample name
        # emu.py uses different columns for species and genus, so combine them to get the full species name
        #full_name = f"{df['species']} {df['genus']}"
        df['Sample'] = sample_name
        df['Sample_date'] = sample_date
        #df = df[df['Rank'] == tax_rank]
        #df = df.drop(columns='Rank')
        for index, row in df.iterrows():
            #add all species to the database with their abundance percentages don't check count when doing 16s as it only outputs percentages
            insert_query = f"""INSERT INTO mmonitor            (taxonomy, abundance, sample_id, project_id) 
            VALUES ('{row['Species']}', {row['Abundance']}, '{sample_name}', '{project_name}')"""
            cursor.execute(insert_query)
        con.commit()
        con.close()

    def update_table_with_kraken_out(self, kraken_out_path: str, tax_rank: str, sample_name: str,
                                     project_name: str, sample_date):
        """
        Update MMonitor data from a file containing kraken output
        """
        con = sqlite3.connect(self._db_path)
        cursor = con.cursor()

        df = pd.read_csv(
            kraken_out_path,
            sep='\t',
            header=None,
            usecols=[1, 3, 5],
            names=['Count', 'Rank', 'Name']
        )
        df = df.sort_values('Count', ascending=False)
        # format name

        # df['Name'] = df['Name'].apply(lambda s: s.strip())
        # add sample name
        df['Sample'] = sample_name
        df['Sample_date'] = sample_date
        df = df[df['Rank'] == tax_rank]
        df = df.drop(columns='Rank')
        for index, row in df.iterrows():
            if row['Count'] > 100:
                check_name_exists_in_sample = f"SELECT EXISTS(SELECT 1 FROM mmonitor WHERE taxonomy='{row['Name']}' AND sample_id='{sample_name}');"
                cursor.execute(check_name_exists_in_sample)
                name_exists = cursor.fetchall()[0][0]
                # if a sample already has a taxon name (e.g. when computation is performed on multiple files for
                # the same sample while sequencing) then add the abundance to the value in the database
                if name_exists == 0:
                    insert_query = f"""INSERT INTO mmonitor
                        (taxonomy, abundance, sample_id, project_id, sample_date) 
                        VALUES 
                        ('{row['Name']}', {row['Count']}, '{sample_name}', '{project_name}','{sample_date}')"""
                    cursor.execute(insert_query)
                # if taxon is new to the data base simply insert it into the table
                elif name_exists == 1:
                    update_query = f"UPDATE mmonitor SET abundance = abundance + {row['Count']} WHERE taxonomy = '{row['Name']}' AND sample_id='{sample_name}'"
                    cursor.execute(update_query)

        con.commit()
        con.close()

    def append_metadata_from_csv(self, csv_file: str):
        """
        Update metadata from a csv
        """
        con = sqlite3.connect(self._db_path)
        df = pd.read_csv(csv_file)
        meta_cols = [col for col in df.columns if col != 'sample_id']
        df['data'] = df.apply(lambda x: dumps({col: x[col] for col in meta_cols}), axis=1)
        df = df.drop(columns=meta_cols)
        df.to_sql('metadata', con, if_exists='append', index=False)
        con.commit()
        con.close()
