import base64
import json
import os
import subprocess
from json import loads
from typing import List, Tuple, Any
import pandas as pd
import requests
from requests.auth import HTTPBasicAuth
from Bio import SeqIO
import keyring
from keyring.errors import PasswordDeleteError
from requests.packages.urllib3.exceptions import InsecureRequestWarning
from requests.exceptions import Timeout, RequestException





# Disable the InsecureRequestWarning
requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

def _parse_dict(x):
    return pd.Series(loads(x))


def _explode_metadata(df):
    return pd.concat([df, df['data'].apply(_parse_dict)], axis=1).drop(columns='data')


"""
This class gets the path to a db_config file. The db_config is has the host ip 'host', the username 'user' and the password 'password'
 
"""


def convert_date_format(date_str):
    parts = date_str.split('.')

    # Check if the date is in the format DD.MM.YYYY or DD.MM.YY
    if len(parts) == 3 and all(part.isdigit() for part in parts):
        day, month, year = parts

        # Convert 2-digit year to 4-digit year assuming it's in the 2000s
        if len(year) == 2:
            year = '20' + year

        return f"{year}-{month}-{day}"
    else:
        # Return the original date if the format is not recognized
        return date_str


class DjangoDBInterface:

    def __init__(self, config_path):
        self.config_path = config_path
        self.username = None
        self.password = None
        self.host = None
        self.port = None
        self.offline_mode = False
        self.load_config()

    def set_offline_mode(self, offline):
        self.offline_mode = offline
        if offline:
            self.host = '127.0.0.1'
            self.port = '8000'
        else:
            # Reset to default online values if needed
            self.host = 'mmonitor.org'
            self.port = '443'
        self.save_config()

    def load_config(self):
        if os.path.exists(self.config_path):
            with open(self.config_path, 'r') as file:
                config = json.load(file)
                self.username = config.get('user')
                self.host = config.get('host')
                self.port = config.get('port', 443)
                self.offline_mode = config.get('offline_mode', False)
                
                stored_password = keyring.get_password("MMonitor", self.username)
                if stored_password:
                    self.password = stored_password
                else:
                    self.password = config.get('password')

    def save_config(self, remember=False):
        config = {
            'user': self.username,
            'host': self.host,
            'port': self.port,
            'offline_mode': self.offline_mode
        }
        if remember:
            keyring.set_password("MMonitor", self.username, self.password)
        else:
            try:
                keyring.delete_password("MMonitor", self.username)
            except PasswordDeleteError:
                pass
        
        with open(self.config_path, 'w') as file:
            json.dump(config, file)

    def login(self, username, password, host, port, remember=False):
        self.username = username
        self.password = password
        self.host = host
        self.port = port
        
        # Set offline_mode to False for online login attempts
        self.offline_mode = False if host != '127.0.0.1' else True
        
        print(f"Attempting to login with: {username}@{host}:{port}")
        print(f"Offline mode: {self.offline_mode}")
        
        if self.verify_credentials():
            self.save_config(remember)
            return True
        return False

    def verify_credentials(self):
        user_id = self.get_user_id(self.username, self.password)
        print(f"Received user_id: {user_id}")
        return user_id is not None

    def get_user_id(self, username: str, password: str):
        protocol = 'http' if self.offline_mode else 'https'
        django_url = f"{protocol}://{self.host}:{self.port}/users/get_user_id/"
        try:
            response = requests.post(
                django_url, 
                data={'username': username, 'password': password}, 
                verify=not self.offline_mode,
                timeout=10
            )
            print(f"Response status code: {response.status_code}")
            print(f"Response content: {response.content}")
            
            if response.status_code == 200:
                try:
                    return response.json().get('user_id')
                except json.JSONDecodeError:
                    print(f"Failed to decode JSON. Raw response: {response.text}")
                    return None
            elif response.status_code == 400:
                print(f"Authentication failed: {response.text}")
                return None
            else:
                print(f"Unexpected status code: {response.status_code}")
                return None
        except Timeout:
            print("Request timed out")
            return None
        except RequestException as e:
            print(f"Request exception: {e}")
            return None
        except Exception as e:
            print(f"Exception occurred while getting user ID: {str(e)}")
            return None

    def get_unique_sample_ids(self):
        protocol = 'http' if self.offline_mode else 'https'
        django_url = f"{protocol}://{self.host}:{self.port}/users/get_unique_sample_ids/"
        response = requests.post(
            django_url,
            data={'username': self.username, 'password': self.password},
            verify=not self.offline_mode
        )
        if response.status_code == 200:
            return response.json().get('sample_ids', [])
        else:
            return None

    def query_to_dataframe(self, query: str) -> pd.DataFrame:
        return pd.read_sql_query(query, self.conn)

    def query_to_list(self, query: str) -> List[Tuple[Any]]:
        self.cursor.execute(query)
        return self.cursor.fetchall()


    def update_django_with_emu_out(self, emu_out_path: str, tax_rank: str, sample_name: str, project_name: str,
                                   sample_date: str, subproject_name: str, overwrite: bool):
        user_id = self.get_user_id(self._db_config['user'], self._db_config['password'])

        if user_id is None:
            print("Invalid user credentials")
            return

        try:
            sample_ids = self.get_unique_sample_ids()
            if sample_name in sample_ids and not overwrite:
                print(
                    f"Skipping sample {sample_name} as it is already in the database. Select overwrite to reprocess a sample.")
                return
        except TypeError as e:
            print(f"Type error {e}")

        abundance_file = f"{emu_out_path}/{sample_name}_rel-abundance-threshold.tsv"
        if not os.path.exists(abundance_file):
            abundance_file = f"{emu_out_path}/{sample_name}_rel-abundance.tsv"
            if not os.path.exists(abundance_file):
                print(f"No abundance file found for sample {sample_name}. Skipping database update.")
                return

        try:
            df = pd.read_csv(
                abundance_file,
                sep='\t',
                header=None,
                usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13],
                names=['Taxid', 'Abundance', 'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Superkingdom',
                       'Clade', 'Subspecies', "count"]
            )
        except Exception as e:
            print(f"Error reading abundance file for sample {sample_name}: {e}")
            return

        df.fillna("Not Available", inplace=True)
        df.sort_values('Abundance', ascending=False, inplace=True)
        df = df.iloc[1:]  # Skipping the first row, assuming it's headers or unwanted data
        df['Sample'] = sample_name
        df['Sample'] = sample_name
        sample_date = convert_date_format(sample_date)  # Convert date format if necessary
        df['Sample_date'] = sample_date

        # Prepare a list of records
        records = []
        for index, row in df.iterrows():
            # print(row)
            records.append({
                "taxonomy": row['Species'],
                "tax_genus": row['Genus'],
                "tax_family": row['Family'],
                "tax_order": row['Order'],
                "tax_class": row['Class'],
                "tax_phylum": row['Phylum'],
                "tax_superkingdom": row['Superkingdom'],
                "tax_clade": row['Clade'],
                "tax_subspecies": row['Subspecies'],
                "abundance": row['Abundance'],
                "count": row["count"],
                "sample_id": sample_name,
                "project_id": project_name,
                "user_id": user_id,
                "subproject": subproject_name,
                "date": sample_date
            })

        # Send all records in one request
        try:
            response = requests.post(
                f"https://{self._db_config['host']}:{self.port}/users/overwrite_nanopore_record/",
                json=records,  # Send the list of records
                auth=HTTPBasicAuth(self._db_config['user'], self._db_config['password'])
            )
            if response.status_code != 201:
                print(f"Failed to add records: {response.content}")
                return False
            else:
                print(f"Records added successfully.")
                return True
        except Exception as e:
            print(e)
            return False

    def send_nanopore_record_centrifuge(self, kraken_out_path: str, sample_name: str, project_id: str, subproject_id: str,
                                        date: str, overwrite: bool):
        """Send Centrifuger results to the database"""
        # First verify user credentials
        user_id = self.get_user_id(self.username, self.password)
        if user_id is None:
            print("Invalid user credentials")
            return False

        try:
            # Read and process the Centrifuger output
            df = pd.read_csv(
                kraken_out_path,
                sep='\t',
                header=None,
                usecols=[0, 1, 3, 5],
                names=["abundance", 'Count', 'Rank', 'Name']
            )

            df = df.sort_values('Count', ascending=False)
            df['Sample'] = sample_name
            df['Sample_date'] = date
            df = df[df['Rank'] == "S"]  # Keep only species-level classifications
            df = df[df['abundance'] > 1]  # Filter by abundance
            df = df.drop(columns='Rank')

            print(f"User id clientside: {user_id}")

            # Check if sample already exists
            sample_ids = self.get_unique_sample_ids()
            print(f"Found samples: {sample_ids}")
            
            if sample_name in sample_ids and not overwrite:
                print(f"Skipping sample {sample_name} as it is already in the database. Select overwrite to reprocess a sample.")
                return False

            # Prepare records for database
            records = []
            for index, row in df.iterrows():
                record_data = {
                    "sample_id": sample_name,
                    "project_id": project_id,
                    "subproject_id": subproject_id,
                    "date": date,
                    "taxonomy": row["Name"].strip(),
                    "abundance": row["abundance"]/100,  # Convert to same format as EMU
                    "count": row["Count"],
                    "project_id": project_id,
                    "subproject": subproject_id,
                    "tax_genus": "empty",
                    "tax_family": "empty",
                    "tax_order": "empty",
                    "tax_class": "empty",
                    "tax_phylum": "empty",
                    "tax_superkingdom": "empty",
                    "tax_clade": "empty",
                    "tax_subspecies": "empty"
                }
                records.append(record_data)

            # Send to database
            protocol = 'http' if self.offline_mode else 'https'
            url = f"{protocol}://{self.host}:{self.port}/users/overwrite_nanopore_record/"
            
            try:
                response = requests.post(
                    url,
                    json=records,
                    auth=HTTPBasicAuth(self.username, self.password),
                    verify=not self.offline_mode
                )
                
                if response.status_code == 201:
                    print(f"Records added successfully for sample {sample_name}")
                    return True
                else:
                    print(f"Failed to add records: {response.content}")
                    return False
                    
            except Exception as e:
                print(f"Error sending data to server: {e}")
                return False
                
        except Exception as e:
            print(f"Error processing Centrifuger output: {e}")
            return False

    import os


    def upload_mag(self, name: str, taxonomy: str, sample_name: str, gff_file_path: str, fasta_file_path: str):
        user_id = self.get_user_id(self._db_config['user'], self._db_config['password'])
        if user_id is None:
            print("Invalid user credentials")
            return

        # Generate the .fai file using samtools faidx
        try:
            subprocess.run(['samtools', 'faidx', fasta_file_path], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error generating .fai file: {e}")
            return

        fai_file_path = f"{fasta_file_path}.fai"

        url = f"https://{self._db_config['host']}:{self.port}/users/upload_mag/"
        auth = HTTPBasicAuth(self._db_config['user'], self._db_config['password'])
        
        # Create the payload for metadata and files
        files = {
            'gff_file': open(gff_file_path, 'rb'),
            'fasta_file': open(fasta_file_path, 'rb'),
            'fai_file': open(fai_file_path, 'rb')
            
        }
        data = {
            "name": name,
            "taxonomy": taxonomy,
            "sample_name": sample_name,
        }

        try:
            response = requests.post(url, files=files, data=data, auth=auth)
            response.raise_for_status()  # Raise an exception for bad status codes
            if response.status_code == 201:
                print(f"MAG {name} uploaded successfully.")
                return response.json()
            else:
                print(f"Failed to upload MAG {name}: {response.content}")
                return None
        except requests.exceptions.RequestException as e:
            print(f"Error uploading MAG {name}: {e}")
            return None
        finally:
            # Close file handlers
            files['gff_file'].close()
            files['fasta_file'].close()
            files['fai_file'].close()









