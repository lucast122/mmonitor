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
from requests.packages.urllib3.exceptions import InsecureRequestWarning # type: ignore
from requests.exceptions import Timeout, RequestException
from build_mmonitor_pyinstaller import ROOT
import threading
import traceback
import datetime




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
        self.host = "127.0.0.1"  # Default to localhost
        self.port = "8000"  # Default port
        self.offline_mode = False
        self.base_url = None
        self.load_config()
        
        # Set base URL based on mode
        self.update_base_url()

    def update_base_url(self):
        """Update base URL based on offline mode"""
        if self.offline_mode:
            self.host = "127.0.0.1"
            self.port = "8000"
            self.base_url = f"http://{self.host}:{self.port}"
        else:
            protocol = "https"
            self.base_url = f"{protocol}://{self.host}"
            if self.port:
                self.base_url += f":{self.port}"

    def set_offline_mode(self, offline=True):
        """Set offline mode and update URLs accordingly"""
        self.offline_mode = offline
        if offline:
            self.username = "offlinemode"
            self.password = "offline123"
            self.host = "127.0.0.1"
            self.port = "8000"
        self.update_base_url()

    def load_config(self):
        """Load configuration from file with proper error handling"""
        if os.path.exists(self.config_path):
            try:
                with open(self.config_path, 'r') as file:
                    content = file.read().strip()
                    if not content:
                        print("Config file is empty")
                        return
                    
                    try:
                        config = json.loads(content)
                        self.username = config.get('user')
                        self.host = config.get('host', self.host)
                        self.port = config.get('port', self.port)
                        self.offline_mode = config.get('offline_mode', False)
                        
                        # Update base URL after loading config
                        self.update_base_url()
                        
                        print(f"Loaded config: host={self.host}, port={self.port}, offline_mode={self.offline_mode}")
                        
                    except json.JSONDecodeError as e:
                        print(f"Error parsing JSON config: {e}")
                        print(f"Content causing error: {content}")
                        
            except Exception as e:
                print(f"Error reading config file: {e}")

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
        """Verify credentials with proper offline mode handling"""
        if self.offline_mode:
            # In offline mode, use offline credentials directly
            user_id = self.get_user_id("offlinemode", "offline123")
        else:
            # In online mode, use provided credentials
            user_id = self.get_user_id(self.username, self.password)
        
        print(f"Verifying credentials: offline_mode={self.offline_mode}, user_id={user_id}")
        return user_id is not None

    def get_user_id(self, username: str, password: str):
        """Get user ID with proper offline mode handling"""
        protocol = 'http' if self.offline_mode else 'https'
        django_url = f"{protocol}://{self.host}:{self.port}/users/get_user_id/"
        
        print(f"Attempting authentication with {username} at {django_url}")
        
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
        """Get list of unique sample IDs from server"""
        try:
            url = f"{self.base_url}/users/get_unique_sample_ids/"
            auth = HTTPBasicAuth(
                "offlinemode" if self.offline_mode else self.username,
                "offline123" if self.offline_mode else self.password
            )
            
            response = requests.post(
                url,
                auth=auth,
                verify=not self.offline_mode
            )
            
            if response.status_code == 200:
                return response.json()
            return None
            
        except Exception as e:
            print(f"Error getting unique sample IDs: {e}")
            return None

    def query_to_dataframe(self, query: str) -> pd.DataFrame:
        return pd.read_sql_query(query, self.conn)

    def query_to_list(self, query: str) -> List[Tuple[Any]]:
        self.cursor.execute(query)
        return self.cursor.fetchall()


    def update_django_with_emu_out(self, emu_out_path: str, tax_rank: str, sample_name: str, project_name: str,
                                   sample_date: str, subproject_name: str, overwrite: bool):
        """Send EMU results to the database with proper error handling"""
        print("\n=== Debug Upload Process ===")
        print(f"Base URL: {self.base_url}")
        print(f"EMU output path: {emu_out_path}")
        
        # Reset credentials if in offline mode
        if self.offline_mode:
            self.username = "offlinemode"
            self.password = "offline123"
            url = "http://127.0.0.1:8000/users/overwrite_nanopore_record/"
        else:
            url = f"{self.base_url}/users/overwrite_nanopore_record/"
        
        print(f"Username: {self.username}")
        print(f"Has password: {bool(self.password)}")
        print(f"Offline mode: {self.offline_mode}")
        print(f"Using URL: {url}")
        
        try:
            # Check for existing samples
            sample_ids = self.get_unique_sample_ids() or []  # Return empty list if None
            if sample_name in sample_ids and not overwrite:
                print(f"Skipping sample {sample_name} as it is already in the database. Select overwrite to reprocess a sample.")
                return False

            # Find the abundance file
            abundance_file = os.path.join(emu_out_path, f"{sample_name}_rel-abundance.tsv")
            if not os.path.exists(abundance_file):
                print(f"No abundance file found at {abundance_file}")
                return False
            
            # First read the header to determine available columns
            with open(abundance_file, 'r') as f:
                header = f.readline().strip().split('\t')
                print(f"Found columns: {header}")
            
            # Read and process the EMU output with proper data types
            df = pd.read_csv(
                abundance_file,
                sep='\t',
                dtype={
                    'abundance': float,  # Ensure abundance is read as float
                    'tax_id': str,
                    'species': str,
                    'genus': str,
                    'family': str,
                    'order': str,
                    'class': str,
                    'phylum': str,
                    'clade': str,
                    'superkingdom': str,
                    'subspecies': str
                }
            )
            print(f"Read {len(df)} records from file")
            
            # Process data
            df.fillna("empty", inplace=True)
            df.sort_values('abundance', ascending=False, inplace=True)
            df = df[df['abundance'] > 0.01]  # Filter by abundance threshold
            print(f"Filtered to {len(df)} records above threshold")
            
            # Convert date to string if it's a datetime object
            if isinstance(sample_date, (datetime.date, datetime.datetime)):
                date_str = sample_date.strftime('%Y-%m-%d')
            else:
                date_str = str(sample_date)
            
            # Prepare records for database
            records = []
            for _, row in df.iterrows():
                record_data = {
                    "sample_id": sample_name,
                    "project_id": project_name,
                    "subproject_id": subproject_name,
                    "date": date_str,
                    "taxonomy": row['species'],
                    "abundance": float(row['abundance']),
                    "count": 1,  # EMU doesn't provide count
                    "project_id": project_name,
                    "subproject": subproject_name,
                    "tax_genus": row['genus'],
                    "tax_family": row['family'],
                    "tax_order": row['order'],
                    "tax_class": row['class'],
                    "tax_phylum": row['phylum'],
                    "tax_superkingdom": row['superkingdom'],
                    "tax_clade": row.get('clade', 'empty'),
                    "tax_subspecies": row.get('subspecies', 'empty')
                }
                records.append(record_data)
                print(f"Prepared record: {record_data}")
            
            print(f"Sending {len(records)} records to server")
            
            # Send to database in a separate thread
            def upload_thread():
                try:
                    # Double check offline mode credentials
                    auth = HTTPBasicAuth(
                        "offlinemode" if self.offline_mode else self.username,
                        "offline123" if self.offline_mode else self.password
                    )
                    
                    # Send the request
                    response = requests.post(
                        url,  # Now url is defined for both offline and online modes
                        json=records,
                        auth=auth,
                        verify=not self.offline_mode,
                        timeout=30
                    )
                    
                    print(f"Server response status: {response.status_code}")
                    print(f"Server response: {response.text}")
                    
                    if response.status_code in [200, 201]:
                        print(f"Successfully uploaded results for sample {sample_name}")
                        return True
                    else:
                        print(f"Upload failed with status {response.status_code}")
                        print(f"Response content: {response.content}")
                        return False
                        
                except Exception as e:
                    print(f"Error sending data to server: {e}")
                    return False
            
            # Start upload thread
            upload_thread = threading.Thread(target=upload_thread)
            upload_thread.start()
            upload_thread.join(timeout=60)  # Wait up to 60 seconds
            
            if upload_thread.is_alive():
                print("Upload timed out")
                return False
                
            return True
                
        except Exception as e:
            print(f"Error processing EMU output: {e}")
            traceback.print_exc()
            return False

    def send_nanopore_record_centrifuger(self, kraken_out_path: str, sample_name: str, project_id: str, subproject_id: str,
                                        date: str, overwrite: bool):
        """Send Centrifuger results to the database"""
        print("\n=== Debug Upload Process ===")
        print(f"Base URL: {self.base_url}")
        print(f"kraken out: {kraken_out_path}")
        
        # Reset credentials if in offline mode
        if self.offline_mode:
            self.username = "offlinemode"
            self.password = "offline123"
        
        print(f"Username: {self.username}")
        print(f"Has password: {bool(self.password)}")
        print(f"Offline mode: {self.offline_mode}")
        
        try:
            # Prepare the URL with trailing slash
            url = f"{self.base_url}/users/overwrite_nanopore_record/"
            
            # Read and process the Centrifuger output
            df = pd.read_csv(
                kraken_out_path,
                sep='\t',
                header=None,
                usecols=[0, 1, 3, 5],
                names=["abundance", 'Count', 'Rank', 'Name']
            )
            print(f"Read {len(df)} records from file")

            df = df.sort_values('Count', ascending=False)
            df['Sample'] = sample_name
            df['Sample_date'] = date
            df = df[df['Rank'] == "S"]  # Keep only species-level classifications
            df = df[df['abundance'] > 0.01]  # Filter by abundance
            df = df.drop(columns='Rank')
            print(f"Filtered to {len(df)} species-level records")

            # Convert date to string if it's a datetime object
            if isinstance(date, (datetime.date, datetime.datetime)):
                date_str = date.strftime('%Y-%m-%d')
            else:
                date_str = str(date)

            # Prepare records for database
            records = []
            for _, row in df.iterrows():
                record_data = {
                    "sample_id": sample_name,
                    "project_id": project_id,
                    "subproject_id": subproject_id,
                    "date": date_str,  # Use string date
                    "taxonomy": row["Name"].strip(),
                    "abundance": float(row["abundance"])/100,  # Convert to same format as EMU
                    "count": int(row["Count"]),
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
                print(f"Prepared record: {record_data}")

            print(f"Sending {len(records)} records to server")
            
            # Send to database in a separate thread
            def upload_thread():
                try:
                    # Double check offline mode credentials
                    auth = HTTPBasicAuth(
                        "offlinemode" if self.offline_mode else self.username,
                        "offline123" if self.offline_mode else self.password
                    )
                    
                    # Send the request - send records as a list, not wrapped in a dict
                    response = requests.post(
                        url,
                        json=records,  # Send records directly as a list
                        auth=auth,
                        verify=not self.offline_mode,
                        timeout=30
                    )
                    
                    print(f"Server response status: {response.status_code}")
                    print(f"Server response: {response.text}")
                    
                    if response.status_code in [200, 201]:
                        print(f"Successfully uploaded results for sample {sample_name}")
                        return True
                    else:
                        print(f"Upload failed with status {response.status_code}")
                        print(f"Response content: {response.content}")
                        return False
                        
                except Exception as e:
                    print(f"Error sending data to server: {e}")
                    return False

            # Start upload thread
            upload_thread = threading.Thread(target=upload_thread)
            upload_thread.start()
            upload_thread.join(timeout=60)  # Wait up to 60 seconds
            
            if upload_thread.is_alive():
                print("Upload timed out")
                return False
                
            return True
                
        except Exception as e:
            print(f"Error processing Centrifuger output: {e}")
            traceback.print_exc()
            return False

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

 






