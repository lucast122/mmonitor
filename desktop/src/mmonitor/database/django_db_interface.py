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
import threading
import traceback
import datetime
import numpy as np

# Get the src directory path
SRC_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

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
        
        print("\n=== Initializing DjangoDBInterface ===")
        
        # Auto-detect offline mode based on config path
        if config_path and ("pipeline_config" in config_path):
            print("Pipeline config detected - setting offline mode by default")
            self.offline_mode = True
            self.username = "offlinemode"
            self.password = "offline123"
        
        # Load configuration from file
        self.load_config()
        
        # Set base URL based on mode
        self.update_base_url()
        
        # Debug info after initialization
        print(f"Initialized with: host={self.host}, port={self.port}")
        print(f"Offline mode: {self.offline_mode}")
        print(f"Username: {self.username}")
        print(f"Base URL: {self.base_url}")

    def update_base_url(self):
        """Update base URL based on offline mode"""
        if self.offline_mode:
            self.host = "127.0.0.1"
            self.port = "8000"
            self.base_url = f"http://{self.host}:{self.port}"
        else:
            # Always use HTTP for localhost or 127.0.0.1 connections
            if self.host == "127.0.0.1" or self.host == "localhost":
                protocol = "http"
            else:
                protocol = "https"
                
            self.base_url = f"{protocol}://{self.host}"
            if self.port:
                self.base_url += f":{self.port}"
        
        print(f"Base URL set to: {self.base_url}")

    def set_offline_mode(self, offline=True):
        """Set offline mode and update URLs accordingly"""
        self.offline_mode = offline
        if offline:
            # Always set these values when in offline mode
            self.username = "offlinemode"
            self.password = "offline123"
            self.host = "127.0.0.1"
            self.port = "8000"
        self.update_base_url()
        print(f"Offline mode set to: {offline}, base URL: {self.base_url}")
        print(f"Using credentials: username='{self.username}', host='{self.host}', port='{self.port}'")

    def load_config(self):
        """Load configuration from file with proper error handling"""
        print(f"Loading config from: {self.config_path}")
        
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
                        
                        # Ensure offline mode credentials are properly set if in offline mode
                        if self.offline_mode:
                            self.username = "offlinemode"
                            self.password = "offline123"
                            self.host = "127.0.0.1"
                            self.port = "8000"
                        
                        # Update base URL after loading config
                        self.update_base_url()
                        
                        print(f"Loaded config: host={self.host}, port={self.port}, offline_mode={self.offline_mode}")
                        print(f"Username={self.username}, base_url={self.base_url}")
                        
                    except json.JSONDecodeError as e:
                        print(f"Error parsing JSON config: {e}")
                        print(f"Content causing error: {content}")
                        
            except Exception as e:
                print(f"Error reading config file: {e}")
        else:
            print(f"Config file not found: {self.config_path}")

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
        
        # Only set offline mode if explicitly logging into localhost with specific credentials
        if host == '127.0.0.1' and username == 'offlinemode':
            self.offline_mode = True
        else:
            # For all other cases, respect user's intent to login normally
            self.offline_mode = False
        
        # Update base_url with the new host and port
        self.update_base_url()
        
        print(f"Attempting to login with: {username}@{host}:{port}")
        print(f"Offline mode: {self.offline_mode}")
        
        if self.verify_credentials():
            self.save_config(remember)
            return True
        return False

    def verify_credentials(self):
        """Verify credentials with proper offline mode handling"""
        # Only force offline mode if explicit offline mode was chosen
        if self.offline_mode and ("127.0.0.1" in self.base_url or "localhost" in self.base_url):
            self.username = "offlinemode"
            self.password = "offline123"
            print("Using offline mode credentials")
            user_id = self.get_user_id("offlinemode", "offline123")
        else:
            # Ensure username is not None
            if self.username is None:
                print("ERROR: Username is None in verify_credentials")
                return False
                
            # In online mode, use provided credentials
            user_id = self.get_user_id(self.username, self.password)
        
        print(f"Verifying credentials: offline_mode={self.offline_mode}, user_id={user_id}")
        return user_id is not None

    def get_user_id(self, username: str, password: str):
        """Get user ID with proper offline mode handling"""
        # Use base_url which already has the correct protocol
        django_url = f"{self.base_url}/users/get_user_id/"
        
        # Safety check - never use None for username/password
        if username is None or password is None:
            print("ERROR: Username or password is None in get_user_id")
            return None
        
        print(f"Attempting authentication with {username} at {django_url}")
        
        try:
            # Determine SSL verification based on host
            verify_ssl = not ("127.0.0.1" in self.base_url or "localhost" in self.base_url)
            
            response = requests.post(
                django_url, 
                data={'username': username, 'password': password}, 
                verify=verify_ssl,
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
            # Force offline mode if using localhost URL
            if "127.0.0.1" in self.base_url or "localhost" in self.base_url:
                self.offline_mode = True
                self.username = "offlinemode"
                self.password = "offline123"
                print("Detected localhost URL - forcing offline mode for sample IDs")
            
            # Ensure we have valid credentials
            username = "offlinemode" if self.offline_mode else self.username
            password = "offline123" if self.offline_mode else self.password
            
            # Final safety check - never use None username
            if username is None:
                print("WARNING: Username is None, forcing offlinemode credentials for sample IDs")
                username = "offlinemode"
                password = "offline123"
                self.offline_mode = True
            
            url = f"{self.base_url}/users/get_unique_sample_ids/"
            auth = HTTPBasicAuth(username, password)
            
            print(f"Getting sample IDs as user: {username}")
            
            # Determine SSL verification based on host
            verify_ssl = not (self.offline_mode or "127.0.0.1" in self.base_url or "localhost" in self.base_url)
            
            response = requests.post(
                url,
                auth=auth,
                verify=verify_ssl
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
        
        # Force offline mode check from current URL
        if "127.0.0.1" in self.base_url or "localhost" in self.base_url:
            self.offline_mode = True
            self.username = "offlinemode"
            self.password = "offline123"
            print("Detected localhost URL - forcing offline mode")
        
        # Double check offline mode credentials are properly set
        if self.offline_mode:
            self.username = "offlinemode"
            self.password = "offline123"
            
        # Force credentials if username is still None
        if self.username is None:
            if self.offline_mode or "127.0.0.1" in self.base_url or "localhost" in self.base_url:
                self.username = "offlinemode"
                self.password = "offline123"
                self.offline_mode = True
                print("Setting offline credentials due to None username")
            else:
                print("WARNING: Username is None but not in offline mode")
        
        # Use base_url which already has the correct protocol
        url = f"{self.base_url}/users/overwrite_nanopore_record/"
        
        print(f"Uploading as user: {self.username}")
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
                    # Force offline mode for localhost URLs
                    if "127.0.0.1" in self.base_url or "localhost" in self.base_url:
                        use_offline = True
                    else:
                        use_offline = self.offline_mode
                    
                    # Ensure we're using the correct credentials for the current mode
                    username = "offlinemode" if use_offline else self.username
                    password = "offline123" if use_offline else self.password
                    
                    # Final safety check - never use None username
                    if username is None:
                        print("WARNING: Username is still None, forcing offlinemode credentials")
                        username = "offlinemode"
                        password = "offline123"
                    
                    # Log which user is being used for the upload
                    print(f"Thread uploading data as user: {username}")
                    
                    # Create auth with appropriate credentials
                    auth = HTTPBasicAuth(username, password)
                    
                    # Determine SSL verification based on host
                    verify_ssl = not (use_offline or "127.0.0.1" in self.base_url or "localhost" in self.base_url)
                    
                    # Send the request
                    response = requests.post(
                        url,
                        json=records,
                        auth=auth,
                        verify=verify_ssl,
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
                    traceback.print_exc()
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
        print("\n=== Debug Centrifuger Upload Process ===")
        print(f"Base URL: {self.base_url}")
        print(f"Centrifuger output path: {kraken_out_path}")
        
        # Force offline mode check from current URL
        if "127.0.0.1" in self.base_url or "localhost" in self.base_url:
            self.offline_mode = True
            self.username = "offlinemode"
            self.password = "offline123"
            print("Detected localhost URL - forcing offline mode")
        
        # Double check offline mode credentials are properly set
        if self.offline_mode:
            self.username = "offlinemode"
            self.password = "offline123"
            
        # Force credentials if username is still None
        if self.username is None:
            if self.offline_mode or "127.0.0.1" in self.base_url or "localhost" in self.base_url:
                self.username = "offlinemode"
                self.password = "offline123"
                self.offline_mode = True
                print("Setting offline credentials due to None username")
            else:
                print("WARNING: Username is None but not in offline mode")
        
        print(f"Uploading as user: {self.username}")
        print(f"Has password: {bool(self.password)}")
        print(f"Offline mode: {self.offline_mode}")
        
        try:
            # Use base_url which already has the correct protocol
            url = f"{self.base_url}/users/overwrite_nanopore_record/"
            print(f"Using URL: {url}")
            
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
                    # Force offline mode for localhost URLs
                    if "127.0.0.1" in self.base_url or "localhost" in self.base_url:
                        use_offline = True
                    else:
                        use_offline = self.offline_mode
                    
                    # Ensure we're using the correct credentials for the current mode
                    username = "offlinemode" if use_offline else self.username
                    password = "offline123" if use_offline else self.password
                    
                    # Final safety check - never use None username
                    if username is None:
                        print("WARNING: Username is still None, forcing offlinemode credentials")
                        username = "offlinemode"
                        password = "offline123"
                    
                    # Log which user is being used for the upload
                    print(f"Thread uploading data as user: {username}")
                    
                    # Create auth with appropriate credentials
                    auth = HTTPBasicAuth(username, password)
                    
                    # Determine SSL verification based on host
                    verify_ssl = not (use_offline or "127.0.0.1" in self.base_url or "localhost" in self.base_url)
                    
                    # Send the request
                    response = requests.post(
                        url,
                        json=records,
                        auth=auth,
                        verify=verify_ssl,
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
                    traceback.print_exc()
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

    def upload_sequencing_statistics(self, data, sample_name):
        """Upload sequencing statistics to database"""
        # Only ensure proper credentials if explicitly in offline mode
        if self.offline_mode:
            self.username = "offlinemode"
            self.password = "offline123"
            
        # Handle None username case gracefully
        if self.username is None:
            print("ERROR: Username is None in upload_sequencing_statistics")
            return None
                
        print(f"Uploading as user: {self.username}")
        print(f"Has password: {bool(self.password)}")
        print(f"Offline mode: {self.offline_mode}")
        
        try:
            # Use base_url which already has the correct protocol
            url = f"{self.base_url}/users/sequencing_statistics/"
            print(f"Using URL: {url}")
            
            # Clean and validate the data
            # Convert numpy values to Python native types to ensure proper JSON serialization
            clean_data = {}
            for key, value in data.items():
                if isinstance(value, (np.integer, np.floating, np.ndarray)):
                    if hasattr(value, 'item'):
                        clean_data[key] = value.item()  # Convert numpy scalar to Python scalar
                    elif hasattr(value, 'tolist'):
                        clean_data[key] = value.tolist()  # Convert numpy array to list
                else:
                    clean_data[key] = value
            
            # Send to database in a separate thread
            def upload_thread():
                try:
                    # Force offline mode for localhost URLs
                    if "127.0.0.1" in self.base_url or "localhost" in self.base_url:
                        use_offline = True
                    else:
                        use_offline = self.offline_mode
                    
                    # Ensure we're using the correct credentials for the current mode
                    username = "offlinemode" if use_offline else self.username
                    password = "offline123" if use_offline else self.password
                    
                    # Final safety check - never use None username
                    if username is None:
                        print("WARNING: Username is still None, forcing offlinemode credentials")
                        username = "offlinemode"
                        password = "offline123"
                    
                    # Log which user is being used for the upload
                    print(f"Thread uploading statistics as user: {username}")
                    
                    # Create auth with appropriate credentials
                    auth = HTTPBasicAuth(username, password)
                    
                    # Determine SSL verification based on host
                    verify_ssl = not (use_offline or "127.0.0.1" in self.base_url or "localhost" in self.base_url)
                    
                    # Send the request
                    response = requests.post(
                        url,
                        json=clean_data,
                        auth=auth,
                        verify=verify_ssl,
                        timeout=30
                    )
                    
                    print(f"Statistics upload response status: {response.status_code}")
                    if response.status_code in [200, 201]:
                        print("Successfully uploaded sequencing statistics")
                        return True
                    else:
                        print(f"Statistics upload failed with status {response.status_code}")
                        print(f"Response content: {response.text}")
                        return False
                        
                except Exception as e:
                    print(f"Error sending statistics to server: {e}")
                    traceback.print_exc()
                    return False
            
            # Start upload thread
            upload_thread = threading.Thread(target=upload_thread)
            upload_thread.start()
            upload_thread.join(timeout=60)  # Wait up to 60 seconds
            
            if upload_thread.is_alive():
                print("Statistics upload timed out")
                return False
                
            return True
                
        except Exception as e:
            print(f"Error processing sequencing statistics: {e}")
            traceback.print_exc()
            return False

    def upload_mag(self, fasta_file_path, gff_file_path, name, taxonomy, sample_name):
        """Upload MAG to the database with proper credentials handling"""
        print(f"\n=== Uploading MAG {name} ===")
        print(f"Base URL: {self.base_url}")
        print(f"Sample name: {sample_name}")
        
        # Only ensure proper credentials if explicitly in offline mode
        if self.offline_mode:
            self.username = "offlinemode"
            self.password = "offline123"
            
        # Handle None username gracefully
        if self.username is None:
            print("ERROR: Username is None in upload_mag")
            return None
            
        print(f"Uploading as user: {self.username}")
        print(f"Has password: {bool(self.password)}")
        print(f"Offline mode: {self.offline_mode}")
        
        # Use base_url which already has the correct protocol
        url = f"{self.base_url}/users/upload_mag/"
        print(f"Using URL: {url}")
        
        # Get user ID - use offline credentials only if in offline mode
        username = "offlinemode" if self.offline_mode else self.username
        password = "offline123" if self.offline_mode else self.password
        
        user_id = self.get_user_id(username, password)
        
        if user_id is None:
            print("Invalid user credentials")
            return None

        # Generate the .fai file using samtools faidx
        try:
            subprocess.run(['samtools', 'faidx', fasta_file_path], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error generating .fai file: {e}")
            return None

        fai_file_path = f"{fasta_file_path}.fai"
        
        # Prepare auth
        auth = HTTPBasicAuth(
            "offlinemode" if self.offline_mode else self.username,
            "offline123" if self.offline_mode else self.password
        )
        
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
            # Determine if SSL verification should be used
            verify_ssl = not (self.offline_mode or "127.0.0.1" in self.base_url or "localhost" in self.base_url)
            
            response = requests.post(
                url, 
                files=files, 
                data=data, 
                auth=auth,
                verify=verify_ssl
            )
            
            print(f"MAG upload response status: {response.status_code}")
            
            if response.status_code in [200, 201]:
                print(f"MAG {name} uploaded successfully.")
                return response.json()
            else:
                print(f"Failed to upload MAG {name}: {response.content}")
                return None
                
        except requests.exceptions.RequestException as e:
            print(f"Error uploading MAG {name}: {e}")
            traceback.print_exc()
            return None
            
        finally:
            # Close file handlers
            for f in files.values():
                f.close()

    def upload_nanopore_abundance(self, kraken_out_path: str, sample_name: str, date: str = None, overwrite: bool = False):
        """Upload Nanopore read abundance from the kraken output file
        
        Args:
            kraken_out_path: Path to the Centrifuge/Kraken output file
            sample_name: Unique sample name
            date: Date the sample was collected
            overwrite: Whether to overwrite existing data
            
        Returns:
            bool: True if successful, False otherwise
        """
        print("\n=== Uploading Nanopore Abundance ===")
        print(f"File: {kraken_out_path}")
        print(f"Sample name: {sample_name}")
        print(f"Base URL: {self.base_url}")
        
        # Only ensure proper credentials if explicitly in offline mode
        if self.offline_mode:
            self.username = "offlinemode"
            self.password = "offline123"
            
        # Handle None username gracefully
        if self.username is None:
            print("ERROR: Username is None in upload_nanopore_abundance")
            return None
            
        print(f"Uploading as user: {self.username}")
        print(f"Has password: {bool(self.password)}")
        print(f"Offline mode: {self.offline_mode}")
        
        try:
            # Use base_url which already has the correct protocol
            url = f"{self.base_url}/users/overwrite_nanopore_record/"
            print(f"Using URL: {url}")
            
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
                    "date": date_str,  # Use string date
                    "taxonomy": row["Name"].strip(),
                    "abundance": float(row["abundance"])/100,  # Convert to same format as EMU
                    "count": int(row["Count"]),
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
                    # Force offline mode for localhost URLs
                    if "127.0.0.1" in self.base_url or "localhost" in self.base_url:
                        use_offline = True
                    else:
                        use_offline = self.offline_mode
                    
                    # Ensure we're using the correct credentials for the current mode
                    username = "offlinemode" if use_offline else self.username
                    password = "offline123" if use_offline else self.password
                    
                    # Final safety check - never use None username
                    if username is None:
                        print("WARNING: Username is still None, forcing offlinemode credentials")
                        username = "offlinemode"
                        password = "offline123"
                    
                    # Log which user is being used for the upload
                    print(f"Thread uploading data as user: {username}")
                    
                    # Create auth with appropriate credentials
                    auth = HTTPBasicAuth(username, password)
                    
                    # Determine SSL verification based on host
                    verify_ssl = not (use_offline or "127.0.0.1" in self.base_url or "localhost" in self.base_url)
                    
                    # Send the request
                    response = requests.post(
                        url,
                        json=records,
                        auth=auth,
                        verify=verify_ssl,
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
                    traceback.print_exc()
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
            print(f"Error processing Nanopore abundance: {e}")
            traceback.print_exc()
            return False

    def upload_emu_abundance(self, emu_out_path: str, sample_name: str, date: str = None, overwrite: bool = False):
        """Upload EMU abundance data to the database
        
        Args:
            emu_out_path: Path to the EMU output directory
            sample_name: Unique sample name
            date: Date the sample was collected
            overwrite: Whether to overwrite existing data
            
        Returns:
            bool: True if successful, False otherwise
        """
        print("\n=== Uploading EMU Abundance ===")
        print(f"Directory: {emu_out_path}")
        print(f"Sample name: {sample_name}")
        print(f"Base URL: {self.base_url}")
        
        # Only ensure proper credentials if explicitly in offline mode
        if self.offline_mode:
            self.username = "offlinemode"
            self.password = "offline123"
            
        # Handle None username gracefully
        if self.username is None:
            print("ERROR: Username is None in upload_emu_abundance")
            return None
            
        print(f"Uploading as user: {self.username}")
        print(f"Has password: {bool(self.password)}")
        print(f"Offline mode: {self.offline_mode}")
        
        # Use base_url which already has the correct protocol
        url = f"{self.base_url}/users/overwrite_nanopore_record/"
        
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
            if isinstance(date, (datetime.date, datetime.datetime)):
                date_str = date.strftime('%Y-%m-%d')
            else:
                date_str = str(date)
            
            # Prepare records for database
            records = []
            for _, row in df.iterrows():
                record_data = {
                    "sample_id": sample_name,
                    "date": date_str,
                    "taxonomy": row['species'],
                    "abundance": float(row['abundance']),
                    "count": 1,  # EMU doesn't provide count
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
                    # Force offline mode for localhost URLs
                    if "127.0.0.1" in self.base_url or "localhost" in self.base_url:
                        use_offline = True
                    else:
                        use_offline = self.offline_mode
                    
                    # Ensure we're using the correct credentials for the current mode
                    username = "offlinemode" if use_offline else self.username
                    password = "offline123" if use_offline else self.password
                    
                    # Final safety check - never use None username
                    if username is None:
                        print("WARNING: Username is still None, forcing offlinemode credentials")
                        username = "offlinemode"
                        password = "offline123"
                    
                    # Log which user is being used for the upload
                    print(f"Thread uploading data as user: {username}")
                    
                    # Create auth with appropriate credentials
                    auth = HTTPBasicAuth(username, password)
                    
                    # Determine SSL verification based on host
                    verify_ssl = not (use_offline or "127.0.0.1" in self.base_url or "localhost" in self.base_url)
                    
                    # Send the request
                    response = requests.post(
                        url,
                        json=records,
                        auth=auth,
                        verify=verify_ssl,
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
                    traceback.print_exc()
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
            print(f"Error processing EMU abundance: {e}")
            traceback.print_exc()
            return False
