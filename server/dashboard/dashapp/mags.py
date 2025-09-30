"""MAGs app for visualizing and analyzing Metagenome-Assembled Genomes."""
import dash
from dash import dcc, html, dash_table
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from django_plotly_dash import DjangoDash
from dash_iconify import DashIconify
import pandas as pd
import logging
import plotly.graph_objects as go
import plotly.express as px
from users.models import MAG
import os
import io
import base64
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from datetime import datetime
import traceback
from dash.dash_table.Format import Format, Scheme
import time
import random
from dash.exceptions import PreventUpdate
from dash.dependencies import Input, Output, State, ClientsideFunction, ALL
import json
import math

# Import pyCirclize and BioPython for circular visualization
from pycirclize import Circos
from Bio.SeqFeature import SeqFeature, FeatureLocation

PYCIRCLIZE_AVAILABLE = True
logger = logging.getLogger(__name__)

class MAGs:
    def __init__(self, user_id):
        """Initialize the MAGs app with the given user ID.
        
        Args:
            user_id: The user's ID to filter MAGs
        """
        self.app = DjangoDash('dashboard_mags')
        self.user_id = user_id
        
        # Initialize with default empty values
        self.mag_data = {}  # Initialize as empty dictionary (not DataFrame)
        self.mag_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mag_data")
        self.mag_files = {}
        
        app_name = "mags"
        try:
            # Ensure the mag dir exists
            os.makedirs(self.mag_dir, exist_ok=True)
            
            # Load MAG data
            self._load_mag_data()
            
            # Create layout
            self.app.layout = self.create_layout()
            
        except Exception as e:
            logger.error(f"Error initializing MAGs app: {str(e)}", exc_info=True)
            # Create a minimal app to avoid breaking the whole dashboard
            self.app = DjangoDash(app_name, suppress_callback_exceptions=True)
            self.app.layout = html.Div("Error loading MAGs app")
    
    def _load_mag_data(self):
        """Load MAG data from the database or create mock data if none exists."""
        try:
            logger.info("Loading MAG data")
            # Get all MAG data using the method that fetches from DB or creates mock data
            self.mag_data = self.get_mag_data()
            
            # Check if we have data
            if not self.mag_data:
                logger.warning("No MAG data found, creating mock data")
                self.mag_data = self._create_mock_data_from_table()
                
            # Make sure we treat mag_data consistently as a dictionary, not a DataFrame
            if isinstance(self.mag_data, pd.DataFrame):
                # Convert DataFrame to dict if it somehow was created as DataFrame
                self.mag_data = self.mag_data.to_dict('records')
            
            logger.info(f"Loaded MAG data with {len(self.mag_data)} entries")
            return True
        except Exception as e:
            logger.error(f"Error loading MAG data: {str(e)}", exc_info=True)
            # Ensure we have at least an empty dict
            self.mag_data = {}
            return False
    
    # Add the update_mag_table method
    def update_mag_table(self, quality_value='all'):
        """Update the MAG table based on quality filter.
        
        Args:
            quality_value (str): The quality filter to apply ('HIGH', 'MEDIUM', 'LOW', or 'all')
            
        Returns:
            tuple: (table_data, selected_rows) - filtered table data and empty selection
        """
        try:
            logger.info(f"Updating MAG table with quality filter: {quality_value}")
            
            # Get table data
            table_data = self._create_table_data()
            
            if not table_data:
                logger.warning("No MAG data available")
                return [], []
            
            # Apply quality filter if specified
            if quality_value and quality_value.upper() != 'ALL':
                quality_value = quality_value.lower()
                filtered_data = [
                    row for row in table_data 
                    if row.get('quality', '').lower() == quality_value
                ]
                
                logger.info(f"Filtered table data: {len(filtered_data)} rows (from {len(table_data)} total)")
                return filtered_data, []
            
            logger.info(f"Returning all table data: {len(table_data)} rows")
            return table_data, []
            
        except Exception as e:
            logger.error(f"Error updating MAG table: {str(e)}", exc_info=True)
            return [], []
    
    def get_mag_data(self, include_fasta=False, force_refresh=False):
        """Get MAG data for the dashboard, with optional FASTA inclusion."""
        start_time = time.time()
        
        try:
            print("\n===== Loading MAG data =====")
            
            # Use the cache if available and force_refresh is not True
            cache_key = f'_mag_data_cache_{include_fasta}'
            if hasattr(self, cache_key) and getattr(self, cache_key) and not force_refresh:
                print(f"Using cached MAG data (include_fasta={include_fasta})")
                return getattr(self, cache_key)
            
            try:
                # Get data from Django models
                print("Querying database for MAG data")
                from users.models import MAG
                
                # Get MAGs with optional optimization for fetching FASTA data
                qs = MAG.objects.filter(user_id=self.user_id)
                
                # Include FASTA data optionally
                if include_fasta:
                    print("Including FASTA data in query")
                    qs = qs.only('name', 'completeness', 'contamination', 
                               'taxonomy', 'fasta_file')
                else:
                    # Exclude FASTA data for better performance
                    print("Excluding FASTA data for better performance")
                    qs = qs.defer('fasta_file')
                
                # Execute query
                mags = list(qs)
                print(f"Found {len(mags)} MAGs in database")
                
                # Create result dictionary with both ID and bin_id as keys
                mag_data = {}
                
                # Process MAGs from database
                for mag in mags:
                    try:
                        # Create simple dictionary for easy serialization
                        mag_dict = {
                            'id': mag.id,
                            'name': mag.name,
                            'bin_id': mag.name,  # Use name as bin_id
                            'completeness': float(mag.completeness or 0),
                            'contamination': float(mag.contamination or 0),
                            'size': 0,  # Initialize to 0, will calculate from sequence data if available
                            'gc_content': 0,  # Initialize to 0, will calculate from sequence if available
                            'n_contigs': 0,  # Initialize to 0, will calculate from sequence data if available 
                            'taxonomy': mag.taxonomy or 'Unknown',
                            'quality': getattr(mag, 'quality', 'unknown'),
                        }
                        
                        # Add to dictionary
                        mag_id = str(mag.id)
                        bin_id = mag.name
                        
                        # Include FASTA data if requested
                        if include_fasta and mag.fasta_file:
                            try:
                                print(f"Parsing FASTA for {bin_id}")
                                mag_dict['fasta_file'] = mag.fasta_file
                                fasta_list = self.parse_fasta_for_viewer(mag.fasta_file) 
                                mag_dict['fasta'] = fasta_list
                                
                                # Calculate metrics from sequence data
                                if fasta_list:
                                    mag_dict['gc_content'] = self._calculate_gc_content(fasta_list)
                                    
                                    # Calculate size as sum of contig lengths
                                    total_size = sum(len(seq.get('sequence', '')) for seq in fasta_list)
                                    mag_dict['size'] = total_size
                                    
                                    # Calculate n_contigs as the number of sequences
                                    mag_dict['n_contigs'] = len(fasta_list)
                                
                                print(f"Parsed {len(fasta_list)} sequences for {bin_id}")
                            except Exception as parsing_error:
                                print(f"Error parsing FASTA for {bin_id}: {str(parsing_error)}")
                                mag_dict['fasta'] = []
                        
                        # Add to dictionary with both keys
                        mag_data[mag_id] = mag_dict
                        mag_data[bin_id] = mag_dict  # Allow lookup by bin ID too
                    except Exception as mag_error:
                        print(f"Error processing MAG {mag.name}: {str(mag_error)}")
                        continue
                
                # Check if bin.13 (a common test MAG) is in the data
                if 'bin.13' in mag_data:
                    print("bin.13 is in the processed MAG data dictionary")
                else:
                    print("WARNING: bin.13 is NOT in the processed MAG data dictionary")
                
                setattr(self, cache_key, mag_data)
                
                elapsed = time.time() - start_time
                print(f"Loaded MAG data in {elapsed:.2f} seconds")
                print(f"===== MAG data loading complete =====\n")
                return mag_data
                
            except Exception as db_error:
                print(f"Database error: {str(db_error)}")
                print(traceback.format_exc())
                print("Falling back to mock data")
                mock_data = self._create_mock_data_from_table()
                print(f"Created mock data with {len(mock_data)} entries")
                return mock_data
            
        except Exception as e:
            print(f"Unexpected error in get_mag_data: {str(e)}")
            print(traceback.format_exc())
            return {}
                    
    def _create_mock_data_from_table(self):
        """Create mock MAG data from the table data when database access fails."""
        start_time = time.time()
        
        try:
            logger.info("Creating mock MAG data from table")
            
            # OPTIMIZATION: Use cached table data if available to avoid circular reference
            table_data = []
            
            if hasattr(self, '_table_data_cache') and self._table_data_cache:
                logger.info("Using cached table data for mock creation")
                table_data = self._table_data_cache
            else:
                # If we're already in _create_table_data, this would cause infinite recursion
                # So we need to create minimal table data without calling that method
                
                # Define some basic MAGs to ensure we always have data
                logger.info("Creating minimal table data to avoid recursion")
                table_data = [
                    {
                        'id': '13',
                        'bin_id': 'bin.13',
                        'name': 'bin.13',
                        'completeness': 0.95,
                        'contamination': 0.02,
                        'taxonomy': 'Bacteria; Proteobacteria; Gammaproteobacteria',
                        'size_mb': 4.5,
                        'num_contigs': 50,
                        'quality': 'High'
                    },
                    {
                        'id': '5',
                        'bin_id': 'bin.5',
                        'name': 'bin.5',
                        'completeness': 0.87,
                        'contamination': 0.05,
                        'taxonomy': 'Bacteria; Firmicutes; Bacilli',
                        'size_mb': 3.8,
                        'num_contigs': 26,
                        'quality': 'Medium'
                    },
                    {
                        'id': '1',
                        'bin_id': 'bin.1',
                        'name': 'bin.1',
                        'completeness': 0.99,
                        'contamination': 0.01,
                        'taxonomy': 'Bacteria; Bacteroidetes',
                        'size_mb': 2.3,
                        'num_contigs': 15,
                        'quality': 'High'
                    }
                ]
                
                # Add a few random MAGs for diversity
                for i in range(4, 13):
                    if i == 5:  # already added bin.5
                        continue
                        
                    quality = random.choice(['High', 'Medium', 'Low'])
                    completeness = random.uniform(0.75, 0.99) if quality == 'High' else random.uniform(0.5, 0.9)
                    contamination = random.uniform(0.01, 0.05) if quality == 'High' else random.uniform(0.05, 0.15)
                    
                    table_data.append({
                        'id': str(i),
                        'bin_id': f'bin.{i}',
                        'name': f'bin.{i}',
                        'completeness': completeness,
                        'contamination': contamination,
                        'taxonomy': 'Bacteria',
                        'size_mb': random.uniform(1.5, 5.0),
                        'num_contigs': random.randint(10, 50),
                        'quality': quality
                    })
            
            if not table_data:
                logger.warning("No table data available for mock creation")
                return {}
            
            logger.info(f"Found {len(table_data)} MAGs in table data")
            
            # OPTIMIZATION: Create a copy of the table data for future use
            self._table_data_cache = table_data
            
            # Convert table data to MAG data format efficiently
            mag_data = {}
            
            # Set an overall processing timeout
            processing_timeout = 5.0  # seconds
            processing_start = time.time()
            processed_count = 0
            
            for row in table_data:
                # Check if processing is taking too long
                if time.time() - processing_start > processing_timeout:
                    logger.warning(f"Mock data processing timeout after {processed_count} MAGs")
                    break
                    
                try:
                    mag_id = row.get('id') or row.get('bin_id')
                    if not mag_id:
                        continue
                        
                    # Check if this is a known large MAG
                    bin_id = row.get('bin_id', mag_id)
                    known_large_mags = ['bin.13', 'bin.5']
                    is_large_mag = bin_id in known_large_mags
                    
                    # Create mock MAG entry
                    name = row.get('name', bin_id)
                    size_in_bp = row.get('size_bp', 0)
                    size_in_mb = row.get('size_mb', 0)
                    
                    # Calculate size in bp from mb if available
                    if not size_in_bp and size_in_mb:
                        size_in_bp = int(size_in_mb * 1000000)
                    
                    # Create the entry with base data
                    mag_entry = {
                        'id': mag_id,
                        'name': name,
                        'bin_id': bin_id,
                        'completeness': row.get('completeness', 0),
                        'contamination': row.get('contamination', 0),
                        'taxonomy': row.get('taxonomy', 'Unknown'),
                        'quality': row.get('quality', 'unknown'),
                        'size': size_in_bp,
                        'gc_content': row.get('gc_content', 0) or row.get('gc', 0) or random.uniform(0.3, 0.7),
                        'n_contigs': row.get('n_contigs', 0) or row.get('num_contigs', 0) or random.randint(10, 30),
                        'fasta_file': '',
                        'fasta': []  # Empty for mock
                    }
                    
                    # For large MAGs, add special handling info
                    if is_large_mag:
                        mag_entry['is_large_mag'] = True
                        mag_entry['max_contigs'] = 20 if bin_id == 'bin.13' else 26
                    
                    # Store in the data dict by multiple keys for easier lookup
                    mag_data[mag_id] = mag_entry
                    
                    # Also make available by bin_id if different from id
                    if bin_id != mag_id:
                        mag_data[bin_id] = mag_entry
                    
                    # And by name if different from id and bin_id
                    if name != mag_id and name != bin_id:
                        mag_data[name] = mag_entry
                        
                    processed_count += 1
                    
                except Exception as e:
                    logger.error(f"Error processing mock MAG {row.get('bin_id', 'unknown')}: {str(e)}")
                    continue
                
            mock_mag_count = len(set(entry.get('id') for entry in mag_data.values() if isinstance(entry, dict)))
            logger.info(f"Created mock data with {mock_mag_count} unique MAGs ({len(mag_data)} total entries including duplicates)")
            
            # Track performance
            elapsed = time.time() - start_time
            logger.info(f"_create_mock_data_from_table completed in {elapsed:.4f}s")
            
            if hasattr(self, 'performance_stats'):
                self.performance_stats['create_mock_data'] = {
                    'total_time': elapsed,
                    'table_rows': len(table_data),
                    'mock_mags': mock_mag_count,
                    'total_entries': len(mag_data)
                }
            
            return mag_data
            
        except Exception as e:
            logger.error(f"Error creating mock MAG data: {str(e)}", exc_info=True)
            elapsed = time.time() - start_time
            logger.error(f"Error after {elapsed:.4f}s of processing")
            
            # Emergency fallback - return at least some data
            emergency_data = {
                'bin.13': {
                    'id': '13',
                    'bin_id': 'bin.13',
                    'name': 'bin.13',
                    'completeness': 0.95,
                    'contamination': 0.02,
                    'taxonomy': 'Bacteria; Emergency Fallback',
                    'quality': 'High',
                    'size': 4500000,
                    'gc_content': 0.55,
                    'n_contigs': 50,
                    'fasta_file': '',
                    'fasta': [],
                    'is_large_mag': True,
                    'max_contigs': 20
                }
            }
            # Also add with alternate key
            emergency_data['13'] = emergency_data['bin.13']
            
            return emergency_data

    def _calculate_gc_content(self, sequences):
        """Calculate GC content for a list of sequences."""
        total_gc = 0
        total_length = 0
        
        for seq in sequences:
            sequence = seq.get('sequence', '')
            if sequence:
                gc_count = sequence.upper().count('G') + sequence.upper().count('C')
                total_gc += gc_count
                total_length += len(sequence)
                
        if total_length > 0:
            return (total_gc / total_length)  # Return as a decimal (0-1) for percentage format
        return 0

    def get_gff_data(self, mag_id):
        """Get GFF data for the specified MAG."""
        try:
            logger.info(f"Getting GFF data for MAG {mag_id}")
            
            # Find the MAG in the dictionary
            mag_info = None
            mag_id_str = str(mag_id).strip()
            
            # Try direct lookup by key
            if mag_id_str in self.mag_data:
                mag_info = self.mag_data[mag_id_str]
            else:
                # Try to find by name or bin_id in values
                for key, info in self.mag_data.items():
                    if isinstance(info, dict):
                        if info.get('name') == mag_id_str or info.get('bin_id') == mag_id_str:
                            mag_info = info
                            break
                
                # Try case-insensitive match if still not found
                if not mag_info:
                    mag_id_lower = mag_id_str.lower()
                    for key, info in self.mag_data.items():
                        if isinstance(info, dict):
                            if str(info.get('name', '')).lower() == mag_id_lower or str(info.get('bin_id', '')).lower() == mag_id_lower:
                                mag_info = info
                                break
            
            if not mag_info:
                logger.error(f"No MAG found with ID or name {mag_id}")
                return None
            
            # Get GFF file content
            gff_file = mag_info.get('gff_file', '')
            
            if not isinstance(gff_file, str) or not gff_file or gff_file == 'empty':
                logger.warning(f"No GFF file found for MAG {mag_id}, checking database...")
                
                # Try to get GFF from database
                try:
                    from users.models import MAG
                    mag_obj = MAG.objects.get(name=mag_id, user_id=self.user_id)
                    gff_file = mag_obj.gff_file
                    
                    if not gff_file:
                        logger.error(f"No GFF file found in database for MAG {mag_id}")
                        return None
                        
                    # Update the mag_info with the GFF file for future use
                    mag_info['gff_file'] = gff_file
                    
                except Exception as db_error:
                    logger.error(f"Error retrieving GFF from database: {str(db_error)}")
                    return None
            
            logger.info(f"GFF file found, length: {len(gff_file)} characters")
            
            # Parse GFF data into feature objects
            try:
                # Get FASTA sequences for this MAG to enable contig mapping
                sequences = []
                if 'fasta' in mag_info and isinstance(mag_info['fasta'], list):
                    sequences = mag_info['fasta']
                elif 'fasta_file' in mag_info and mag_info['fasta_file']:
                    sequences = self.parse_fasta_for_viewer(mag_info['fasta_file'])
                
                if not sequences:
                    logger.warning(f"No sequence data found for MAG {mag_id}, GFF parsing may not work properly")
                
                # Parse the GFF content into annotations
                annotations = self._parse_gff_for_annotations(gff_file, sequences)
                logger.info(f"Successfully parsed {len(annotations)} annotations from GFF file")
                
                return annotations
            except Exception as parsing_error:
                logger.error(f"Error parsing GFF file: {str(parsing_error)}", exc_info=True)
                return None
                
        except Exception as e:
            logger.error(f"Error getting GFF data for MAG {mag_id}: {str(e)}", exc_info=True)
            return None

    def _assess_quality(self, row):
        """Assess the quality of a MAG based on completeness and contamination."""
        try:
            # Handle None and convert to float safely
            completeness = float(row.get('completeness', 0) or 0)
            contamination = float(row.get('contamination', 0) or 0)
            
            # Check if values are already in 0-1 range or 0-100 range
            completeness_for_comparison = completeness * 100 if completeness <= 1 else completeness
            contamination_for_comparison = contamination * 100 if contamination <= 1 else contamination
            
            if completeness_for_comparison >= 90 and contamination_for_comparison <= 5:
                return 'high'
            elif completeness_for_comparison >= 70 and contamination_for_comparison <= 10:
                return 'medium'
            else:
                return 'low'
        except Exception as e:
            logger.error(f"Error assessing MAG quality: {str(e)}")
            return 'unknown'

    def _create_initial_figures(self, highlight_mag=None):
        """Create initial figures for the MAGs app."""
        figures = {}
        
        if not self.mag_data:
            logger.warning("No MAG data available for creating figures")
            empty_fig = go.Figure()
            empty_fig.update_layout(
                title="No MAG Data Available",
                annotations=[{
                    'text': "No MAG data found in database",
                    'xref': "paper",
                    'yref': "paper",
                    'showarrow': False,
                    'font': {'size': 20}
                }]
            )
            return {
                'quality_plot': empty_fig,
                'rna_content_plot': empty_fig,
                'kegg_plot': empty_fig
            }
        
        # Quality assessment scatter plot
        figures['quality_plot'] = self._create_quality_scatter(highlight_mag)
        
        # Gene content distribution
        # figures['gene_content_plot'] = self._create_gene_content_plot(highlight_mag)
        
        # RNA counts
        figures['rna_content_plot'] = self._create_rna_counts_plot(highlight_mag)
        
        # KEGG pathway overview
        figures['kegg_plot'] = self._create_kegg_overview(highlight_mag)
        
        return figures
    
    def _create_quality_scatter(self, highlight_mag=None):
        """Create scatter plot of MAG quality."""
        try:
            # Get MAG data
            mag_data = self.get_mag_data()
            if not mag_data:
                logger.warning("No MAG data for quality scatter plot")
                return go.Figure()

            # Create figure
            fig = go.Figure()

            # Organize MAGs by quality category
            mags_by_quality = {
                'High': [],
                'Medium': [],
                'Low': [],
                'Unknown': []
            }
            
            # Process MAG data (skip duplicates by tracking processed names)
            processed_names = set()
            
            for mag_id, mag_info in mag_data.items():
                if not isinstance(mag_info, dict) or 'name' not in mag_info:
                    continue
                    
                # Skip if we've already processed this MAG name
                mag_name = mag_info.get('name', '')
                if mag_name in processed_names:
                    continue
                processed_names.add(mag_name)
                
                # Get MAG quality data
                completeness = float(mag_info.get('completeness', 0) or 0)
                contamination = float(mag_info.get('contamination', 0) or 0)
                taxonomy = mag_info.get('taxonomy', 'Unknown')
                size = float(mag_info.get('size', 1000000))  # Use genome size for bubble size
                
                # Assess quality if not already present
                quality = mag_info.get('quality', '').capitalize()
                if not quality or quality.lower() == 'unknown':
                    quality = self._assess_quality({
                        'completeness': completeness,
                        'contamination': contamination
                    }).capitalize()
                
                # Add to appropriate quality category
                mags_by_quality[quality].append({
                    'name': mag_name,
                    'completeness': completeness * 100 if completeness <= 1 else completeness,  # Show as percentage
                    'contamination': contamination * 100 if contamination <= 1 else contamination,  # Show as percentage
                    'taxonomy': taxonomy,
                    'size': size,  # For marker size
                    'quality': quality,
                    'is_highlighted': mag_name == highlight_mag
                })

            # Plot points for each quality level
            for quality, mags in mags_by_quality.items():
                if not mags:
                    continue
                    
                color = {
                    'High': '#2ecc71',
                    'Medium': '#f1c40f',
                    'Low': '#e74c3c',
                    'Unknown': '#95a5a6'
                }.get(quality, '#95a5a6')
                
                # Split into highlighted and non-highlighted
                highlighted = [m for m in mags if m['is_highlighted']]
                regular = [m for m in mags if not m['is_highlighted']]
                
                # Add regular MAGs
                if regular:
                    fig.add_trace(go.Scatter(
                        x=[m['completeness'] for m in regular],
                        y=[m['contamination'] for m in regular],
                        mode='markers',
                        name=quality,
                        marker=dict(
                            color=color,
                            size=10,
                            sizemode='diameter',
                            sizeref=2.0 * max([m['size'] for m in regular]) / (40**2),
                            sizemin=5
                        ),
                        text=[m['name'] for m in regular],
                        customdata=[[m['name'], m['taxonomy']] for m in regular],
                        hovertemplate="<b>%{text}</b><br>" +
                                    "Taxonomy: %{customdata[1]}<br>" +
                                    "Completeness: %{x:.1f}%<br>" +
                                    "Contamination: %{y:.1f}%<br>" +
                                    "Quality: " + quality + "<extra></extra>"
                    ))
                
                # Add highlighted MAGs with distinct appearance
                if highlighted:
                    fig.add_trace(go.Scatter(
                        x=[m['completeness'] for m in highlighted],
                        y=[m['contamination'] for m in highlighted],
                        mode='markers',
                        name=f"{quality} (Selected)",
                        marker=dict(
                            color='#e74c3c',
                            size=15,
                            line=dict(width=2, color='black'),
                            sizemode='diameter',
                            sizeref=2.0 * max([m['size'] for m in highlighted]) / (40**2),
                            sizemin=8
                        ),
                        text=[m['name'] for m in highlighted],
                        customdata=[[m['name'], m['taxonomy']] for m in highlighted],
                        hovertemplate="<b>%{text}</b><br>" +
                                    "Taxonomy: %{customdata[1]}<br>" +
                                    "Completeness: %{x:.1f}%<br>" +
                                    "Contamination: %{y:.1f}%<br>" +
                                    "Quality: " + quality + "<extra></extra>"
                    ))

            # Update layout
            fig.update_layout(
                title="MAG Quality Assessment",
                xaxis_title="Completeness (%)",
                yaxis_title="Contamination (%)",
                showlegend=True,
                height=400,
                plot_bgcolor='rgba(0,0,0,0)',
                paper_bgcolor='rgba(0,0,0,0)',
                legend=dict(
                    yanchor="bottom",
                    y=1.02,
                    xanchor="right",
                    x=1
                ),
                # Add reference lines for high quality thresholds
                shapes=[
                    # Completeness threshold (90%)
                    dict(
                        type="line",
                        x0=90, x1=90,
                        y0=0, y1=100,
                        line=dict(
                            color="rgba(46, 204, 113, 0.3)",
                            width=2,
                            dash="dash"
                        )
                    ),
                    # Contamination threshold (5%)
                    dict(
                        type="line",
                        x0=0, x1=100,
                        y0=5, y1=5,
                        line=dict(
                            color="rgba(231, 76, 60, 0.3)",
                            width=2,
                            dash="dash"
                        )
                    )
                ],
                # Add annotations explaining thresholds
                annotations=[
                    dict(
                        x=90,
                        y=2.5,
                        xref="x",
                        yref="y",
                        text="High Quality >90%",
                        showarrow=False,
                        font=dict(
                            size=10,
                            color="rgba(46, 204, 113, 1)"
                        )
                    ),
                    dict(
                        x=50,
                        y=5,
                        xref="x",
                        yref="y",
                        text="Low Contamination <5%",
                        showarrow=False,
                        font=dict(
                            size=10,
                            color="rgba(231, 76, 60, 1)"
                        )
                    )
                ]
            )

            return fig

        except Exception as e:
            logger.error(f"Error creating quality scatter plot: {str(e)}", exc_info=True)
            # Return empty plot with error message
            fig = go.Figure()
            fig.update_layout(
                title="Error in MAG Quality Plot",
                annotations=[
                    dict(
                        text=f"Error loading MAG quality data: {str(e)}",
                        showarrow=False,
                        font=dict(size=12)
                    )
                ]
            )
            return fig
    
    def _create_rna_counts_plot(self, selected_rows=None, table_data=None):
        """Create grouped bar plot of RNA gene counts."""
        try:
            if not self.mag_data:
                logger.warning("No MAG data available for creating RNA counts plot")
                return go.Figure()

            # If selected rows provided, highlight that MAG
            highlight_mag = None
            if selected_rows is not None and table_data is not None and len(selected_rows) > 0:
                highlight_mag = table_data[selected_rows[0]]['name']

            # Get RNA counts from annotation dictionary
            self.mag_data['trna_count'] = self.mag_data['annotations_dict'].apply(lambda x: len([g for g in x.get('genes', []) if g.get('type') == 'tRNA']) if isinstance(x, dict) else 0)
            self.mag_data['rrna_count'] = self.mag_data['annotations_dict'].apply(lambda x: len([g for g in x.get('genes', []) if g.get('type') == 'rRNA']) if isinstance(x, dict) else 0)
            
            df = self.mag_data.sort_values('rrna_count', ascending=False).head(20)
            
            # Create grouped bar plot
            fig = go.Figure()
            
            # Add bars for each RNA type
            for idx, row in df.iterrows():
                is_highlighted = highlight_mag is not None and row['name'] == highlight_mag
                opacity = 1.0 if is_highlighted else 0.7
                
                fig.add_trace(go.Bar(
                    name=row['name'],
                    x=['tRNA', 'rRNA'],
                    y=[row['trna_count'], row['rrna_count']],
                    showlegend=False,
                    marker_color=['#3498db', '#2ecc71'] if not is_highlighted else ['#e74c3c', '#c0392b'],
                    opacity=opacity,
                    hovertemplate='<b>' + row['name'] + '</b><br>' +
                                '%{x}: %{y}<br>' +
                                '<extra></extra>'
                ))
            
            # Update layout
            fig.update_layout(
                title='RNA Gene Content',
                xaxis_title='RNA Type',
                yaxis_title='Count',
                barmode='group',
                showlegend=False,
                hovermode='closest'
            )
            
            return fig
            
        except Exception as e:
            logger.error(f"Error creating RNA counts plot: {str(e)}", exc_info=True)
            return go.Figure()
    
    def _create_table_data(self):
        """Create data for the MAGs table with extended columns."""
        try:
            print("\n=== Starting Table Data Creation ===")
            # Get MAG data from database or other source
            print("Fetching MAG data for table creation")
            mag_data = self.get_mag_data()
            
            if not mag_data:
                logger.warning("No MAG data available for creating table data")
                return []

            # Create table data with extended columns
            table_data = []
            
            # Process only unique MAGs by name (avoid duplicates)
            processed_names = set()
            processed_count = 0
            
            print(f"Processing {len(mag_data)} MAG entries for table display")
            
            # Log the first entry's keys for debugging
            if len(mag_data) > 0:
                first_key = next(iter(mag_data))
                first_entry = mag_data[first_key]
                if isinstance(first_entry, dict):
                    print(f"Sample MAG data keys: {list(first_entry.keys())}")
                    
            for mag_id, mag_info in mag_data.items():
                try:
                    # Skip if it's a duplicate name reference (we created name-key references previously)
                    if not isinstance(mag_info, dict) or 'name' not in mag_info:
                        print(f"Skipping MAG ID {mag_id}: not a valid MAG data dictionary")
                        continue
                        
                    # Skip if we've already processed this MAG
                    mag_name = mag_info.get('name', '')
                    if mag_name in processed_names:
                        print(f"Skipping duplicate MAG name: {mag_name}")
                        continue
                    processed_names.add(mag_name)
                    
                    # Extract basic info
                    bin_id = mag_info.get('bin_id', mag_id)
                    
                    # Get numeric values with safe defaults
                    try:
                        completeness = float(mag_info.get('completeness', 0))
                    except (ValueError, TypeError):
                        completeness = 0.0
                        print(f"Warning: Invalid completeness value for MAG {mag_name}")
                        
                    try:
                        contamination = float(mag_info.get('contamination', 0))
                    except (ValueError, TypeError):
                        contamination = 0.0
                        print(f"Warning: Invalid contamination value for MAG {mag_name}")
                        
                    try:
                        size_mb = float(mag_info.get('size', 0)) / 1000000
                    except (ValueError, TypeError):
                        size_mb = 0.0
                        print(f"Warning: Invalid size value for MAG {mag_name}")
                        
                    try:
                        gc_content = float(mag_info.get('gc_content', 0))
                    except (ValueError, TypeError):
                        gc_content = 0.0
                        print(f"Warning: Invalid GC content value for MAG {mag_name}")
                        
                    try:
                        n_contigs = int(mag_info.get('n_contigs', 0))
                    except (ValueError, TypeError):
                        n_contigs = 0
                        print(f"Warning: Invalid contig count for MAG {mag_name}")
                        
                    try:
                        n50 = int(mag_info.get('n50', 0))
                    except (ValueError, TypeError):
                        n50 = 0
                        print(f"Warning: Invalid N50 value for MAG {mag_name}")
                        
                    try:
                        longest_contig = int(mag_info.get('longest_contig', 0))
                    except (ValueError, TypeError):
                        longest_contig = 0
                        print(f"Warning: Invalid longest contig value for MAG {mag_name}")
                    
                    # Get taxonomy directly from the mag_info
                    taxonomy = mag_info.get('taxonomy', "Unknown")
                    if not taxonomy or taxonomy.lower() in ['unknown', 'unclassified', 'none', '']:
                        taxonomy = "Unknown"
                    
                    # Make sure taxonomy is displayed as a string
                    if isinstance(taxonomy, dict):
                        # Format taxonomy dictionary nicely
                        tax_parts = []
                        for level, name in taxonomy.items():
                            if name and name.lower() not in ['unknown', 'unclassified', 'none']:
                                tax_parts.append(f"{name}")
                        taxonomy = "; ".join(tax_parts) if tax_parts else "Unknown"
                    
                    # Get quality assessment
                    quality = mag_info.get('quality', 'unknown')
                    if quality == 'unknown':
                        quality = self._assess_quality({
                            'completeness': completeness,
                            'contamination': contamination
                        })
                    
                    # Add all available data
                    record = {
                        "id": mag_id,
                        "name": mag_name,
                        "bin_id": bin_id,
                        "taxonomy": taxonomy,
                        "completeness": completeness / 100 if completeness > 1 else completeness,  # Ensure percentage format (0-1)
                        "contamination": contamination / 100 if contamination > 1 else contamination,  # Ensure percentage format (0-1)
                        "size_mb": size_mb,
                        "gc_content": gc_content / 100 if gc_content > 1 else gc_content,  # Ensure percentage format (0-1)
                        "num_contigs": n_contigs,
                        "n50": n50,
                        "longest_contig": longest_contig,
                        "quality": quality,
                        # Include raw data for callbacks
                        "fasta_file": mag_info.get('fasta_file', '')
                    }
                    table_data.append(record)
                    processed_count += 1
                    
                    # Print progress every 10 MAGs
                    if processed_count % 10 == 0:
                        print(f"Processed {processed_count} unique MAGs...")
                    
                except Exception as e:
                    print(f"Error processing MAG {mag_id} for table display: {str(e)}")
                    continue
            
            # Sort by completeness (highest first)
            if table_data:
                table_data.sort(key=lambda x: x.get('completeness', 0), reverse=True)
                print(f"\nSuccessfully created table data for {processed_count} MAGs")
                
                # Print quality distribution
                quality_counts = {}
                for mag in table_data:
                    quality = mag.get('quality', 'unknown')
                    quality_counts[quality] = quality_counts.get(quality, 0) + 1
                    
                print("\nQuality distribution in final table data:")
                for quality, count in quality_counts.items():
                    print(f"- {quality.capitalize()}: {count} MAGs")
            else:
                print("No MAG data was successfully processed for table display")
            
            print("\n=== Table Data Creation Complete ===")
            return table_data

        except Exception as e:
            print(f"\n[ERROR] Error creating table data: {str(e)}")
            print(traceback.format_exc())
            # Return empty list with proper structure
            return []

    def parse_fasta_for_viewer(self, fasta_data):
        """Parse FASTA data into a format suitable for the sequence viewer."""
        print("\n===== Parsing FASTA Data =====")
        
        if not fasta_data:
            print("Empty FASTA data provided to parser")
            return []
        
        try:
            sequences = []
            current_header = None
            current_sequence = ""
            
            # Log first 100 chars of FASTA content for debugging
            fasta_sample = fasta_data[:100] + "..." if len(fasta_data) > 100 else fasta_data
            print(f"FASTA content sample: {fasta_sample}")
            print(f"FASTA content length: {len(fasta_data)}")
            
            # Check if the FASTA file is properly formatted
            if '>' not in fasta_data:
                print("WARNING: FASTA data doesn't contain '>' character - may not be properly formatted")
                # If this appears to be raw sequence data without headers, create a single sequence
                if len(fasta_data) > 1000:  # Only for reasonably sized sequences
                    print("Treating data as a raw sequence without FASTA headers")
                    seq_id = f"seq_1"
                    sequences.append({
                        'header': f'>{seq_id}',
                        'sequence': fasta_data.strip(),
                        'length': len(fasta_data.strip()),
                        'id': seq_id
                    })
                    print(f"Created sequence entry from raw data, length: {len(fasta_data.strip())}")
                    return sequences
            
            line_count = 0
            sequence_count = 0
            
            # Normal FASTA parsing
            for line in fasta_data.split('\n'):
                line_count += 1
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('>'):
                    # Process previous sequence if it exists
                    if current_header and current_sequence:
                        sequence_count += 1
                        seq_id = current_header.split()[0].replace('>', '')
                        sequences.append({
                            'header': current_header,
                            'sequence': current_sequence,
                            'length': len(current_sequence),
                            'id': seq_id
                        })
                        
                        # Log every 10th sequence for debugging
                        if sequence_count % 10 == 0:
                            print(f"Processed {sequence_count} sequences, latest: {seq_id} ({len(current_sequence)} bp)")
                    
                    current_header = line
                    current_sequence = ""
                else:
                    current_sequence += line
            
            # Add the last sequence
            if current_header and current_sequence:
                seq_id = current_header.split()[0].replace('>', '')
                sequences.append({
                    'header': current_header,
                    'sequence': current_sequence,
                    'length': len(current_sequence),
                    'id': seq_id
                })
                sequence_count += 1
            
            # If no sequences were found but we have substantial content
            if not sequences and len(fasta_data) > 1000:
                print("No sequences parsed but data is substantial - attempting alternate parsing")
                # Check if this might be a single sequence without proper formatting
                if '>' in fasta_data:
                    # Try to extract a sequence with a cleaner approach
                    parts = fasta_data.split('>', 1)
                    if len(parts) > 1:
                        header_line = parts[1].split('\n', 1)
                        if len(header_line) > 1:
                            header = header_line[0]
                            sequence = header_line[1].replace('\n', '')
                            
                            seq_id = header.split()[0]
                            sequences.append({
                                'header': f'>{header}',
                                'sequence': sequence,
                                'length': len(sequence),
                                'id': seq_id
                            })
                            print(f"Successfully extracted sequence with alternate method, length: {len(sequence)}")
                else:
                    # Treat whole content as one sequence
                    seq_id = "sequence_1"
                    sequences.append({
                        'header': f'>{seq_id}',
                        'sequence': fasta_data.strip(),
                        'length': len(fasta_data.strip()),
                        'id': seq_id
                    })
                    print(f"Created single sequence from content, length: {len(fasta_data.strip())}")
                
            # Sort sequences by length (largest first)
            sequences.sort(key=lambda x: x['length'], reverse=True)
            
            # Log sequence stats
            if sequences:
                total_length = sum(s['length'] for s in sequences)
                max_length = max(s['length'] for s in sequences)
                min_length = min(s['length'] for s in sequences)
                print(f"Parsed {len(sequences)} sequences from {line_count} lines")
                print(f"Total sequence length: {total_length:,} bp")
                print(f"Sequence length range: {min_length:,} - {max_length:,} bp")
                
                # Log IDs of largest sequences
                if len(sequences) > 0:
                    largest_seq_ids = [s['id'] for s in sequences[:3]]
                    print(f"Largest sequences: {largest_seq_ids}")
            else:
                print("No sequences found in FASTA data")
                
            print(f"===== FASTA Parsing Complete =====\n")
            return sequences
        except Exception as e:
            print(f"Error parsing FASTA data: {str(e)}")
            print(traceback.format_exc())
            
            # Create an emergency sequence if we have substantial data but parsing failed
            if fasta_data and len(fasta_data) > 1000:
                print("Exception during parsing - creating emergency sequence from available data")
                seq_id = "emergency_sequence"
                # Clean up the data as best we can
                sequence_data = fasta_data
                if '>' in sequence_data:
                    # Remove everything before the first '>' and try to extract a sequence
                    parts = sequence_data.split('>', 1)
                    if len(parts) > 1:
                        try:
                            header_line = parts[1].split('\n', 1)
                            if len(header_line) > 1:
                                sequence_data = header_line[1].replace('\n', '')
                                seq_id = header_line[0].split()[0]
                        except:
                            # If that fails, just use the raw data without headers
                            sequence_data = sequence_data.replace('>', '').replace('\n', '')
                return [{
                    'header': f'>{seq_id}',
                    'sequence': sequence_data.strip(),
                    'length': len(sequence_data.strip()),
                    'id': seq_id
                }]
            
            return []

    def _calculate_n50(self, sequences):
        """Calculate N50 for a list of sequences."""
        if not sequences:
            return 0
            
        # Get sequence lengths sorted in descending order
        seq_lengths = sorted([seq.get('length', 0) for seq in sequences], reverse=True)
        total_length = sum(seq_lengths)
        
        # Calculate N50
        running_sum = 0
        for length in seq_lengths:
            running_sum += length
            if running_sum >= total_length / 2:
                return length
                
        return 0
        
    def _get_longest_contig(self, sequences):
        """Get the length of the longest contig."""
        if not sequences:
            return 0
            
        return max([seq.get('length', 0) for seq in sequences])

    def create_layout(self):
        """Create the layout for the MAGs app.
        
        Returns:
            A Dash layout component
        """
        try:
            print("Creating MAGs layout")
            
            # Create the mag table data
            table_data = self._create_table_data()
            
            # Create initial figures
            initial_figures = self._create_initial_figures()
            
            # Create quality metrics scatter plot
            quality_scatter = self._create_quality_scatter()
            
            # Create content section
            mag_content = html.Div([
                html.H5("MAG Visualization", className="mb-3"),
                
                html.Hr(),
                html.Div([
                    html.H5("Genome Visualization", id="circos-plot-title"),
                    html.Iframe(
                        id="circos-plot-iframe",
                        srcDoc="""
                        <html>
                        <body style="text-align: center; font-family: Arial, sans-serif; background-color: #f9f9f9;">
                            <div style="padding: 40px 20px;">
                                <div style="font-size: 18px; margin-bottom: 15px;">Select a MAG from the table to view genome visualization</div>
                                <div style="color: #666; font-size: 14px;">The circular plot will show contigs and gene annotations when a MAG is selected</div>
                            </div>
                        </body>
                        </html>
                        """,
                        style={
                            'width': '100%',
                            'height': '800px',
                            'border': 'none',
                            'padding': '0',
                            'margin': '0',
                            'position': 'relative',
                            'z-index': 1
                        },
                        sandbox='allow-scripts allow-same-origin allow-popups'
                    )
                ]),
                
                # MAG info section (initially hidden)
                html.Div([
                    html.H5("MAG Sequence Information", id="mag-info-title"),
                    html.Div(id="mag-genome-warning", style={"display": "none"}),
                    
                    # Main info container
                    html.Div([
                        # Dropdown sequence selector
                        html.Div([
                            html.Label("Select Sequence"),
                            dcc.Dropdown(
                                id="sequence-selector",
                                options=[],
                                placeholder="Select a sequence",
                            ),
                        ], className="mb-3"),
                        
                        # Sequence viewer
                        html.Div(id="sequence-viewer-container", className="mb-4"),
                        
                        # Gene search
                        html.Div([
                            html.Label("Search Genes"),
                            html.Div([
                                dcc.Input(
                                    id="gene-search",
                                    type="text",
                                    placeholder="Search by gene name, product, or function",
                                    className="form-control",
                                    style={"width": "100%", "marginRight": "10px"}
                                ),
                                html.Button(
                                    "Search", 
                                    id="gene-search-button",
                                    className="btn btn-primary",
                                    style={"marginTop": "10px"}
                                )
                            ], style={"display": "flex", "flexDirection": "column"}),
                            html.Div(id="gene-search-results", className="mt-3"),
                        ], className="mb-4"),
                    ], id="mag-info-section", style={"display": "none"})
                ], style={
                    'marginBottom': '30px', 
                    'padding': '10px', 
                    'backgroundColor': '#f9f9f9',
                    'borderRadius': '5px'
                }),
                
                
                # Debug information (hidden initially)
                html.Div(id="mags-debug-info", style={"fontSize": "12px", "margin": "20px 0", "display": "none"}),
                
                # Hidden components
                html.Div(id="mag-content", style={"display": "none"}),
                html.Div(id="sequence-data", style={"display": "none"}),
                html.Div(id="fasta-data", style={"display": "none"}),
                html.Div(id="copy-alert", style={"display": "none"}),
                
            ], id="mag-content-container")
            
            # Create the full layout including the sidebar
            layout = dmc.Container([
                dmc.Grid([
                    # Left column with controls and table
                    dmc.GridCol([
                        html.H4("MAGs Viewer"),
                        html.P("This tool allows you to explore Metagenome-Assembled Genomes (MAGs) in your dataset.",
                              style={"color": "#666"}),
                        
                        # Quality filter control
                        dmc.Paper([
                            dmc.Stack([
                                dmc.Text("Quality Filter", fw=500),
                                dmc.RadioGroup(
                                    id="quality-filter",
                                    value="all",
                                    children=[
                                        dmc.Radio(label="All", value="all"),
                                        dmc.Radio(label="High", value="high"),
                                        dmc.Radio(label="Medium", value="medium"),
                                        dmc.Radio(label="Low", value="low"),
                                    ],
                                ),
                            ]),
                        ], p="md", withBorder=True, style={"marginBottom": "20px"}),
                        
                        # MAG table - using Dash DataTable for better interaction
                        dmc.Paper([
                            html.Div([
                                dash_table.DataTable(
                                    id="mag-table",
                                    columns=[
                                        {"name": "Name", "id": "name"},
                                        {"name": "Completeness", "id": "completeness", "type": "numeric", "format": Format(precision=1, scheme=Scheme.percentage)},
                                        {"name": "Contamination", "id": "contamination", "type": "numeric", "format": Format(precision=1, scheme=Scheme.percentage)},
                                        {"name": "Quality", "id": "quality"},
                                        {"name": "Contigs", "id": "num_contigs", "type": "numeric", "format": Format(group=True)},
                                        {"name": "Size (Mb)", "id": "size_mb", "type": "numeric", "format": Format(precision=2)},
                                        {"name": "GC Content", "id": "gc_content", "type": "numeric", "format": Format(precision=1, scheme=Scheme.percentage)},
                                        {"name": "Taxonomy", "id": "taxonomy"},
                                    ],
                                    data=table_data,
                                    style_table={"overflowX": "auto"},
                                    style_cell={
                                        "textAlign": "left",
                                        "overflow": "hidden",
                                        "textOverflow": "ellipsis",
                                    },
                                    style_cell_conditional=[
                                        {"if": {"column_id": "completeness"}, "textAlign": "right"},
                                        {"if": {"column_id": "contamination"}, "textAlign": "right"},
                                        {"if": {"column_id": "num_contigs"}, "textAlign": "right"},
                                        {"if": {"column_id": "size_mb"}, "textAlign": "right"},
                                        {"if": {"column_id": "gc_content"}, "textAlign": "right"},
                                    ],
                                    style_data_conditional=[
                                        {
                                            "if": {"filter_query": "{quality} = 'high'"},
                                            "backgroundColor": "rgba(46, 204, 113, 0.1)"
                                        },
                                        {
                                            "if": {"filter_query": "{quality} = 'medium'"},
                                            "backgroundColor": "rgba(241, 196, 15, 0.1)"
                                        },
                                        {
                                            "if": {"filter_query": "{quality} = 'low'"},
                                            "backgroundColor": "rgba(231, 76, 60, 0.1)"
                                        },
                                    ],
                                    style_header={
                                        "backgroundColor": "rgb(230, 230, 230)",
                                        "fontWeight": "bold"
                                    },
                                    row_selectable="single",
                                    selected_rows=[],  # No rows selected by default
                                    page_action="native",
                                    page_size=20,
                                    sort_action="native",
                                    sort_mode="multi",
                                ),
                                # Store the full table data for reference
                                html.Div(id="mag-table-data", children=json.dumps(table_data), style={"display": "none"}),
                            ])
                        ], p="md", withBorder=True, style={"marginBottom": "20px"}),
                        
                        # Quality scatter plot
                        dmc.Paper([
                            dmc.Stack([
                                dmc.Text("MAG Quality Assessment", fw=500),
                                dcc.Graph(
                                    id="quality-scatter",
                                    figure=quality_scatter,
                                    style={"height": "400px"},
                                ),
                            ]),
                        ], p="md", withBorder=True),
                        
                    ], span=5, p="md"),
                    
                    # Right column with MAG info and visualizations
                    dmc.GridCol([
                        # MAG Details Section
                        html.Div([
                            # Circos Plot - this is an OUTPUT of the callback
                            html.Div(id="circos-plot", children=[
                                html.Div(
                                    "Select a MAG from the table to view genome visualization",
                                    style={
                                        "textAlign": "center",
                                        "padding": "40px 20px",
                                        "fontSize": "18px",
                                        "color": "#666"
                                    }
                                )
                            ]),
                            
                            # Sequence Information - this is an OUTPUT of the callback
                            html.Div(id="mag-sequence-info"),
                            
                            # Gene Content - this is an OUTPUT of the callback
                            html.Div(id="mag-gene-content")
                        ], style={"marginTop": "20px"}),
                        
                        # MAG content panel (holds circos plot and sequence info)
                        dmc.Paper([
                            dmc.Stack([
                                dmc.Text("MAG Content", fw=500),
                                dmc.Divider(),
                                mag_content,
                            ]),
                        ], p="md", withBorder=True),
                        
                    ], span=7, p="md"),
                ]),
            ], fluid=True, p=0, style={"marginTop": "20px"})
            
            return layout
            
        except Exception as e:
            print(f"Error in MAGs layout creation: {str(e)}")
            print(traceback.format_exc())
            return html.Div(
                dmc.Alert(
                    children=[
                        html.H5("Error Creating MAGs Viewer"),
                        html.P(f"An error occurred: {str(e)}")
                    ],
                    title="Error",
                    color="red"
                )
            )

    def _create_test_circos_plot(self, mag_id=None):
        """Create an NGCircos plot showing MAG/contig data."""
        try:
            # Clear logging to help debug
            logger.info(f"Creating circos plot for MAG ID: {mag_id if mag_id else 'None'}")
            
            # If no MAG ID is provided, don't try to create a plot with data
            if not mag_id:
                logger.info("No MAG ID provided, returning placeholder message")
                return """
                <html>
                <body style="text-align: center; font-family: Arial, sans-serif; background-color: #f9f9f9;">
                    <div style="padding: 40px 20px;">
                        <div style="font-size: 18px; margin-bottom: 15px;">Select a MAG from the table to view genome visualization</div>
                        <div style="color: #666; font-size: 14px;">The circular plot will show contigs and gene annotations when a MAG is selected</div>
                    </div>
                </body>
                </html>
                """
            
            # Get MAG data including sequence information
            mag_info = None
            sequences = []
            
            # Default empty annotations
            gene_annotations = []
            
            # Force refresh the mag_data to avoid using cached data
            # This is critical to ensure we're not using stale data for different MAGs
            mag_data = self.get_mag_data(include_fasta=True, force_refresh=True)
            
            try:
                # Get the MAG data for the specific MAG ID
                logger.info(f"Retrieving MAG data for {mag_id}")
                
                if not mag_data:
                    logger.warning("No MAG data returned from get_mag_data")
                    return """
                    <html>
                    <body style="text-align: center; font-family: Arial, sans-serif; background-color: #f9f9f9;">
                        <div style="padding: 40px 20px;">
                            <div style="font-size: 18px; color: #e74c3c; margin-bottom: 15px;">No MAG data available</div>
                            <div style="color: #666; font-size: 14px;">Please check database connection or data availability</div>
                        </div>
                    </body>
                    </html>
                    """
                    
                if mag_id in mag_data:
                    mag_info = mag_data[mag_id]
                    logger.info(f"Found MAG info for {mag_id}")
                else:
                    # Try case-insensitive match if not found
                    mag_id_lower = str(mag_id).lower()
                    for key, info in mag_data.items():
                        if isinstance(info, dict) and str(info.get('bin_id', '')).lower() == mag_id_lower:
                            mag_info = info
                            logger.info(f"Found MAG info for {mag_id} via case-insensitive match")
                            break
                    
                    if not mag_info:
                        logger.warning(f"MAG ID {mag_id} not found in MAG data dictionary")
                        return f"""
                        <html>
                        <body style="text-align: center; font-family: Arial, sans-serif; background-color: #f9f9f9;">
                            <div style="padding: 40px 20px;">
                                <div style="font-size: 18px; color: #e74c3c; margin-bottom: 15px;">MAG ID {mag_id} not found</div>
                                <div style="color: #666; font-size: 14px;">Please select a different MAG</div>
                            </div>
                        </body>
                        </html>
                        """
                    
                # Get sequences if we found the MAG
                if mag_info:
                    # Get sequences
                    if 'fasta' in mag_info and isinstance(mag_info['fasta'], list):
                        sequences = mag_info['fasta']
                        logger.info(f"Found {len(sequences)} parsed sequences for MAG {mag_id}")
                    elif 'fasta_file' in mag_info and mag_info['fasta_file']:
                        sequences = self.parse_fasta_for_viewer(mag_info['fasta_file'])
                        logger.info(f"Parsed {len(sequences)} sequences from FASTA file for MAG {mag_id}")
                    else:
                        logger.warning(f"No FASTA data found for MAG {mag_id}")
                        return f"""
                        <html>
                        <body style="text-align: center; font-family: Arial, sans-serif; background-color: #f9f9f9;">
                            <div style="padding: 40px 20px;">
                                <div style="font-size: 18px; color: #e74c3c; margin-bottom: 15px;">No sequence data found for MAG {mag_id}</div>
                                <div style="color: #666; font-size: 14px;">This MAG doesn't have any FASTA sequence data</div>
                            </div>
                        </body>
                        </html>
                        """
            except Exception as e:
                logger.error(f"Error retrieving MAG data: {str(e)}", exc_info=True)
                return f"""
                <html>
                <body style="text-align: center; font-family: Arial, sans-serif; background-color: #f9f9f9;">
                    <div style="padding: 40px 20px;">
                        <div style="font-size: 18px; color: #e74c3c; margin-bottom: 15px;">Error retrieving MAG data</div>
                        <div style="color: #666; font-size: 14px;">Error details: {str(e)}</div>
                    </div>
                </body>
                </html>
                """
            
            # Try to get GFF data for annotations
            try:
                logger.info(f"Retrieving GFF data for {mag_id}")
                
                # First check if GFF data is directly available in the mag_info
                if mag_info and 'gff_file' in mag_info and mag_info['gff_file']:
                    logger.info(f"Found GFF data in mag_info for {mag_id}, parsing...")
                    gff_content = mag_info['gff_file']
                    gene_annotations = self._parse_gff_for_annotations(gff_content, sequences)
                    logger.info(f"Parsed {len(gene_annotations)} annotations from GFF data")
                else:
                    # Try the database method
                    logger.info(f"Using get_gff_data method for MAG {mag_id}")
                    gff_data = self.get_gff_data(mag_id)
                    if gff_data:
                        logger.info(f"Found GFF data for MAG {mag_id}")
                        if isinstance(gff_data, str):
                            # Parse the GFF content if it's a string
                            gene_annotations = self._parse_gff_for_annotations(gff_data, sequences)
                            logger.info(f"Parsed {len(gene_annotations)} annotations from GFF string")
                        elif isinstance(gff_data, list):
                            # GFF data is already a list of feature dictionaries
                            gene_annotations = gff_data
                            logger.info(f"Using {len(gene_annotations)} annotations from GFF list data")
                        else:
                            logger.warning(f"Unexpected GFF format for MAG {mag_id}: {type(gff_data)}")
                    else:
                        logger.warning(f"No GFF data found for MAG {mag_id}")
                
                # If no annotations were found, try direct GFF parsing from database
                if not gene_annotations:
                    try:
                        from users.models import MAG
                        logger.info(f"Trying to get GFF directly from database for MAG {mag_id}")
                        mag_obj = MAG.objects.get(name=mag_id, user_id=self.user_id)
                        
                        if mag_obj.gff_file:
                            logger.info(f"Found GFF in database for MAG {mag_id}, parsing...")
                            gene_annotations = self._parse_gff_for_annotations(mag_obj.gff_file, sequences)
                            logger.info(f"Parsed {len(gene_annotations)} annotations directly from database GFF")
                        else:
                            logger.warning(f"No GFF file found in database for MAG {mag_id}")
                    except Exception as db_error:
                        logger.error(f"Error retrieving GFF from database: {str(db_error)}")
                
                # If we have annotations, sort them by position to ensure better distribution
                if gene_annotations:
                    # First group by contig
                    annotations_by_contig = {}
                    for ann in gene_annotations:
                        contig = ann.get('contig', 'unknown')
                        if contig not in annotations_by_contig:
                            annotations_by_contig[contig] = []
                        annotations_by_contig[contig].append(ann)
                    
                    # Then sort each contig's annotations by position
                    for contig in annotations_by_contig:
                        annotations_by_contig[contig].sort(key=lambda x: x.get('start', 0))
                    
                    # Flatten back to a single list, ensuring we get genes from all positions
                    gene_annotations = []
                    # Determine how many to take from each contig in each pass
                    max_passes = 10  # Number of passes through all contigs
                    for _ in range(max_passes):
                        for contig in annotations_by_contig:
                            if annotations_by_contig[contig]:
                                chunk_size = max(1, len(annotations_by_contig[contig]) // max_passes)
                                gene_annotations.extend(annotations_by_contig[contig][:chunk_size])
                                annotations_by_contig[contig] = annotations_by_contig[contig][chunk_size:]
                    
                    # Add any remaining annotations
                    for contig in annotations_by_contig:
                        gene_annotations.extend(annotations_by_contig[contig])
                
                # As a last resort, create placeholder annotations
                if not gene_annotations and sequences:
                    logger.info(f"No annotations found for MAG {mag_id}, creating placeholder annotations")
                    gene_annotations = self._create_placeholder_annotations(sequences)
                    logger.info(f"Created {len(gene_annotations)} placeholder annotations")
                    
                # Log annotation details
                if gene_annotations:
                    feature_types = {}
                    for ann in gene_annotations[:50]:  # Only check a sample for large datasets
                        ann_type = ann.get('type', 'unknown')
                        feature_types[ann_type] = feature_types.get(ann_type, 0) + 1
                    logger.info(f"Annotation types: {feature_types}")
                    
                    # Log a sample of annotations
                    logger.info("Sample annotations:")
                    for ann in gene_annotations[:3]:
                        logger.info(f"  {ann.get('type', 'unknown')} - {ann.get('name', 'no_name')} ({ann.get('product', 'no_product')[:50]}...)")
                else:
                    logger.warning(f"Still no annotations available for MAG {mag_id}")
                    
            except Exception as e:
                logger.error(f"Error processing GFF data: {str(e)}", exc_info=True)
                # Continue with empty annotations - we'll show the contigs at least
                gene_annotations = []
            
            # Sort sequences by length and get top N contigs
            genome_data = []
            if sequences:
                try:
                    # Sort by sequence length (if available)
                    sequences.sort(key=lambda x: x.get('length', 0), reverse=True)
                    # Limit to top 24 sequences for better visualization
                    sequences = sequences[:24]
                    logger.info(f"Sorted and limited to {len(sequences)} sequences")
                    
                    # Format the contig data for NGCircos
                    for seq in sequences:
                        seq_id = seq.get('id', '').split()[0]
                        length = seq.get('length', 0)
                        if seq_id and length > 0:
                            genome_data.append([self._sanitize_text(seq_id), length])
                    
                    # Log the lengths to help with debugging
                    for contig_name, contig_length in genome_data[:5]:
                        logger.info(f"Contig {contig_name} length: {contig_length}")
                    
                    # Check if any features are beyond 500kb
                    if gene_annotations:
                        max_pos = max([ann.get('end', 0) for ann in gene_annotations])
                        logger.info(f"Maximum feature position found: {max_pos}")
                    
                    logger.info(f"Created genome_data with {len(genome_data)} contigs")
                except Exception as e:
                    logger.error(f"Error processing sequences: {str(e)}", exc_info=True)
            
            # If no valid genome data, use sample data
            if not genome_data:
                logger.info("Using sample genome data as fallback")
                genome_data = [
                    ["contig_1", 5000000],  # Use larger genome sizes
                    ["contig_2", 4000000],
                    ["contig_3", 3000000],
                    ["contig_4", 2000000],
                    ["contig_5", 1500000],
                    ["contig_6", 1000000]
                ]
            
            # Generate a unique ID based on the MAG ID to avoid caching issues
            # Sanitize the MAG ID to create a valid CSS selector by replacing periods with underscores
            safe_mag_id = str(mag_id).replace('.', '_').replace(':', '_').replace(' ', '_') if mag_id else 'sample'
            unique_id = f"ngcircos-{safe_mag_id}-{int(time.time())}"
            
            # Test plot initialization JavaScript
            html_content = ("""
                <!DOCTYPE html>
                <html>
                <head>
                    <meta charset="utf-8">
                    <meta name="viewport" content="width=device-width, initial-scale=1">
                    <title>MAG Visualization</title>
                    <!-- Load scripts with defer attribute to ensure they load efficiently -->
                    <script defer src="https://code.jquery.com/jquery-2.1.3.js"></script>
                    <script defer src="https://d3js.org/d3.v3.js"></script>
                    <script defer src="https://cdn.jsdelivr.net/gh/YaCui/NG-Circos@master/NGCircos.js"></script>
                </head>
                <body>
                    <div id="NGCircos-""" + unique_id + """" style="width:100%; height:600px;"></div>
                    <div id="debug-info-""" + unique_id + """" style="background:#f8f9fa; padding:10px; margin-top:10px; border:1px solid #ddd; font-family:monospace; font-size:12px;">
                        <strong>MAG ID:</strong> """ + str(safe_mag_id) + """<br>
                        <strong>Contigs:</strong> """ + str(len(genome_data)) + """<br>
                        <strong>Annotations:</strong> """ + str(len(gene_annotations)) + """<br>
                        <strong>Time:</strong> """ + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + """<br>
                    </div>
                    <script>
                        // Wait for window load to ensure all resources are available
                        window.onload = function() {
                            console.log('Window loaded, initializing circos plot for """ + str(safe_mag_id) + """');
                            
                            // Ensure libraries are loaded before proceeding
                            function checkLibraries() {
                                if (window.jQuery && window.d3 && window.NGCircos) {
                                    console.log('All libraries loaded, creating plot...');
                                    createTestPlot();
                                } else {
                                    console.log('Waiting for libraries to load...');
                                    setTimeout(checkLibraries, 100);
                                }
                            }
                            
                            // Check for required libraries
                            checkLibraries();

                            // Debug function to help troubleshoot
                            function showDebugInfo(message) {
                                console.log(message);
                                var debugDiv = document.getElementById('debug-info-""" + unique_id + """');
                                if (debugDiv) {
                                    debugDiv.innerHTML += message + '<br>';
                                }
                            }

                            function createTestPlot() {
                                try {
                                    showDebugInfo('Creating MAG visualization for """ + str(safe_mag_id) + """...');
                                    
                                    // Load genome data
                                    var genomeData = """ + json.dumps(genome_data) + """;
                                    
                                    // Format for NGCircos
                                    var NGCircosGenome = [genomeData];
                                    showDebugInfo('Using ' + genomeData.length + ' contigs');
                                    
                                    // Add gene annotations if available
                                    if (""" + str(len(gene_annotations) > 0).lower() + """) {
                                        showDebugInfo('Adding ' + """ + str(len(gene_annotations)) + """ + ' gene annotations');
                                        
                                        // Create GENE01 data array in NGCircos format
                                        var GENE01 = [ "GENE01" , {
                                            outerRadius: 300,  // Place genes inside the genome track
                                            innerRadius: 200,  // Make the track wider
                                            arrow: true,
                                            arrowGap: 6,
                                            arrowColor: "black",
                                            arrowSize: "10px",
                                            cdsColor: "#3498db",
                                            cdsStrokeColor: "#2980b9",
                                            cdsStrokeWidth: 0.5,
                                            utrWidth: 6,
                                            utrColor: "#e74c3c",
                                            utrStrokeColor: "#c0392b",
                                            utrStrokeWidth: 0.5,
                                            GENEAnimationDisplay: true,
                                            GENEAnimationType: "linear",
                                            GENEAnimationTime: 2000,
                                            GENEAnimationDelay: 0.5,
                                        } , [
                                            """ + ",\n                                            ".join(["""
                                            {chr: "%s", strand: "%s", start: "%s", end: "%s", type: "%s", name: "%s", product: "%s", link: ""}
                                            """ % (
                                                self._sanitize_text(ann.get('contig', 'unknown')), 
                                                '+' if ann.get('strand', '+') == '+' or ann.get('strand', 1) == 1 else '-',
                                                ann.get('start', 0), 
                                                ann.get('end', 0),
                                                self._sanitize_text(ann.get('type', 'gene')),
                                                self._sanitize_text(ann.get('name', '') or ann.get('gene', '') or ann.get('locus_tag', '') or f"{ann.get('type', 'gene')}_{ann.get('start', 0)}"),
                                                self._sanitize_text(ann.get('product', '') or ann.get('function', '') or '')
                                            ) for ann in gene_annotations[:5000]]) + """
                                        ]];
                                        
                                        // Add a GC content track
                                        var GC_CONTENT = [ "GC_CONTENT" , {
                                            maxRadius: 170,
                                            minRadius: 140,
                                            LineColor: "#2ecc71", 
                                            LineWidth: 2,
                                            maxGap: 1000,
                                            direction: "out"
                                        } , [""" + 
                                        ",".join([
                                            """
                                            {chr: "%s", pos: "%s", value: "%.2f"}
                                            """ % (
                                                self._sanitize_text(seq.get('id', '').split()[0]),
                                                int(seq.get('length', 1000) / 2),
                                                seq.get('gc_content', 50.0) / 100.0
                                            ) for seq in sequences[:24] if seq.get('id') and seq.get('gc_content')
                                        ]) + """
                                        ]];
                                        
                                        // Add a GC Skew track
                                        var GC_SKEW = [ "GC_SKEW" , {
                                            maxRadius: 130,
                                            minRadius: 100,
                                            LineColor: "#e74c3c",
                                            LineWidth: 1.5,
                                            axisBorderColor: "#95a5a6",
                                            axisBorderWidth: 0.5
                                        }, [""" + 
                                        ",".join([
                                            """
                                            {chr: "%s", pos: "%s", value: "%s"}
                                            """ % (
                                                self._sanitize_text(seq.get('id', '').split()[0]),
                                                int(seq.get('length', 1000) / 2),
                                                random.uniform(-0.2, 0.2)
                                            ) for seq in sequences[:24] if seq.get('id')
                                        ]) + """
                                        ]];
                                    }
                                    
                                    // Initialize NGCircos with configuration
                                    var NGCircos1 = new NGCircos(""" + ('GENE01, GC_CONTENT, GC_SKEW, ' if len(gene_annotations) > 0 else '') + """NGCircosGenome, {
                                        target: "NGCircos-""" + unique_id + """",
                                        svgWidth: 1000,
                                        svgHeight: 1000,
                                        svgClassName: "gene",
                                        chrPad: 0.04,
                                        innerRadius: 320,
                                        outerRadius: 370,
                                        zoom: true,
                                        genomeFillColor: ["#34495e", "#2c3e50", "#1abc9c", "#16a085", "#3498db", "#2980b9", "#9b59b6", "#8e44ad", "#f1c40f"],
                                        ticks: {
                                            display: true,
                                            len: 5,
                                            color: "#2c3e50",
                                            textSize: 10,
                                            textColor: "#2c3e50",
                                            scale: 3110456,
                                            interval: 6,
                                            realLength: true,
                                        },
                                        genomeLabel: {
                                            display: true,
                                            textSize: 14,
                                            textColor: "#2c3e50",
                                            dx: 0.028,
                                            dy: "0.55em",
                                        },
                                        genomeBorder: {
                                            display: true,
                                            borderColor: "#2c3e50",
                                            borderSize: 0.8,
                                        },
                                        GENEMouseEvent: true,
                                        GENEMouseClickDisplay: true,
                                        GENEMouseClickColor: "none",
                                        GENEMouseClickArcOpacity: 1.0,
                                        GENEMouseClickArcStrokeColor: "#e67e22",
                                        GENEMouseClickArcStrokeWidth: 2,
                                        GENEMouseClickTextFromData: "sixth",

                                        GENEMouseClickTextOpacity: 1,
                                        GENEMouseClickTextColor: "#2c3e50",
                                        GENEMouseClickTextSize: 12,
                                        GENEMouseClickTextPostionX: 0,
                                        GENEMouseClickTextPostionY: 0,
                                        GENEMouseClickTextDrag: true,

                                        GENEMouseOutDisplay: true,
                                        GENEMouseOutAnimationTime: 500,
                                        GENEMouseOutColor: "none",
                                        GENEMouseOutArcOpacity: 1.0,
                                        GENEMouseOutArcStrokeColor: "none",
                                        GENEMouseOutArcStrokeWidth: 0,

                                        GENEMouseOverDisplay: true,
                                        GENEMouseOverColor: "none",
                                        GENEMouseOverArcOpacity: 1.0,
                                        GENEMouseOverArcStrokeColor: "#f39c12",
                                        GENEMouseOverArcStrokeWidth: 2,
                                        
                                        // Fixed tooltip configurations to match the data structure
                                        GENEMouseOverTooltipsHtml01: "<strong>Contig:</strong> ",  // chr
                                        GENEMouseOverTooltipsHtml02: "<br><strong>Strand:</strong> ", // strand
                                        GENEMouseOverTooltipsHtml03: "<br><strong>Start:</strong> ", // start
                                        GENEMouseOverTooltipsHtml04: "<br><strong>End:</strong> ", // end
                                        GENEMouseOverTooltipsHtml05: "<br><strong>Type:</strong> ", // type
                                        GENEMouseOverTooltipsHtml06: "<br><strong>Name:</strong> ", // name
                                        GENEMouseOverTooltipsHtml07: "<br><strong>Product:</strong> ", // product
                                        GENEMouseOverTooltipsHtml08: "",

                                        GENEMouseOverTooltipsPosition: "absolute",
                                        GENEMouseOverTooltipsBackgroundColor: "rgba(250, 250, 250, 0.95)",
                                        GENEMouseOverTooltipsBorderStyle: "solid",
                                        GENEMouseOverTooltipsBorderWidth: 1,
                                        GENEMouseOverTooltipsPadding: "10px",
                                        GENEMouseOverTooltipsBorderRadius: "5px",
                                        GENEMouseOverTooltipsOpacity: 0.95,
                                        GENEMouseOverTooltipsBorderColor: "#3498db",
                                        GENEMouseOverTooltipsBoxShadow: "2px 2px 4px rgba(0, 0, 0, 0.2)",
                                    });
                                    
                                    // Draw genome and tracks
                                    NGCircos1.draw_genome(NGCircos1.genomeLength);
                                    
                                    showDebugInfo('Visualization complete!');
                                } catch (error) {
                                    showDebugInfo('Error creating plot: ' + error.message);
                                    console.error('Error details:', error);
                                }
                            }
                        };
                    </script>
                </body>
                </html>
            """)
            logger.info("Created circos plot HTML content successfully")
            print(html_content)
            return html_content
        except Exception as e:
            logger.error(f"Error creating circos plot: {str(e)}", exc_info=True)
            return f"""
                <html>
                <body>
                    <div style="color: red; padding: 20px; border: 1px solid #ccc;">
                        <h3>Error creating plot</h3>
                        <p>{str(e)}</p>
                        <pre>{traceback.format_exc()}</pre>
                    </div>
                </body>
                </html>
            """

    def _create_sequence_viewer_content(self, sequence_details):
        """Create content for the sequence viewer with enhanced features."""
        logger.info("Creating sequence viewer content")
        
        if not sequence_details:
            logger.warning("No sequence details provided to sequence viewer")
            return html.Div("No sequence details available", 
                           style={'padding': '20px', 'textAlign': 'center'})
        
        sequence = sequence_details.get('sequence', '')
        header = sequence_details.get('header', '')
        seq_id = sequence_details.get('id', '')
        length = sequence_details.get('length', 0)
        gc_content = sequence_details.get('gc_content', 0)
        
        logger.info(f"Creating viewer for sequence: id={seq_id}, length={length}, gc={gc_content:.1f}%")
        logger.debug(f"Sequence header: {header}")
        logger.debug(f"Sequence first 50 bases: {sequence[:50]}...")
        
        try:
            # Format sequence with line breaks and position numbers
            formatted_sequence = []
            line_length = 60
            
            # Limit to first 10 lines for performance
            max_lines = 10
            total_lines = (len(sequence) + line_length - 1) // line_length  # Ceiling division
            
            for i in range(0, min(max_lines * line_length, len(sequence)), line_length):
                line = sequence[i:i+line_length]
                line_num = i + 1
                formatted_sequence.append(
                    html.Div([
                        html.Span(f"{line_num:>8} ", 
                                 style={'fontFamily': 'monospace', 'color': '#666', 'marginRight': '10px'}),
                        html.Span(line, 
                                 style={'fontFamily': 'monospace', 'letterSpacing': '0.5px'})
                    ], style={'marginBottom': '5px'})
                )
            
            # Add a "more lines" indicator if we truncated the sequence
            if total_lines > max_lines:
                formatted_sequence.append(
                    html.Div([
                        html.Span(f"       ", style={'fontFamily': 'monospace', 'color': '#666', 'marginRight': '10px'}),
                        html.Span(f"... {total_lines - max_lines} more lines not shown (use Copy buttons to get full sequence)",
                                 style={'fontStyle': 'italic', 'color': '#666'})
                    ], style={'marginBottom': '5px', 'marginTop': '10px'})
                )
            
            logger.info(f"Created {len(formatted_sequence)} formatted sequence lines (showing {max_lines} of {total_lines} total)")
            
            # Create a hidden div with the full sequence data for copying
            full_sequence_div = html.Div(sequence, id='sequence-data', style={'display': 'none'})
            full_fasta_div = html.Div(f">{header}\n{sequence}", id='fasta-data', style={'display': 'none'})
            
            return html.Div([
                # Hidden data for copy operations
                full_sequence_div,
                full_fasta_div,
                
                # Sequence header and metadata
                html.Div([
                    html.H5(header, style={'wordBreak': 'break-all', 'marginBottom': '15px'}),
                    html.Div([
                        html.Span(f"Length: {length:,} bp", style={'marginRight': '20px'}),
                        html.Span(f"GC Content: {gc_content:.1f}%")
                    ], style={'marginBottom': '15px'})
                ], style={'borderBottom': '1px solid #ddd', 'paddingBottom': '10px', 'marginBottom': '15px'}),
                
                # Action buttons
                html.Div([
                    dmc.Button(
                        "Copy Sequence",
                        id="copy-sequence-btn",
                        leftSection=DashIconify(icon="mdi:content-copy"),
                        style={'marginRight': '10px'}
                    ),
                    dmc.Button(
                        "Copy FASTA",
                        id="copy-fasta-btn",
                        leftSection=DashIconify(icon="mdi:dna"),
                        variant="outline"
                    ),
                    # Clipboard success message
                    dmc.Alert(
                        "Copied to clipboard!",
                        id="copy-alert",
                        color="green",
                        withCloseButton=True,
                        style={'display': 'none', 'marginTop': '10px'}
                    )
                ], style={'marginBottom': '20px'}),
                
                # Preview notice
                html.Div([
                    html.P("Preview only - showing first 10 lines. Use copy buttons for full sequence.", 
                          style={'color': '#666', 'fontStyle': 'italic', 'marginBottom': '10px'})
                ]),
                
                # Sequence content with highlighting
                html.Div([
                    html.H6("Sequence:", style={'marginBottom': '10px'}),
                    html.Div(
                        formatted_sequence,
                        id='sequence-text',
                        style={
                            'backgroundColor': '#f8f9fa',
                            'padding': '15px',
                            'border': '1px solid #ddd',
                            'borderRadius': '4px',
                            'overflowX': 'auto',
                            'maxHeight': '400px',
                            'overflowY': 'auto',
                            'fontFamily': 'monospace'
                        }
                    )
                ])
            ])
        except Exception as e:
            logger.error(f"Error creating sequence viewer content: {str(e)}", exc_info=True)
            return html.Div([
                html.H5("Error Rendering Sequence"),
                html.P(f"An error occurred while rendering the sequence: {str(e)}"),
                html.Hr(),
                html.P("Sequence details:"),
                html.Pre(f"ID: {seq_id}\nLength: {length}\nHeader: {header[:100]}...")
            ], className="alert alert-danger")

    def _get_sequence_details(self, sequence_id, mag_id=None):
        """Get details for a specific sequence.
        
        Args:
            sequence_id (str): ID of the sequence to retrieve
            mag_id (str, optional): ID of the MAG containing the sequence. If provided, 
                                    only this MAG's data will be loaded.
        
        Returns:
            dict or None: Sequence details or None if not found
        """
        print(f"\n===== Getting Sequence Details for {sequence_id} =====")
        
        try:
            if not sequence_id:
                print("No sequence_id provided")
                return None
            
            # If we know which MAG to check, optimize by only loading that MAG's data
            if mag_id:
                print(f"Loading only data for MAG {mag_id}")
                # Get MAG data with fasta content, but only for the specified MAG
                mag_data = {}
                try:
                    from users.models import MAG
                    mag = MAG.objects.get(name=mag_id, user_id=self.user_id)
                    
                    if mag.fasta_file:
                        # Parse the MAG's FASTA file
                        sequences = self.parse_fasta_for_viewer(mag.fasta_file)
                        
                        # Create MAG data entry
                        mag_data[mag_id] = {
                            'id': mag.id,
                            'name': mag.name,
                            'bin_id': mag.name,
                            'taxonomy': mag.taxonomy or {},
                            'completeness': mag.completeness or 0,
                            'contamination': mag.contamination or 0,
                            'fasta': sequences,
                            'fasta_file': mag.fasta_file
                        }
                        
                        print(f"Loaded single MAG {mag_id} with {len(sequences)} sequences")
                    else:
                        print(f"MAG {mag_id} has no FASTA file")
                        return None
                except MAG.DoesNotExist:
                    print(f"MAG {mag_id} not found in database")
                    return None
                except Exception as e:
                    print(f"Error loading single MAG {mag_id}: {str(e)}")
                    print(traceback.format_exc())
                    return None
            else:
                # Get all MAG data with fasta content (this is slower)
                print("Fetching MAG data with FASTA content for all MAGs")
                mag_data = self.get_mag_data(include_fasta=True)
            
            if not mag_data:
                print("No MAG data available")
                return None
            
            print(f"MAG data contains {len(mag_data)} entries")
            # Log a sample of MAG keys
            if len(mag_data) > 0:
                sample_keys = list(mag_data.keys())[:3]
                print(f"Sample MAG keys: {sample_keys}")
            
            # Handle the case where sequence_id is 'complete_sequence' (fallback option)
            if sequence_id == 'complete_sequence' or sequence_id == 'sequence_1' or sequence_id == 'emergency_sequence':
                print(f"Handling fallback sequence request: {sequence_id}")
                # Find first MAG with fasta_file content
                for mag_id, mag in mag_data.items():
                    if isinstance(mag, dict) and 'fasta_file' in mag and mag['fasta_file']:
                        fasta_content = mag['fasta_file']
                        # Create a sequence with the entire content
                        cleaned_sequence = fasta_content
                        if '>' in cleaned_sequence:
                            # Remove FASTA headers to get only sequence data
                            parts = []
                            for line in fasta_content.split('\n'):
                                if not line.startswith('>'):
                                    parts.append(line.strip())
                            cleaned_sequence = ''.join(parts)
                        
                        # Create sequence details
                        return {
                            'id': sequence_id,
                            'header': f'>{sequence_id}',
                            'sequence': cleaned_sequence,
                            'length': len(cleaned_sequence),
                            'gc_content': self._calculate_gc_content_for_sequence(cleaned_sequence)
                        }
        
            # Search through all MAGs for this sequence ID
            for mag_id, mag in mag_data.items():
                if not isinstance(mag, dict):
                    continue
                    
                if 'fasta' in mag and isinstance(mag['fasta'], list):
                    for seq in mag['fasta']:
                        if seq.get('id') == sequence_id:
                            sequence_info = seq
                            parent_mag = mag
                            print(f"Found sequence in MAG {mag.get('bin_id', mag_id)}")
                            break
                    
                    if 'sequence_info' in locals() and sequence_info:
                        break
                        
                elif 'fasta_file' in mag and mag['fasta_file']:
                    # Parse FASTA content
                    sequences = self.parse_fasta_for_viewer(mag['fasta_file'])
                    for seq in sequences:
                        if seq.get('id') == sequence_id:
                            sequence_info = seq
                            parent_mag = mag
                            print(f"Found sequence in parsed FASTA for MAG {mag.get('bin_id', mag_id)}")
                            break
                    
                    if 'sequence_info' in locals() and sequence_info:
                        break
        
            if 'sequence_info' not in locals() or not sequence_info:
                print(f"Sequence {sequence_id} not found in any MAG")
                return None
        
            # Calculate GC content if not already present
            if 'gc_content' not in sequence_info and 'sequence' in sequence_info:
                sequence = sequence_info['sequence']
                if sequence:
                    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
                    total = len(sequence)
                    sequence_info['gc_content'] = (gc_count / total) * 100 if total > 0 else 0
                    print(f"Calculated GC content: {sequence_info['gc_content']:.4f}")
            
            # Add parent MAG info to sequence details
            if 'parent_mag' in locals() and parent_mag:
                sequence_info['parent_mag'] = parent_mag.get('bin_id', parent_mag.get('name', 'Unknown'))
                sequence_info['taxonomy'] = parent_mag.get('taxonomy', 'Unknown')
            
            print(f"Returning sequence details with keys: {list(sequence_info.keys())}")
            return sequence_info
        
        except Exception as e:
            print(f"Error getting sequence details: {str(e)}")
            print(traceback.format_exc())
            return None

    def _create_placeholder_annotations(self, sequences):
        """Create placeholder gene annotations for visualization when no GFF data is available."""
        annotations = []
        
        print("Creating placeholder annotations for sequences")
        
        for seq in sequences:
            seq_id = seq.get('id', '')
            seq_len = seq.get('length', 0)
            
            if not seq_id or not seq_len:
                continue
                
            # Create evenly distributed gene features along the sequence
            gene_count = min(20, max(5, int(seq_len / 5000)))  # 1 gene per ~5kb but at least 5, max 20
            interval = seq_len / gene_count
            
            for i in range(gene_count):
                # Calculate gene position
                pos = int(i * interval)
                gene_len = min(1500, int(interval * 0.7))  # Gene length ~70% of interval, max 1.5kb
                
                # Alternate gene features between + and - strand
                strand = '+' if i % 2 == 0 else '-'
                
                # Determine gene type - mostly CDS with some RNA genes
                feature_type = 'cds'
                if i % 10 == 0:
                    feature_type = 'rRNA'
                elif i % 7 == 0:
                    feature_type = 'tRNA'
                
                # Map to NGCircos types
                if feature_type in ['rRNA', 'tRNA']:
                    ng_type = 'utr'  # RNA features are 'utr' type in NGCircos
                else:
                    ng_type = 'cds'  # Protein-coding genes are 'cds' type
                
                # Set color based on feature type
                color = '#1f77b4'  # Blue for CDS
                if feature_type == 'rRNA':
                    color = '#2ca02c'  # Green for rRNA
                elif feature_type == 'tRNA':
                    color = '#d62728'  # Red for tRNA
                
                # Create gene annotation
                gene_name = f"gene_{i+1}"
                function = f"Hypothetical protein - placeholder gene {i+1}"
                if feature_type == 'tRNA':
                    function = f"Transfer RNA for amino acid transport - placeholder {i+1}"
                elif feature_type == 'rRNA':
                    function = f"Ribosomal RNA {['16S', '23S', '5S'][i % 3]} - structural RNA - placeholder {i+1}"
                
                # Create annotation in the format expected by NGCircos
                annotations.append({
                    'contig': seq_id,
                    'start': pos,
                    'end': pos + gene_len,
                    'strand': strand,
                    'type': ng_type,  # Use NGCircos type
                    'length': gene_len,
                    'product': function,
                    'gene': gene_name,
                    'locus_tag': f"loc_{seq_id}_{i+1}",
                    'name': gene_name,
                    'color': color,
                    'function': function,
                    'original_type': feature_type  # Keep original type for reference
                })
            
            # Add an operon - genes close together on same strand
            operon_start = int(seq_len * 0.7)
            gene_gap = 100
            genes_in_operon = 4
            for j in range(genes_in_operon):
                gene_len = 600
                pos = operon_start + (j * (gene_len + gene_gap))
                if pos + gene_len < seq_len:
                    annotations.append({
                        'contig': seq_id,
                        'start': pos,
                        'end': pos + gene_len,
                        'strand': '+',  # All in same orientation for operon
                        'type': 'cds',
                        'length': gene_len,
                        'product': f"Operon protein {j+1} - related function",
                        'gene': f"operon_gene_{j+1}",
                        'locus_tag': f"op_{seq_id}_{j+1}",
                        'name': f"op_gene_{j+1}",
                        'color': '#ff7f0e',  # Orange for operon genes
                        'function': f"Operon protein {j+1} - related function",
                        'original_type': 'CDS'
                    })
        
        print(f"Created {len(annotations)} placeholder annotations across {len(sequences)} sequences")
        return annotations

    def _parse_gff_for_annotations(self, gff_content, sequences):
        """Parse GFF file to extract gene annotations for visualization."""
        if not gff_content:
            return []
        
        print("\n===== Parsing GFF for Annotations =====")
        
        # Get FASTA sequence IDs
        fasta_ids = [seq.get('id', '') for seq in sequences]
        print(f"FASTA sequence IDs: {', '.join(fasta_ids[:5])}...")
        print(f"Total sequences in FASTA: {len(sequences)}")
        
        # Create a mapping from GFF contig IDs to FASTA IDs
        contig_map = self._create_contig_mapping(gff_content, sequences)
        
        # Process the GFF file with our mapping
        annotations = []
        line_count = 0
        feature_count = 0
        feature_types = set()
        
        try:
            # First, extract the GFF features
            gff_features = []
            for line in gff_content.split('\n'):
                line_count += 1
                
                # Skip comments and empty lines
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.split('\t')
                if len(parts) < 8:
                    continue
                
                gff_contig = parts[0]
                feature_type = parts[2].lower()
                
                # Map GFF contig to FASTA contig
                fasta_contig = contig_map.get(gff_contig)
                if not fasta_contig:
                    continue
                
                feature_types.add(feature_type)
                
                # Extract common feature attributes
                try:
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = parts[6]
                    length = end - start + 1
                    
                    # Skip very large regions (likely whole contigs)
                    if feature_type == 'region' and length > 100000:
                        continue
                    
                    # Parse attributes
                    attributes = {}
                    if len(parts) > 8:
                        for attr in parts[8].split(';'):
                            if '=' in attr:
                                key, value = attr.split('=', 1)
                                attributes[key] = value.replace('%2C', ',').replace('%2F', '/').replace('%20', ' ')
                    
                    # Store the feature for later processing
                    gff_features.append({
                        'contig': fasta_contig,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'type': feature_type,
                        'length': length,
                        'attributes': attributes
                    })
                    
                except Exception as e:
                    print(f"Error parsing GFF line {line_count}: {str(e)}")
                    continue
            
            print(f"Extracted {len(gff_features)} features from GFF file")
            
            # Convert GFF features to NGCircos annotation format
            for feature in gff_features:
                # Only include features we care about
                if feature['type'] not in ['cds', 'gene', 'trna', 'rrna', 'ncRNA', 'repeat_region', 
                                      'tmrna', 'srna', 'oric', 'oriv', 'crispr']:
                    continue
                
                # Get feature details
                attributes = feature['attributes']
                
                # Get identifiers from various attribute fields
                gene_id = attributes.get('ID', '')
                locus_tag = attributes.get('locus_tag', '')
                gene_name = attributes.get('gene', '')
                product = attributes.get('product', '') or attributes.get('Name', '')
                function = attributes.get('function', '') or attributes.get('description', '') or ''
                
                # Determine feature name (prioritize informative identifiers)
                feature_name = None
                for key in ['gene', 'locus_tag', 'ID', 'Name']:
                    if key in attributes and attributes[key]:
                        feature_name = attributes[key]
                        break
                
                if not feature_name:
                    feature_name = f"{feature['type']}_{feature['start']}"
                
                # Standardize feature type for NGCircos
                ng_type = feature['type']
                if ng_type in ['trna', 'rrna', 'srna', 'ncRNA', 'tmrna']:
                    ng_type = 'utr'  # RNA features displayed as 'utr' in NGCircos
                elif ng_type in ['cds', 'gene']:
                    ng_type = 'cds'  # Protein-coding genes as 'cds'
                
                # Combine function descriptions for tooltip
                function_desc = ''
                if product:
                    function_desc = product
                if function and function != product:
                    function_desc = f"{function_desc} - {function}" if function_desc else function
                if gene_name and gene_name not in function_desc:
                    function_desc = f"{function_desc} ({gene_name})" if function_desc else gene_name
                    
                if not function_desc:
                    function_desc = gene_id or locus_tag or feature_name or feature['type']
                
                # Add annotation
                annotation = {
                    'contig': feature['contig'],
                    'start': feature['start'],
                    'end': feature['end'],
                    'strand': '+' if feature['strand'] == '+' else '-',
                    'type': ng_type,
                    'length': feature['length'],
                    'product': function_desc,
                    'gene': gene_name or locus_tag or feature_name,
                    'locus_tag': locus_tag,
                    'name': feature_name,
                    'color': self._get_feature_color(feature['type']),
                    'function': function_desc
                }
                annotations.append(annotation)
                feature_count += 1
                
                # Log progress
                if feature_count % 500 == 0:
                    print(f"Processed {feature_count} annotations...")
            
            print(f"Created {len(annotations)} annotations from {len(gff_features)} GFF features")
            
            # Only limit features if absolutely necessary for performance
            max_features = 50000  # Much larger limit
            if len(annotations) > max_features:
                print(f"Warning: Too many features ({len(annotations)}), limiting to {max_features} for performance")
                
                # Instead of sorting by length (which might exclude distant genes),
                # take an evenly distributed sample across the genome
                step = len(annotations) // max_features
                if step < 1:
                    step = 1
                
                # Get a sampling that preserves position distribution
                sampled_annotations = []
                for i in range(0, len(annotations), step):
                    sampled_annotations.append(annotations[i])
                
                # Add a few from the end to ensure we have representation of all regions
                if len(sampled_annotations) < max_features and len(annotations) > max_features:
                    remaining = max_features - len(sampled_annotations)
                    end_sample = annotations[-remaining:]
                    sampled_annotations.extend(end_sample)
                
                annotations = sampled_annotations
                print(f"Sampled down to {len(annotations)} annotations")
            
            return annotations
            
        except Exception as e:
            print(f"Error parsing GFF: {str(e)}")
            print(traceback.format_exc())
            return []

    def _get_feature_color(self, feature_type):
        """Get color for different feature types."""
        # Convert to lowercase for consistent comparison
        feature_type = feature_type.lower()
        
        color_map = {
            'cds': '#1f77b4',      # Blue
            'gene': '#ff7f0e',     # Orange
            'rrna': '#2ca02c',     # Green
            'trna': '#d62728',     # Red
            'ncrna': '#9467bd',    # Purple
            'repeat_region': '#8c564b',   # Brown
            'mobile_element': '#e377c2',  # Pink
            'tmrna': '#7f7f7f',    # Gray
            'srp_rna': '#bcbd22',  # Olive
            'oric': '#17becf',     # Cyan
            'oriv': '#17becf',     # Cyan
            'crispr': '#e6ab02',   # Gold
            'srna': '#a6761d',     # Dark orange
            'ncrna': '#9467bd',    # Purple
            'exon': '#ff7f0e',     # Orange
            'mirna': '#9edae5',    # Light blue
            'pseudogene': '#c49c94' # Light brown
        }
        
        # Standardize common feature type variations
        if 'rna' in feature_type:
            if 'rrna' in feature_type:
                return color_map['rrna']
            elif 'trna' in feature_type:
                return color_map['trna']
            elif 'ncrna' in feature_type or 'nc_rna' in feature_type:
                return color_map['ncrna']
            elif 'srna' in feature_type or 'small' in feature_type:
                return color_map['srna']
            elif 'tmrna' in feature_type:
                return color_map['tmrna']
            else:
                return color_map['ncrna']  # Default for any RNA
        
        # Return the color if we have it, otherwise use gray as default
        return color_map.get(feature_type, '#7f7f7f')
    
    def _calculate_gc_skew(self, sequence, window_size=1000):
        """Calculate GC skew for a DNA sequence."""
        if not sequence:
            return []
        
        gc_skew = []
        sequence = sequence.upper()
        for i in range(0, len(sequence), window_size):
            window = sequence[i:i+window_size]
            g_count = window.count('G')
            c_count = window.count('C')
            
            if g_count + c_count > 0:
                skew = (g_count - c_count) / (g_count + c_count)
            else:
                skew = 0
                
            gc_skew.append({
                'position': i,
                'skew': skew
            })
            
        return gc_skew
    
    def _calculate_gc_content_for_sequence(self, sequence):
        """Calculate GC content for a given sequence."""
        if not sequence:
            return 0.0
            
        sequence = sequence.upper()
        gc_count = sequence.count('G') + sequence.count('C')
        total_bases = len(sequence)
        return (gc_count / total_bases) * 100 if total_bases > 0 else 0

    def _create_kegg_overview(self, highlight_mag=None):
        """Create a KEGG pathway overview plot."""
        try:
            if not self.mag_data:
                logger.warning("No MAG data available for KEGG overview")
                return go.Figure()

            # Get KEGG pathways from annotation dictionary
            self.mag_data['kegg_pathways'] = self.mag_data['annotations_dict'].apply(
                lambda x: [pathway for pathway in x.get('kegg_pathways', []) if isinstance(pathway, str)] if isinstance(x, dict) else []
            )
            
            # Flatten the list of pathways
            all_pathways = [pathway for pathways in self.mag_data['kegg_pathways'] for pathway in pathways]
            
            # Count pathways
            pathway_counts = {}
            for pathway in all_pathways:
                pathway_counts[pathway] = pathway_counts.get(pathway, 0) + 1
            
            # Create bar plot
            fig = go.Figure()
            
            fig.add_trace(go.Bar(
                x=list(pathway_counts.keys()),
                y=list(pathway_counts.values()),
                marker_color='#1f77b4',
                hovertemplate='<b>%{x}</b><br>Pathway count: %{y}<extra></extra>'
            ))
            
            # Update layout
            fig.update_layout(
                title='KEGG Pathway Overview',
                xaxis_title='KEGG Pathway',
                yaxis_title='Count',
                height=500,
                margin=dict(l=10, r=10, t=40, b=10),
                yaxis={'categoryorder': 'total ascending'},
                showlegend=False
            )
            
            return fig
        except Exception as e:
            logger.error(f"Error creating KEGG overview plot: {str(e)}", exc_info=True)
            return go.Figure()

        def sync_table_data(data, page_size):
            """Sync the original table data with our DMC table."""
            if not data:
                return [], 1
                
            # Calculate total pages
            page_size = int(page_size) if page_size else 10
            total_pages = (len(data) + page_size - 1) // page_size
            
            return data, total_pages

   
    def _create_contig_mapping(self, gff_content, sequences):
        """Create a mapping from GFF contig IDs to FASTA contig IDs."""
        if not gff_content or not sequences:
            return {}
            
        # Get FASTA sequence IDs and lengths
        fasta_seq_info = {seq.get('id', ''): seq.get('length', 0) for seq in sequences if seq.get('id')}
        fasta_ids = list(fasta_seq_info.keys())
        
        # Extract GFF contig IDs
        gff_contig_ids = set()
        for line in gff_content.split('\n'):
            if not line.startswith('#') and '\t' in line:
                parts = line.split('\t')
                if len(parts) > 0:
                    gff_contig_ids.add(parts[0])
        
        print(f"GFF contig IDs found: {', '.join(list(gff_contig_ids)[:10])}{'...' if len(gff_contig_ids) > 10 else ''}")
        
        # Create mapping dictionary
        contig_map = {}
        
        # Check for Bakta-style renaming (contig_1, contig_2, etc.)
        bakta_style_contigs = [cid for cid in gff_contig_ids if cid.startswith('contig_') and cid.split('_')[1].isdigit()]
        is_bakta_style = len(bakta_style_contigs) > 0
        
        if is_bakta_style and len(sequences) > 0:
            print("Detected Bakta-style contig naming - using positional mapping")
            
            # Sort FASTA sequences by length (largest first)
            sorted_fasta = sorted(sequences, key=lambda x: x.get('length', 0), reverse=True)
            sorted_fasta_ids = [seq.get('id', '') for seq in sorted_fasta]
            
            # Create mapping assuming gff contig_N maps to Nth largest fasta contig
            for i in range(1, len(sorted_fasta_ids) + 1):
                gff_id = f"contig_{i}"
                if i <= len(sorted_fasta_ids):
                    fasta_id = sorted_fasta_ids[i-1]
                    contig_map[gff_id] = fasta_id
                    print(f"Mapping {gff_id} → {fasta_id}")
        else:
            # Try direct matching first
            for gff_id in gff_contig_ids:
                if gff_id in fasta_ids:
                    contig_map[gff_id] = gff_id
            
            # If few matches, try more advanced mapping
            if len(contig_map) < min(5, len(gff_contig_ids) // 2):
                print("Few direct matches found, trying advanced mapping")
                
                # Try case-insensitive matching
                for gff_id in gff_contig_ids:
                    gff_id_lower = gff_id.lower()
                    for fasta_id in fasta_ids:
                        if fasta_id.lower() == gff_id_lower:
                            contig_map[gff_id] = fasta_id
                
                # Special case for single contig
                if len(sequences) == 1 and len(gff_contig_ids) > 0:
                    print("Single contig case - mapping all GFF contigs to the single FASTA sequence")
                    fasta_id = fasta_ids[0]
                    for gff_id in gff_contig_ids:
                        if gff_id not in contig_map:
                            contig_map[gff_id] = fasta_id
        
        print(f"Created {len(contig_map)} contig mappings for {len(gff_contig_ids)} GFF contigs")
        return contig_map

    
    def _calculate_gc_content_windows(self, sequence, window_size=1000, step=500):
        """Calculate GC content for a DNA sequence using a sliding window."""
        if not sequence or len(sequence) < window_size:
            return []
            
        gc_content_data = []
        sequence = sequence.upper()
        
        for i in range(0, len(sequence) - window_size + 1, step):
            window = sequence[i : i + window_size]
            g_count = window.count('G')
            c_count = window.count('C')
            gc_total = g_count + c_count
            total_bases = len(window)
            
            if total_bases > 0:
                gc_percent = (gc_total / total_bases) * 100
                position = i + window_size // 2  # Use midpoint of window
                gc_content_data.append((position, gc_percent))
                
        return gc_content_data

    def _calculate_gc_skew_windows(self, sequence, window_size=1000, step=500):
        """Calculate GC skew (G-C)/(G+C) for a DNA sequence using a sliding window.
        
        Args:
            sequence: DNA sequence string
            window_size: Size of the sliding window
            step: Step size for the sliding window
            
        Returns:
            List of tuples (position, skew_value) where position is the midpoint of the window
        """
        if not sequence or len(sequence) < window_size:
            return []
            
        gc_skew_data = []
        sequence = sequence.upper()
        
        for i in range(0, len(sequence) - window_size + 1, step):
            window = sequence[i : i + window_size]
            g_count = window.count('G')
            c_count = window.count('C')
            
            # Calculate skew only if we have G or C bases
            if g_count + c_count > 0:
                # Calculate GC skew: (G-C)/(G+C)
                skew = (g_count - c_count) / (g_count + c_count)
                # Use midpoint of window for position
                position = i + window_size // 2
                gc_skew_data.append((position, skew))
            else:
                # If no G or C bases, use 0 skew
                position = i + window_size // 2
                gc_skew_data.append((position, 0))
                
        return gc_skew_data

    def _sanitize_text(self, text):
        """Sanitize text to be safe for JavaScript strings."""
        if not text:
            return ""
        # Replace problematic characters with safe alternatives
        text = str(text)
        # Replace quotes and backslashes for JS string safety
        text = text.replace('\\', '\\\\').replace('"', '\\"').replace("'", "\\'")
        # Replace other potentially problematic characters
        text = text.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
        # Handle other special characters by replacing with ASCII equivalents or removing
        # Map common non-ASCII characters to ASCII equivalents
        char_map = {
            '⊗': 'x',  # Replace with 'x'
            '→': '->',  # Right arrow
            '←': '<-',  # Left arrow
            '↔': '<->',  # Bidirectional arrow
            '≥': '>=',  # Greater than or equal
            '≤': '<=',  # Less than or equal
            '±': '+/-', # Plus-minus
            '°': 'deg', # Degree
            '×': 'x',   # Multiplication
            '÷': '/',   # Division
            '≈': '~',   # Approximately equal
            '≠': '!=',  # Not equal
            '∞': 'inf', # Infinity
            '²': '2',   # Superscript 2
            '³': '3',   # Superscript 3
            'α': 'alpha',
            'β': 'beta',
            'γ': 'gamma',
            'δ': 'delta',
            'ε': 'epsilon',
            'µ': 'mu',
            'π': 'pi',
            'σ': 'sigma',
            'τ': 'tau',
            'φ': 'phi'
        }
        for char, replacement in char_map.items():
            text = text.replace(char, replacement)
        
        # Replace any remaining non-ASCII characters with underscores
        text = ''.join(c if ord(c) < 128 else '_' for c in text)
        return text

    def _create_cgview_plot(self, mag_id=None):
        """Create a CGView.js circular genome visualization.
        
        Args:
            mag_id (str): ID of the MAG to visualize
            
        Returns:
            str: HTML content with CGView.js visualization
        """
        try:
            # Clear logging to help debug
            logger.info(f"Creating CGView plot for MAG ID: {mag_id if mag_id else 'None'}")
            
            # If no MAG ID is provided, don't try to create a plot with data
            if not mag_id:
                logger.info("No MAG ID provided, returning placeholder message")
                return """
                <html>
                <body style="text-align: center; font-family: Arial, sans-serif; background-color: #f9f9f9;">
                    <div style="padding: 40px 20px;">
                        <div style="font-size: 18px; margin-bottom: 15px;">Select a MAG from the table to view genome visualization</div>
                        <div style="color: #666; font-size: 14px;">The circular plot will show contigs and gene annotations when a MAG is selected</div>
                    </div>
                </body>
                </html>
                """
            
            # Get MAG data including sequence information
            mag_info = None
            sequences = []
            
            # Force refresh the mag_data to avoid using cached data
            mag_data = self.get_mag_data(include_fasta=True, force_refresh=True)
            
            try:
                # Get the MAG data for the specific MAG ID
                logger.info(f"Retrieving MAG data for {mag_id}")
                
                if not mag_data:
                    logger.warning("No MAG data returned from get_mag_data")
                    return """
                    <html>
                    <body style="text-align: center; font-family: Arial, sans-serif; background-color: #f9f9f9;">
                        <div style="padding: 40px 20px;">
                            <div style="font-size: 18px; color: #e74c3c; margin-bottom: 15px;">No MAG data available</div>
                            <div style="color: #666; font-size: 14px;">Please check database connection or data availability</div>
                        </div>
                    </body>
                    </html>
                    """
                    
                if mag_id in mag_data:
                    mag_info = mag_data[mag_id]
                    logger.info(f"Found MAG info for {mag_id}")
                else:
                    # Try case-insensitive match if not found
                    mag_id_lower = str(mag_id).lower()
                    for key, info in mag_data.items():
                        if isinstance(info, dict) and str(info.get('bin_id', '')).lower() == mag_id_lower:
                            mag_info = info
                            logger.info(f"Found MAG info for {mag_id} via case-insensitive match")
                            break
                    
                    if not mag_info:
                        logger.warning(f"MAG ID {mag_id} not found in MAG data dictionary")
                        return f"""
                        <html>
                        <body style="text-align: center; font-family: Arial, sans-serif; background-color: #f9f9f9;">
                            <div style="padding: 40px 20px;">
                                <div style="font-size: 18px; color: #e74c3c; margin-bottom: 15px;">MAG ID {mag_id} not found</div>
                                <div style="color: #666; font-size: 14px;">Please select a different MAG</div>
                            </div>
                        </body>
                        </html>
                        """
                    
                # Get sequences if we found the MAG
                if mag_info:
                    # Get sequences
                    if 'fasta' in mag_info and isinstance(mag_info['fasta'], list):
                        sequences = mag_info['fasta']
                        logger.info(f"Found {len(sequences)} parsed sequences for MAG {mag_id}")
                    elif 'fasta_file' in mag_info and mag_info['fasta_file']:
                        sequences = self.parse_fasta_for_viewer(mag_info['fasta_file'])
                        logger.info(f"Parsed {len(sequences)} sequences from FASTA file for MAG {mag_id}")
                    else:
                        logger.warning(f"No FASTA data found for MAG {mag_id}")
                        return f"""
                        <html>
                        <body style="text-align: center; font-family: Arial, sans-serif; background-color: #f9f9f9;">
                            <div style="padding: 40px 20px;">
                                <div style="font-size: 18px; color: #e74c3c; margin-bottom: 15px;">No sequence data found for MAG {mag_id}</div>
                                <div style="color: #666; font-size: 14px;">This MAG doesn't have any FASTA sequence data</div>
                            </div>
                        </body>
                        </html>
                        """
            except Exception as e:
                logger.error(f"Error retrieving MAG data: {str(e)}", exc_info=True)
                return f"""
                <html>
                <body style="text-align: center; font-family: Arial, sans-serif; background-color: #f9f9f9;">
                    <div style="padding: 40px 20px;">
                        <div style="font-size: 18px; color: #e74c3c; margin-bottom: 15px;">Error retrieving MAG data</div>
                        <div style="color: #666; font-size: 14px;">Error details: {str(e)}</div>
                    </div>
                </body>
                </html>
                """
        
            # Try to get GFF data for annotations
            gene_annotations = []
            try:
                logger.info(f"Retrieving GFF data for {mag_id}")
                
                # First check if GFF data is directly available in the mag_info
                if mag_info and 'gff_file' in mag_info and mag_info['gff_file']:
                    logger.info(f"Found GFF data in mag_info for {mag_id}, parsing...")
                    gff_content = mag_info['gff_file']
                    gene_annotations = self._parse_gff_for_annotations(gff_content, sequences)
                    logger.info(f"Parsed {len(gene_annotations)} annotations from GFF data")
                else:
                    # Try the database method
                    logger.info(f"Using get_gff_data method for MAG {mag_id}")
                    gff_data = self.get_gff_data(mag_id)
                    if gff_data:
                        logger.info(f"Found GFF data for MAG {mag_id}")
                        if isinstance(gff_data, str):
                            # Parse the GFF content if it's a string
                            gene_annotations = self._parse_gff_for_annotations(gff_data, sequences)
                            logger.info(f"Parsed {len(gene_annotations)} annotations from GFF string")
                        elif isinstance(gff_data, list):
                            # GFF data is already a list of feature dictionaries
                            gene_annotations = gff_data
                            logger.info(f"Using {len(gene_annotations)} annotations from GFF list data")
                        else:
                            logger.warning(f"Unexpected GFF format for MAG {mag_id}: {type(gff_data)}")
                    else:
                        logger.warning(f"No GFF data found for MAG {mag_id}")
                
                # If no annotations were found, try direct GFF parsing from database
                if not gene_annotations:
                    try:
                        from users.models import MAG
                        logger.info(f"Trying to get GFF directly from database for MAG {mag_id}")
                        mag_obj = MAG.objects.get(name=mag_id, user_id=self.user_id)
                        
                        if mag_obj.gff_file:
                            logger.info(f"Found GFF in database for MAG {mag_id}, parsing...")
                            gene_annotations = self._parse_gff_for_annotations(mag_obj.gff_file, sequences)
                            logger.info(f"Parsed {len(gene_annotations)} annotations directly from database GFF")
                        else:
                            logger.warning(f"No GFF file found in database for MAG {mag_id}")
                    except Exception as db_error:
                        logger.error(f"Error retrieving GFF from database: {str(db_error)}")
                
                # As a last resort, create placeholder annotations
                if not gene_annotations and sequences:
                    logger.info(f"No annotations found for MAG {mag_id}, creating placeholder annotations")
                    gene_annotations = self._create_placeholder_annotations(sequences)
                    logger.info(f"Created {len(gene_annotations)} placeholder annotations")
                    
            except Exception as e:
                logger.error(f"Error processing GFF data: {str(e)}", exc_info=True)
                # Continue with empty annotations - we'll show the contigs at least
                gene_annotations = []
            
            # Sort sequences by length and get top N contigs
            if sequences:
                try:
                    # Sort by sequence length
                    sequences.sort(key=lambda x: x.get('length', 0), reverse=True)
                    # Limit to top 24 sequences for better visualization
                    sequences = sequences[:24]
                    logger.info(f"Sorted and limited to {len(sequences)} sequences")
                except Exception as e:
                    logger.error(f"Error processing sequences: {str(e)}", exc_info=True)
            
            # If no valid sequences, use sample data
            if not sequences:
                logger.info("Using sample sequences as fallback")
                # Create some fake sequence data
                sequences = [
                    {
                        'id': 'contig_1',
                        'header': '>contig_1',
                        'length': 5000000,
                        'sequence': 'A' * 100  # Just the first 100 bases for demo
                    },
                    {
                        'id': 'contig_2',
                        'header': '>contig_2',
                        'length': 4000000,
                        'sequence': 'G' * 100
                    }
                ]
            
            # Generate a unique ID for this visualization
            safe_mag_id = str(mag_id).replace('.', '_').replace(':', '_').replace(' ', '_') if mag_id else 'sample'
            unique_id = f"cgview-{safe_mag_id}-{int(time.time())}"
            
            # Create a CGView.js compatible data structure
            cgview_data = self._prepare_cgview_data(sequences, gene_annotations, safe_mag_id, mag_info)
            
            # Generate HTML with embedded CGView.js
            html_content = self._generate_cgview_html(unique_id, cgview_data, safe_mag_id)
            
            return html_content
            
        except Exception as e:
            logger.error(f"Error creating CGView plot: {str(e)}", exc_info=True)
            return f"""
                <html>
                <body>
                    <div style="color: red; padding: 20px; border: 1px solid #ccc;">
                        <h3>Error creating plot</h3>
                        <p>{str(e)}</p>
                        <pre>{traceback.format_exc()}</pre>
                    </div>
                </body>
                </html>
            """

    def _prepare_cgview_data(self, sequences, gene_annotations, mag_id, mag_info=None):
        """Prepare data for CGView.js format.
        
        Args:
            sequences: List of sequence dictionaries
            gene_annotations: List of annotation dictionaries
            mag_id: MAG identifier
            mag_info: Additional MAG metadata
            
        Returns:
            dict: Data structure for CGView.js
        """
        logger.info(f"Preparing CGView data for {mag_id} with {len(sequences)} sequences and {len(gene_annotations)} annotations")
        
        # Initialize a CGView data structure
        cgview_data = {
            "name": f"MAG {mag_id}",
            "id": mag_id,
            "width": 800,
            "height": 800,
            "settings": {
                "format": "circular",
                "backgroundColor": "white",
                "showShading": True,
                "arrowHeadLength": 0.3,
                "initialMapThicknessProportion": 0.1,
                "maxMapThicknessProportion": 0.5
            },
            "sequence": {
                "name": mag_id,
                "font": "SansSerif, plain, 14",
                "color": "black",
                "visible": True,
                "contigs": []
            },
            "legend": {
                "position": "top-right",
                "anchor": "auto",
                "defaultFont": "SansSerif, plain, 12",
                "defaultFontColor": "black",
                "textAlignment": "left",
                "backgroundColor": "white",
                "on": "canvas",
                "visible": True,
                "items": []
            },
            "backbone": {
                "thickness": 5,
                "color": "grey",
                "colorAlternate": "rgb(200,200,200)",
                "decoration": "arrow" if len(sequences) > 1 else "arc",
                "visible": True
            },
            "ruler": {
                "font": "sans-serif, plain, 10",
                "color": "black",
                "visible": True
            },
            "dividers": {
                "color": "grey",
                "thickness": 1,
                "spacing": 1,
                "visible": True
            },
            "annotation": {
                "font": "monospace, plain, 12",
                "color": "black",
                "onlyDrawFavorites": False,
                "labelPlacement": "default", 
                "visible": True
            },
            "tracks": [],
            "features": [],
            "captions": [
                {
                    "name": f"MAG: {mag_id}",
                    "position": "top-center",
                    "anchor": "middle-center",
                    "font": "SansSerif, plain, 16",
                    "fontColor": "#2c3e50",
                    "textAlignment": "center",
                    "backgroundColor": "white",
                    "on": "canvas",
                    "visible": True
                }
            ]
        }
        
        # Add metadata caption if available
        if mag_info:
            metadata_caption = {
                "name": f"Completeness: {float(mag_info.get('completeness', 0)) * 100:.1f}%, "
                       f"Contamination: {float(mag_info.get('contamination', 0)) * 100:.1f}%, "
                       f"Contigs: {len(sequences)}",
                "position": "bottom-center",
                "anchor": "middle-center",
                "font": "SansSerif, plain, 12",
                "fontColor": "#7f8c8d",
                "textAlignment": "center",
                "backgroundColor": "white",
                "on": "canvas",
                "visible": True
            }
            cgview_data["captions"].append(metadata_caption)
            
            # Add taxonomy caption if available
            if mag_info.get('taxonomy'):
                taxonomy = mag_info.get('taxonomy', 'Unknown')
                if isinstance(taxonomy, dict):
                    tax_parts = []
                    for level, name in taxonomy.items():
                        if name and name.lower() not in ['unknown', 'unclassified', 'none']:
                            tax_parts.append(f"{name}")
                    taxonomy = "; ".join(tax_parts) if tax_parts else "Unknown"
                
                cgview_data["captions"].append({
                    "name": f"Taxonomy: {taxonomy}",
                    "position": "top-center",
                    "anchor": "middle-center",
                    "font": "SansSerif, plain, 12",
                    "fontColor": "#7f8c8d",
                    "textAlignment": "center",
                    "backgroundColor": "white",
                    "on": "canvas",
                    "visible": True
                })
        
        # Add contigs
        total_length = 0
        for i, seq in enumerate(sequences):
            seq_id = seq.get('id', '').split()[0]
            length = seq.get('length', 0)
            total_length += length
            
            # Calculate GC content if not already present
            if 'gc_content' not in seq and 'sequence' in seq:
                sequence = seq.get('sequence', '')
                if sequence:
                    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
                    total = len(sequence)
                    seq['gc_content'] = (gc_count / total) * 100 if total > 0 else 0
            
            # Add contig to sequence
            contig = {
                "name": seq_id,
                "length": length,
                "orientation": "+",
                "color": self._get_contig_color(i),
                "visible": True
            }
            cgview_data["sequence"]["contigs"].append(contig)
        
        # Update sequence length
        cgview_data["sequence"]["length"] = total_length
        
        # Create legend items for different feature types
        legend_items = {
            "CDS": {
                "name": "CDS",
                "fontColor": "black",
                "swatchColor": "#1f77b4",
                "decoration": "arrow"
            },
            "rRNA": {
                "name": "rRNA",
                "fontColor": "black",
                "swatchColor": "#2ca02c",
                "decoration": "arrow"
            },
            "tRNA": {
                "name": "tRNA",
                "fontColor": "black",
                "swatchColor": "#d62728",
                "decoration": "arrow"
            },
            "ncRNA": {
                "name": "ncRNA",
                "fontColor": "black",
                "swatchColor": "#9467bd",
                "decoration": "arrow"
            },
            "repeat": {
                "name": "Repeat",
                "fontColor": "black",
                "swatchColor": "#8c564b",
                "decoration": "arrow"
            },
            "other": {
                "name": "Other",
                "fontColor": "black",
                "swatchColor": "#7f7f7f",
                "decoration": "arrow"
            }
        }
        
        # Add legend items
        for key, item in legend_items.items():
            cgview_data["legend"]["items"].append(item)
        
        # Create a "forward strand" track
        forward_track = {
            "name": "Forward Strand",
            "dataType": "feature",
            "dataMethod": "source",
            "dataKeys": "forward",
            "position": "outside",
            "separateFeaturesBy": "none",
            "thicknessRatio": 1,
            "visible": True
        }
        
        # Create a "reverse strand" track
        reverse_track = {
            "name": "Reverse Strand",
            "dataType": "feature",
            "dataMethod": "source",
            "dataKeys": "reverse",
            "position": "inside",
            "separateFeaturesBy": "none",
            "thicknessRatio": 1,
            "visible": True
        }
        
        # Add tracks
        cgview_data["tracks"].append(forward_track)
        cgview_data["tracks"].append(reverse_track)
        
        # Create a GC content track if we have sequence data
        if any('sequence' in seq for seq in sequences):
            # Create a GC content plot
            gc_plot = {
                "name": "GC Content",
                "legend": "GC Content",
                "baseline": 50,
                "axisMax": 70,
                "axisMin": 30,
                "positions": [],
                "scores": []
            }
            
            # Add plot to a new track
            gc_track = {
                "name": "GC Content",
                "dataType": "plot",
                "dataMethod": "source",
                "dataKeys": "gc_content",
                "position": "outside",
                "thicknessRatio": 0.5,
                "visible": True
            }
            
            cgview_data["tracks"].append(gc_track)
        
        # Process gene annotations and add features
        if gene_annotations:
            feature_count = 0
            for annotation in gene_annotations:
                try:
                    # Some basic validation
                    if not isinstance(annotation, dict):
                        continue
                    
                    contig = annotation.get('contig')
                    if not contig:
                        continue
                    
                    start = annotation.get('start')
                    end = annotation.get('end')
                    if not start or not end:
                        continue
                    
                    # Determine strand
                    strand = annotation.get('strand')
                    strand_string = '+' if strand in ['+', 1, '1'] else '-'
                    source = "forward" if strand_string == '+' else "reverse"
                    
                    # Determine feature type and legend
                    feature_type = annotation.get('type', 'unknown').lower()
                    feature_type_orig = annotation.get('original_type', feature_type).lower()
                    
                    # Map feature type to legend
                    legend = "CDS"
                    if 'rrna' in feature_type_orig:
                        legend = "rRNA"
                    elif 'trna' in feature_type_orig:
                        legend = "tRNA"
                    elif 'ncrna' in feature_type_orig or 'nc_rna' in feature_type_orig:
                        legend = "ncRNA"
                    elif 'repeat' in feature_type_orig:
                        legend = "repeat"
                    elif feature_type_orig not in ['cds', 'gene']:
                        legend = "other"
                    
                    # Get feature name
                    feature_name = annotation.get('name', '') or annotation.get('gene', '') or annotation.get('locus_tag', '')
                    if not feature_name:
                        feature_name = f"{feature_type}_{start}"
                    
                    # Create feature object
                    feature = {
                        "name": feature_name,
                        "type": feature_type,
                        "legend": legend,
                        "source": source,
                        "contig": contig,
                        "start": start,
                        "stop": end,
                        "strand": strand_string,
                        "visible": True
                    }
                    
                    # Add function/product if available
                    product = annotation.get('product', '') or annotation.get('function', '')
                    if product:
                        feature["meta"] = {
                            "product": product
                        }
                    
                    # Add feature
                    cgview_data["features"].append(feature)
                    feature_count += 1
                    
                    # Log progress occasionally
                    if feature_count % 500 == 0:
                        logger.info(f"Processed {feature_count} features")
                    
                except Exception as e:
                    logger.error(f"Error processing annotation: {str(e)}")
                    continue
        
        logger.info(f"Prepared CGView data with {len(cgview_data['features'])} features and {len(cgview_data['tracks'])} tracks")
        return cgview_data

    def _generate_cgview_html(self, unique_id, cgview_data, mag_id):
        """Generate HTML content with embedded CGView.js.
        
        Args:
            unique_id: Unique identifier for the plot
            cgview_data: CGView.js compatible data structure
            mag_id: MAG identifier
            
        Returns:
            str: HTML content
        """
        # Escape the data for embedding in JavaScript
        import json
        cgview_data_json = json.dumps(cgview_data)
        
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <meta charset="utf-8">
            <meta name="viewport" content="width=device-width, initial-scale=1">
            <title>MAG {mag_id} Visualization</title>
            <style>
                body {{
                    margin: 0;
                    padding: 0;
                    font-family: Arial, sans-serif;
                }}
                .container {{
                    width: 100%;
                    height: 100vh;
                    display: flex;
                    flex-direction: column;
                }}
                .header {{
                    padding: 10px;
                    text-align: center;
                    background-color: #f8f9fa;
                    border-bottom: 1px solid #ddd;
                }}
                .plot-container {{
                    flex: 1;
                    position: relative;
                    overflow: hidden;
                }}
                #cgview-container {{
                    width: 100%;
                    height: 100%;
                }}
                .loading {{
                    position: absolute;
                    top: 50%;
                    left: 50%;
                    transform: translate(-50%, -50%);
                    text-align: center;
                    color: #666;
                }}
                .debug-info {{
                    padding: 10px;
                    font-size: 12px;
                    color: #666;
                    background-color: #f8f9fa;
                    border-top: 1px solid #ddd;
                    max-height: 150px;
                    overflow: auto;
                }}
            </style>
            <!-- Load D3.js and CGView.js -->
            <script src="https://d3js.org/d3.v7.min.js"></script>
            <script src="https://cdn.jsdelivr.net/npm/cgview/dist/cgview.min.js"></script>
            <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/cgview/dist/cgview.css">
        </head>
        <body>
            <div class="container">
                <div class="header">
                    <h2>MAG: {mag_id}</h2>
                </div>
                <div class="plot-container">
                    <div id="loading-{unique_id}" class="loading">
                        <p>Loading genome visualization...</p>
                    </div>
                    <div id="{unique_id}" class="cgview-container"></div>
                </div>
                <div id="debug-{unique_id}" class="debug-info">
                    <p>Initializing CGView.js...</p>
                </div>
            </div>
            <script>
                // Initialize CGView when the page is loaded
                document.addEventListener('DOMContentLoaded', function() {{
                    const debugEl = document.getElementById('debug-{unique_id}');
                    const loadingEl = document.getElementById('loading-{unique_id}');
                    
                    function logDebug(message) {{
                        console.log(message);
                        debugEl.innerHTML += '<p>' + message + '</p>';
                    }}
                    
                    try {{
                        logDebug('Loading CGView.js for MAG {mag_id}...');
                        
                        // Set up the viewer
                        const cgvContainer = document.getElementById('{unique_id}');
                        const cgvData = {cgview_data_json};
                        
                        logDebug('CGView container: ' + cgvContainer);
                        logDebug('Creating CGView instance...');
                        
                        // Create the CGView instance
                        const cgv = new CGV.Viewer(cgvContainer, cgvData);
                        
                        // Listen for loading errors
                        cgv.on('error', function(error) {{
                            logDebug('Error loading CGView: ' + error.message);
                            loadingEl.innerHTML = '<p style="color: red;">Error loading visualization: ' + error.message + '</p>';
                        }});
                        
                        // When fully loaded, hide the loading indicator
                        cgv.on('ready', function() {{
                            logDebug('CGView loaded successfully!');
                            loadingEl.style.display = 'none';
                        }});
                        
                        // Additional event listeners for debugging
                        cgv.on('features-add', function(features) {{
                            logDebug('Features added: ' + features.length);
                        }});
                        
                        cgv.on('tracks-add', function(tracks) {{
                            logDebug('Tracks added: ' + tracks.length);
                        }});
                        
                        // If loading takes too long, provide additional information
                        setTimeout(function() {{
                            if (loadingEl.style.display !== 'none') {{
                                logDebug('Still loading after 5 seconds. You may need to check browser console for errors.');
                            }}
                        }}, 5000);
                        
                    }} catch (error) {{
                        logDebug('Error initializing CGView: ' + error.message);
                        loadingEl.innerHTML = '<p style="color: red;">Error initializing visualization: ' + error.message + '</p>';
                        console.error(error);
                    }}
                }});
            </script>
        </body>
        </html>
        """
        
        return html

    def _get_contig_color(self, index):
        """Get a color for a contig based on its index."""
        colors = [
            "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
            "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
            "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"
        ]
        return colors[index % len(colors)]

  