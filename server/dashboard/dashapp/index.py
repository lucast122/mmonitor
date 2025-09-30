import base64
import io
import json
import re
from typing import Any, List, Union
import time
import traceback
import functools
import logging
from functools import lru_cache
import hashlib
import cProfile
import pstats
import base64
import io
import json
import re
from typing import Any, List, Union
import time
import sys
import traceback
import functools
import logging
from functools import lru_cache
import pstats
from functools import wraps
import asyncio
from datetime import datetime
import numpy as np
import dash
from dash import html, dcc, dash_table
import dash_mantine_components as dmc
from dash.dependencies import Input, Output, State, ALL
from dash.exceptions import PreventUpdate
from django_plotly_dash import DjangoDash
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import pandas as pd
from dash_iconify import DashIconify
from scipy.ndimage import gaussian_filter
from scipy.stats import zscore
from statsmodels.nonparametric.smoothers_lowess import lowess
import random
from functools import wraps
from concurrent.futures import ThreadPoolExecutor
import asyncio
from datetime import datetime, timedelta
import numpy as np
import dash
from dash import html, dcc, dash_table
import dash_mantine_components as dmc
# import dash_bio as dashbio  # Temporarily commented out
from dash.dependencies import Input, Output, State, ALL
from dash.exceptions import PreventUpdate
from django_plotly_dash import DjangoDash
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots
import pandas as pd
from dash_iconify import DashIconify
from scipy.ndimage import gaussian_filter
from scipy.stats import zscore
from statsmodels.nonparametric.smoothers_lowess import lowess
import random
from dash.dependencies import Input, Output, State 
from dash.exceptions import PreventUpdate
from dash import no_update
import logging
import traceback

# Try to import pyCirclize for MAG visualization
try:
    from pycirclize import Circos
    PYCIRCLIZE_AVAILABLE = True
except ImportError:
    PYCIRCLIZE_AVAILABLE = False

from users.models import NanoporeRecord, Metadata, User, SequencingStatistics, MAG
from . import taxonomy, correlations, qc, diversity, horizon, mags

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def get_current_user(request=None):
    """Get the current user from the Django request context."""
    if request and hasattr(request, 'user'):
        return request.user
    return None

def benchmark(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        duration = time.time() - start
        logger.info(f"Function {func.__name__} took {duration:.2f} seconds to execute")
        return result
    return wrapper

def profile(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        profiler = cProfile.Profile()
        try:
            return profiler.runcall(func, *args, **kwargs)
        finally:
            stats = pstats.Stats(profiler)
            stats.sort_stats('cumulative')
            logger.info(f"\n{'='*20} Profiling for {func.__name__} {'='*20}")
            stats.print_stats(20)  # Print top 20 time-consuming operations
    return wrapper



class Index:
    """
    The Index class is the main class of the application.
    It contains the main layout and callbacks.
    """

    def __init__(self, user_id):
        print("Initializing Index app")
        self.user_id = user_id
        print(sys.executable)  # Add this to a Django view temporarily to check
        # Initialize the app with suppress_callback_exceptions=True to handle cross-app callbacks
        self.app = DjangoDash('dashboard_index', add_bootstrap_links=True, suppress_callback_exceptions=True)
        
        # Initialize all sub-apps
        print("Initializing Correlations app")
        self.correlations_app = correlations.Correlations(user_id)
        
        # Initialize other apps...
        
        self.taxonomy_app = taxonomy.Taxonomy(user_id)
        self.diversity_app = diversity.Diversity(user_id)
        self.correlations_app = correlations.Correlations(user_id)
        self.horizon_app = horizon.Horizon(user_id)
        self.qc_app = qc.QC(user_id)
        self.mags_app = mags.MAGs(user_id)


        
        self._apps = {
            '/dashapp/taxonomy': {
                'name': 'Taxonomy',
                'app': self.taxonomy_app.app,
                'instance': self.taxonomy_app
            },
            '/dashapp/horizon': {
                'name': 'Horizon',
                'app': self.horizon_app.app,
                'instance': self.horizon_app
            },
            '/dashapp/diversity': {
                'name': 'Diversity',
                'app': self.diversity_app.app,
                'instance': self.diversity_app
            },
            '/dashapp/correlations': {
                'name': 'Correlations',
                'app': self.correlations_app.app,
                'instance': self.correlations_app
            },
            '/dashapp/qc': {
                'name': 'QC',
                'app': self.qc_app.app,
                'instance': self.qc_app
            },
            '/dashapp/mags': {
                'name': 'MAGs Viewer',
                'app': self.mags_app.app,
                'instance': self.mags_app
            }
        }
        
        # Ensure all apps have suppress_callback_exceptions=True
        for app_path, app_info in self._apps.items():
            if hasattr(app_info['app'], '_suppress_callback_exceptions'):
                app_info['app']._suppress_callback_exceptions = True

        self._init_data()
        
        print("Initializing layout")
        self._init_layout()
        print("Initializing callbacks")
        self._init_callbacks()
        print("Index app initialization complete")

    def _init_data(self):
        """Initialize data with optimized database queries"""
        if hasattr(self, 'df'):
            return
        
        logger.info("Starting data initialization...")
        
        try:
            # Use select_related to optimize database queries
            records = (NanoporeRecord.objects
                       .filter(user_id=self.user_id)
                    #    .select_related()
                       .values())
            
            # Create DataFrame with optimized dtypes
            self.df = pd.DataFrame.from_records(records)
            self.df_sorted = self.df.sort_values(by=['project_id', 'subproject', 'sample_id'])
            
            # Initialize cache with size limit
            self._cache = {}
            self._cache_timeout = 3600  # 1 hour cache timeout
            self._cache_size_limit = 100  # Maximum number of items in cache
            
            if len(self.df) > 0:
                # Optimize memory usage
                self.df = self._optimize_dataframe(self.df)
                
                # Pre-calculate common aggregations
                self._precalculate_aggregations()
            
            logger.info(f"Data initialization complete. DataFrame shape: {self.df.shape}")
            
        except Exception as e:
            logger.error(f"Error initializing data: {str(e)}")
            logger.error(traceback.format_exc())
            self.df = pd.DataFrame()  # Initialize empty DataFrame on error
            self.df_sorted = pd.DataFrame()
            self._cache = {}
            self._cache_timeout = 3600  # 1 hour cache timeout
            self._cache_size_limit = 100  # Maximum number of items in cache
            logger.info("Initialized empty dataframe due to error")
        
        # Initialize cache with size limit
        self._cache = {}
        self._cache_timeout = 3600  # 1 hour cache timeout
        self._cache_size_limit = 100  # Maximum number of items in cache
        
        
        
        if len(self.df) > 0:
            # Optimize memory usage
            self.df = self._optimize_dataframe(self.df)
            
            # Pre-calculate common aggregations
            self._precalculate_aggregations()
        

        logger.info(f"Data initialization complete. DataFrame shape: {self.df.shape}")

    def _get_icon_for_app(self, app_name):
        icon_map = {
            'Taxonomy': 'fluent-emoji-high-contrast:microbe',
            'Horizon': 'mdi:chart-timeline-variant',
            'Diversity': 'tabler:chart-scatter-3d',
            'Correlations': 'carbon:qq-plot',
            'QC': 'mdi:check-circle-outline',
            'MAGs Viewer': 'mdi:dna'
        }
        return icon_map.get(app_name, 'mdi:application')

    def _init_layout(self):
        
        hidden_components = html.Div([
            # Data components 
            dcc.Store(id='session-id', storage_type='session'),
            dcc.Location(id='redirect', refresh=True),
            dcc.Download(id="download-svg"),
            dcc.Download(id="download-csv"),
            dcc.Download(id="download-legend-svg"),
            dcc.Download(id="horizon-plot-download"),
            dcc.Download(id="download-pdf"),
            
            # Empty container for horizon plot
            html.Div(id="horizon-plot-container", style={'display': 'none'}),
            
            # Use div instead of iframe to avoid security issues
            html.Div(id="horizon-plot-content", style={'display': 'none'}),
            
            # Store components
            dcc.Store(id="initialization-store"),
            dcc.Store(id="data-store"),
            dcc.Store(id="ngcircos-status", data="unchecked"),
            
            # Add a hidden div for the format template fix
            html.Div(id='fix-format-template-issue', style={'display': 'none'}),
            
            # Error container
            html.Div(id="global-error-container", style={'display': 'none'}),
            
            # CSS for animations and circos plots - using dcc.Markdown instead of html.Style
            html.Div([
                dcc.Markdown('''
                <style>
                @keyframes spin {
                    0% { transform: rotate(0deg); }
                    100% { transform: rotate(360deg); }
                }
                
                .circos-container {
                    width: 100%;
                    height: 600px;
                    position: relative;
                    overflow: visible !important;
                }
                
                .circos-plot {
                    width: 100%;
                    height: 100%;
                    text-align: center;
                }
                
                .circos-svg {
                    max-width: 100%;
                    max-height: 100%;
                }
                
                /* Fix for NGCircos and D3 SVG sizing */
                #circos-plot svg {
                    width: 100% !important;
                    height: 100% !important;
                    overflow: visible !important;
                }
                
                .dash-bio-circos {
                    width: 100% !important;
                    height: 600px !important;
                }
                
                .ngcircos-container {
                    width: 100%;
                    height: 600px;
                    position: relative;
                    border: 1px solid #ddd;
                    background: #fcfcfc;
                }
                
                .ngcircos-zoom-controls {
                    position: absolute;
                    top: 10px;
                    left: 10px;
                    z-index: 10;
                    background: rgba(255,255,255,0.8);
                    padding: 5px;
                    border-radius: 3px;
                    border: 1px solid #ddd;
                }
                
                .ngcircos-zoom-btn {
                    margin: 2px;
                    padding: 5px 10px;
                    border: 1px solid #ddd;
                    background: white;
                    cursor: pointer;
                    font-weight: bold;
                }
                
                .ngcircos-status {
                    position: absolute;
                    top: 50%;
                    left: 50%;
                    transform: translate(-50%, -50%);
                    background: rgba(255,255,255,0.9);
                    padding: 20px;
                    border-radius: 5px;
                    border: 1px solid #ddd;
                    text-align: center;
                    z-index: 5;
                }
                
                .gene-legend {
                    margin-top: 10px;
                    padding: 5px;
                    display: flex;
                    flex-wrap: wrap;
                    justify-content: center;
                }
                
                .gene-legend-item {
                    display: flex;
                    align-items: center;
                    margin: 5px 10px;
                }
                
                .gene-legend-color {
                    width: 20px;
                    height: 20px;
                    margin-right: 5px;
                    border-radius: 3px;
                }
                
                .gene-legend-label {
                    font-size: 12px;
                }
                
                .loader {
                    border: 16px solid #f3f3f3;
                    border-top: 16px solid #3498db;
                    border-radius: 50%;
                    width: 120px;
                    height: 120px;
                    animation: spin 2s linear infinite;
                }
                </style>
                ''', dangerously_allow_html=True)
            ], style={'display': 'block'}),
            
            # Hidden divs to suppress "ID not found in layout" warnings for other apps' components
            html.Div([
                dcc.Graph(id="gene-content-plot", style={'display': 'none'}),
                dcc.Dropdown(id="sequence-selector", style={'display': 'none'}),
                dcc.Graph(id="heatmap-dendrogram", style={'display': 'none'}),
                html.Div(id="notification-output", style={'display': 'none'}),
                html.Div(id="ngcircos-status-display", style={'display': 'none'}),
                dash_table.DataTable(
                    id="mag-table",
                    data=[],
                    columns=[],
                    selected_rows=[],
                    style_table={'display': 'none'}
                ),
                dcc.Graph(id="average-read-length-plot", style={'display': 'none'}),
                dcc.Graph(id="read-length-dist", style={'display': 'none'}),
                dcc.Graph(id="gc-dist", style={'display': 'none'}),
                dcc.Graph(id="base-composition", style={'display': 'none'}),
                html.Iframe(id="horizon-plot-iframe", style={'display': 'none'}),
                html.Div(id="mag-genome-warning", style={'display': 'none'}),
                html.Div(id="mag-info-section", style={'display': 'none'}),
                html.Div(id="mags-debug-info", style={'display': 'none'}),
                html.Div(id="circos-plot", style={'display': 'none'}),
                html.Div(id="gene-search-results", style={'display': 'none'}),
                html.Div(id="sequence-viewer-container", style={'display': 'none'}),
                html.Div(id="analysis-output", style={'display': 'none'}),
                dcc.DatePickerRange(id="date-range-picker", style={'display': 'none'}),
                html.Div(id="copy-alert", style={'display': 'none'}),
                dcc.Graph(id="pcoa_plot_container", style={'display': 'none'}),
                dcc.Graph(id="rarefaction_plot", style={'display': 'none'}),
                dcc.Download(id="download-diversity-csv"),
                dcc.Download(id="download-corr-csv"),
                dcc.Download(id="download-svg-diversity"),
                dcc.Graph(id="graph1", style={'display': 'none'}),
                dcc.Download(id="download-plot")
            ], id='suppress-warnings', style={'display': 'none'})
        ], id='hidden-components', style={'display': 'none'})

        
        
        location = dcc.Location(id='url', refresh=True)

        nav_links = [
            dmc.NavLink(
                label=values['name'],
                style={'font-size': '1.2rem'},
                leftSection=DashIconify(icon=self._get_icon_for_app(values['name']), width=35),
                href=url,
            ) for url, values in self._apps.items()
        ]

        navbar = dmc.AppShellNavbar(
            p="xs",
            children=[
                dmc.ScrollArea(
                    type="hover",
                    children=[
                        dmc.Stack([
                            dmc.Text(id="app-subtitle", size="xs", c="dimmed"),
                            dmc.Divider(mb=20),
                            *nav_links
                        ])
                    ]
                )
            ],
        )

        header = dmc.AppShellHeader(
            p="md",
            children=[
                dmc.Group(
                    [
                        dmc.Text(id="header-title", size="xl", fw=700),
                    ],
                    gap="xl",
                )
            ],
        )

        initial_content = html.Div("Please select an app from the menu.", id='initial-content')
        page_content = html.Div(id='page-content', children=[initial_content], style={'width': '100%', 'padding': '0'})

        
        layout = dmc.AppShell(
            [
                hidden_components,
                location,
                header,
                navbar,
                dmc.AppShellMain(
                    children=[
                        dmc.Container(
                            id="page-content",
                            fluid=True,
                            style={
                                'padding': '0',
                                'height': 'calc(100vh - 20px)',
                                'overflow': 'auto',
                                'minHeight': '0',
                                'maxHeight': '100vh'
                            }
                        )
                    ]
                ),
                dmc.Text(id="header-title", style={"display": "none"})
            ],
            style={
                'margin': '0',
                'padding': '0',
                'width': '100%',
                'maxWidth': '100%',
                'height': '100vh',
                'overflow': 'hidden',
                'display': 'flex',
                'flexDirection': 'column',
                'minHeight': '0'  # Added minHeight
            }
        )

        container = dmc.MantineProvider(
            theme={
                "colorScheme": "light",
                "primaryColor": "indigo",
                "fontFamily": "Arial, sans-serif",
                "spacing": {"xs": 4, "sm": 8, "md": 12, "lg": 16, "xl": 20},
            },
            children=[layout]
        )

        self.app.layout = container
        print("Initialized Index layout")

    def _get_app_icon(self, app_name: str) -> str:
        """Get the icon name for a given app."""
        # For custom image apps, use a special prefix to indicate custom icons
        custom_icons = {
            'Taxonomy': 'custom:microorganism.png',
            'MAGs': 'custom:circ-chromosome.png'
        }
        
        if app_name in custom_icons:
            return custom_icons[app_name]
            
        # Standard iconify icons
        icon_map = {
            'Horizon': 'mdi:chart-timeline',
            'Diversity': 'tabler:chart-scatter-3d',
            'Correlations': 'carbon:qq-plot',
            'QC': 'material-symbols:fact-check-outline'
        }
        return icon_map.get(app_name, 'mdi:application')

    def _get_cache_key(self, prefix, *args, **kwargs):
        """Generate a cache key from arguments"""
        key_dict = {
            'args': args,
            'kwargs': kwargs,
            'user_id': self.user_id
        }
        key_str = json.dumps(key_dict, sort_keys=True)
        return f"{prefix}:{hashlib.md5(key_str.encode()).hexdigest()}"

    @lru_cache(maxsize=128)
    def _get_cached_data(self):
        """Get cached data with performance optimizations"""
        start = time.time()
        try:
            records = NanoporeRecord.objects.filter(user_id=self.user_id).values()
            df = pd.DataFrame.from_records(records)
            duration = time.time() - start
            logger.info(f"Data fetching and processing took {duration:.2f} seconds")
            return df
        except Exception as e:
            logger.error(f"Error fetching data: {str(e)}")
            return pd.DataFrame()

    def _optimize_dataframe(self, df):
        """Optimize DataFrame memory usage"""
        start = time.time()
        
        # Convert date columns to datetime
        date_cols = [col for col in df.columns if 'date' in col.lower()]
        for col in date_cols:
            df[col] = pd.to_datetime(df[col], errors='ignore')
        
        # Optimize numeric columns
        for col in df.select_dtypes(include=['float64']).columns:
            df[col] = pd.to_numeric(df[col], downcast='float')
        
        for col in df.select_dtypes(include=['int64']).columns:
            df[col] = pd.to_numeric(df[col], downcast='integer')
        
        # Categorize string columns with low cardinality
        for col in df.select_dtypes(include=['object']).columns:
            unique_ratio = df[col].nunique() / len(df)
            if unique_ratio < 0.5:  # If less than 50% unique values
                df[col] = df[col].astype('category')
        
        duration = time.time() - start
        memory_usage = df.memory_usage(deep=True).sum() / 1024**2
        logger.info(f"DataFrame optimization took {duration:.2f} seconds. Memory usage: {memory_usage:.2f} MB")
        
        return df

    

    def _init_callbacks(self):
        """Initialize all callbacks for the index page."""
        
        logger = logging.getLogger(__name__)
        logger.info("Initializing callbacks...")
        
        # Add a global error handler to display in the debug console
        @self.app.callback(
            Output("global-error-container", "children"),
            [Input("url", "pathname")],
            prevent_initial_call=True
        )
        def handle_global_errors(pathname):
            try:
                # This function is just for error handling setup
                return []
            except Exception as e:
                logger.error(f"Global error: {str(e)}")
                return html.Div(f"Error: {str(e)}")
        
        # Register the clientside callback to show the loading indicator on navigation
        

        @self.app.callback(
            Output('mag-info-panel-content', 'children'),
            [Input('mag-table-selected-rows', 'data')],
            [State('mag-table-data', 'children')],
            prevent_initial_call=True, allow_duplicate=True
        )
        def update_mag_info_panel(selected_rows, table_data_json):
            """Update the MAG info panel with details of the selected MAG."""
            if not selected_rows or len(selected_rows) == 0:
                return dmc.Alert(
                    "Please select a MAG from the table to view detailed information.",
                    color="blue",
                    variant="light",
                    title="No MAG Selected"
                )
            
            try:
                # Parse the JSON table data
                table_data = json.loads(table_data_json)
                selected_row = table_data[selected_rows[0]]
                
                # Format percentages and values correctly
                completeness = selected_row.get('completeness', 0)
                contamination = selected_row.get('contamination', 0)
                gc_content = selected_row.get('gc_content', 0)
                
                # Convert to percentages if needed
                completeness_display = f"{completeness*100:.1f}%" if completeness <= 1 else f"{completeness:.1f}%"
                contamination_display = f"{contamination*100:.1f}%" if contamination <= 1 else f"{contamination:.1f}%"
                gc_content_display = f"{gc_content*100:.1f}%" if gc_content <= 1 else f"{gc_content:.1f}%"
                
                # Format size
                size_bp = selected_row.get('size', 0)
                if size_bp == 0 and 'size_mb' in selected_row:
                    size_bp = selected_row.get('size_mb', 0) * 1000000
                size_display = f"{size_bp:,.0f} bp"
                
                # Format contigs, N50, longest contig
                num_contigs = selected_row.get('num_contigs', 'N/A')
                n50 = selected_row.get('n50', 0)
                longest_contig = selected_row.get('longest_contig', 0)
                
                # Create a list of MAG properties to display
                properties = [
                    ("Name", selected_row.get('name', 'N/A')),
                    ("Bin ID", selected_row.get('bin_id', 'N/A')),
                    ("Taxonomy", selected_row.get('taxonomy', 'N/A')),
                    ("Completeness", completeness_display),
                    ("Contamination", contamination_display),
                    ("Size", size_display),
                    ("GC Content", gc_content_display),
                    ("Number of Contigs", f"{num_contigs:,}" if isinstance(num_contigs, (int, float)) else num_contigs),
                    ("N50", f"{n50:,.0f} bp" if n50 else "N/A"),
                    ("Longest Contig", f"{longest_contig:,.0f} bp" if longest_contig else "N/A")
                ]
                
                # Create a list of property rows
                property_rows = [
                    dmc.Group([
                        dmc.Text(prop[0], fw=500, size="sm", style={"width": "150px"}),
                        dmc.Text(prop[1], size="sm")
                    ], justify="apart", style={"marginBottom": "8px"})
                    for prop in properties
                ]
                
                return dmc.Stack(property_rows, gap="xs")
                
            except Exception as e:
                logger.error(f"Error updating MAG info panel: {str(e)}", exc_info=True)
                return dmc.Alert(
                    "Error loading MAG information. Please try again.",
                    color="red",
                    variant="light",
                    title="Error"
                )

        
        @self.app.callback(
            Output('heatmap-dendrogram', 'figure'),
            [Input('apply-correlation-settings', 'n_clicks')],
            [State('correlation-method', 'value'),
             State('pvalue-cutoff', 'value'),
             State('filter-significant', 'checked')],
            prevent_initial_call=True
        )
        def update_heatmap(n_clicks, correlation_method, p_value_cutoff, filter_significant):
            if n_clicks is None:
                raise PreventUpdate
                
            # Log the inputs for debugging
            logger.info(f"Updating heatmap with method={correlation_method}, p_value={p_value_cutoff}, filter={filter_significant}")
                
            # Update correlation method if changed
            if correlation_method != self.correlations_app.corr_method:
                logger.info(f"Correlation method changed from {self.correlations_app.corr_method} to {correlation_method}")
                self.correlations_app.corr_method = correlation_method
                # Recalculate correlations with new method
                self.correlations_app.correlation_scores = self.correlations_app.compute_correlations_for_taxonomies(
                    self.correlations_app.nanopore_df, 
                    self.correlations_app.meta_df, 
                    correlation_method,  # Use the new method directly
                    'taxonomy'
                )
            
            # Update the heatmap with the current settings
            return self.correlations_app.create_heatmap(p_value_cutoff=p_value_cutoff, filter_significant=filter_significant)
        
        # Add callback for initial load to set default values properly
        @self.app.callback(
            Output('heatmap-dendrogram', 'figure', allow_duplicate=True),
            [Input('correlation-method', 'value')],
            prevent_initial_call=True
        )
        def initial_heatmap(correlation_method):
            # Don't recompute on first load
            try:
                return self.correlations_app.create_heatmap(p_value_cutoff=0.05, filter_significant=False)
            except Exception as e:
                logger.error(f"Error in initial heatmap: {str(e)}")
                logger.error(traceback.format_exc())
                empty_fig = go.Figure()
                empty_fig.add_annotation(
                    text=f"Error: {str(e)}",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5, showarrow=False,
                    font=dict(size=20)
                )
                return empty_fig
        
        @self.app.callback(
            Output('notification-output', 'children'),
            [Input('apply-correlation-settings', 'n_clicks')],
            [State('correlation-method', 'value'),
             State('pvalue-cutoff', 'value'),
             State('filter-significant', 'checked')],
            prevent_initial_call=True
        )
        def show_applied_settings(n_clicks, correlation_method, p_value_cutoff, filter_significant):
            if n_clicks is None:
                raise PreventUpdate
                
            try:
                logger.info(f"Applying settings: method={correlation_method}, p_value={p_value_cutoff}, filter={filter_significant}")
                    
                notification_text = f"Applied settings: {correlation_method.capitalize()} correlation"
                if filter_significant:
                    notification_text += f", showing correlations with p < {p_value_cutoff}"
                else:
                    notification_text += ", showing all correlations"
                
                return dmc.Alert(
                    id="settings-notification",
                    title="Settings Applied",
                    children=notification_text,
                    color="green",
                    withCloseButton=True,
                    variant="filled",
                    style={"marginBottom": "15px"}
                )
            except Exception as e:
                logger.error(f"Error showing applied settings: {str(e)}")
                logger.error(traceback.format_exc())
                return dmc.Alert(
                    id="settings-notification",
                    title="Error",
                    children=f"Error applying settings: {str(e)}",
                    color="red",
                    withCloseButton=True,
                    variant="filled",
                    style={"marginBottom": "15px"}
                )

        
        
        
        @self.app.callback(
            Output('sequence-selector', 'options', allow_duplicate=True),
            [Input('mag-table', 'selected_rows')],
            [State('mag-table', 'data')],
            prevent_initial_call=True
        )
        def update_sequence_selector(selected_rows, table_data):
            """Update the sequence selector dropdown with sequences from the selected MAG."""
            from users.models import MAG
            import logging
            import traceback
            
            logger = logging.getLogger(__name__)
            
            if not selected_rows or not table_data or len(selected_rows) == 0:
                return []
            
            try:
                selected_row = table_data[selected_rows[0]]
                bin_id = selected_row.get('bin_id')
                if not bin_id:
                    logger.warning(f"No bin_id found in selected row: {selected_row}")
                    return []
                
                logger.info(f"Fetching sequence data for MAG: {bin_id}")
                
                try:
                    # Get MAG from database
                    mag = MAG.objects.get(name=bin_id, user_id=self.user_id)
                    if not mag.fasta_file:
                        logger.warning(f"No FASTA file found for MAG {bin_id}")
                        return []
                    
                    # Parse FASTA file to get sequences
                    sequences = []
                    current_header = None
                    current_seq = ""
                    
                    for line in mag.fasta_file.split('\n'):
                        line = line.strip()
                        if not line:
                            continue
                        
                        if line.startswith('>'):
                            if current_header and current_seq:
                                # Get ID from header (first word after >)
                                seq_id = current_header.split()[0]
                                sequences.append({
                                    'id': seq_id,
                                    'header': f">{current_header}",
                                    'length': len(current_seq),
                                    'sequence': current_seq
                                })
                            current_header = line[1:]  # Remove '>' from start
                            current_seq = ""
                        else:
                            current_seq += line
                    
                    # Add the last sequence
                    if current_header and current_seq:
                        seq_id = current_header.split()[0]
                        sequences.append({
                            'id': seq_id,
                            'header': f">{current_header}",
                            'length': len(current_seq),
                            'sequence': current_seq
                        })
                    
                    # Sort by length (largest first)
                    sequences.sort(key=lambda x: x['length'], reverse=True)
                    
                    # Limit to top 50 sequences to avoid performance issues
                    sequences = sequences[:50]
                    
                    # Create dropdown options
                    options = [{'label': f"{s['id']} ({s['length']:,} bp)", 'value': s['id']} for s in sequences]
                    
                    logger.info(f"Found {len(options)} sequences for MAG {bin_id}")
                    return options
                
                except MAG.DoesNotExist:
                    logger.error(f"MAG {bin_id} not found in database")
                    return []
                except Exception as e:
                    logger.error(f"Error fetching MAG data: {str(e)}")
                    logger.error(traceback.format_exc())
                    return []
                
            except Exception as e:
                logger.error(f"Error updating sequence selector: {str(e)}")
                logger.error(traceback.format_exc())
                return []
        
        # Add callback to update gene content plot
        @self.app.callback(
            Output('gene-content-plot', 'figure'),
            [Input('mag-table', 'selected_rows'),
             Input('mag-table', 'data')]
        )
        def update_gene_content_plot(selected_rows, table_data):
            """Update gene content plot showing gene counts from MAGs."""
            try:
                # Default empty figure
                empty_fig = go.Figure()
                empty_fig.update_layout(
                    title='Gene Content Distribution (Top 20 MAGs)',
                    xaxis_title='MAG',
                    yaxis_title='Gene Count',
                    showlegend=False,
                    height=400
                )
                
                # If no data, return empty figure
                if not table_data or len(table_data) == 0:
                    return empty_fig
                
                # Get the mag data from the mags app
                mag_data = self.mags_app.mag_data
                if mag_data is None or mag_data.empty:
                    return empty_fig
                
                # Add gene counts based on CDS count since gene_count field doesn't exist
                df = mag_data.copy()
                df['gene_count'] = df['cds_count']
                
                # Get highlight mag if available
                highlight_mag = None
                if selected_rows and len(selected_rows) > 0 and selected_rows[0] < len(table_data):
                    highlight_mag = table_data[selected_rows[0]].get('name')
                
                # Sort by gene count and get top 20
                df = df.sort_values('gene_count', ascending=False).head(20)
                
                # Create figure
                fig = go.Figure()
                
                # Add bars for each MAG
                for idx, row in df.iterrows():
                    mag_name = row['name']
                    is_highlighted = highlight_mag is not None and mag_name == highlight_mag
                    
                    fig.add_trace(go.Bar(
                        x=[mag_name],
                        y=[row['gene_count']],
                        name=mag_name,
                        showlegend=False,
                        marker_color='#3498db' if not is_highlighted else '#e74c3c',
                        opacity=1.0 if is_highlighted else 0.7,
                        hovertemplate='<b>%{x}</b><br>' +
                                     'Gene Count: %{y}<br>' +
                                     '<extra></extra>'
                    ))
                
                # Update layout
                fig.update_layout(
                    title='Gene Content Distribution (Top 20 MAGs)',
                    xaxis_title='MAG',
                    yaxis_title='Gene Count',
                    xaxis_tickangle=45,
                    showlegend=False,
                    hovermode='closest',
                    height=400
                )
                
                return fig
                
            except Exception as e:
                logger.error(f"Error updating gene content plot: {str(e)}", exc_info=True)
                # Return empty figure if there's an error
                empty_fig = go.Figure()
                empty_fig.update_layout(
                    title='Gene Content Distribution (Error occurred)',
                    height=400
                )
                return empty_fig
        
        @self.app.callback(
            [Output('horizon-plot-iframe', 'srcDoc'),
             Output('horizon-plot-container', 'style')],
            [Input('project-select-horizon', 'value'),
             Input('subproject-select-horizon', 'value'),
             Input('taxa-count', 'value')],
            prevent_initial_call=True
        )
        def update_horizon_plot(selected_project, selected_subproject, taxa_count):
            """Update the horizon plot based on project selection and taxa count."""
            try:
                logger.info(f"update_horizon_plot callback triggered")
                
                # Default fallback content and style
                default_html = "<p>Loading horizon plot...</p>"
                
                # Ensure taxa_count is a valid number
                try:
                    taxa_count = int(taxa_count) if taxa_count else 20
                except (ValueError, TypeError):
                    taxa_count = 20
                
                # Calculate the container height based on taxa count
                container_height = 150 + (50 * taxa_count)
                container_style = {
                    "width": "100%", 
                    "margin": "0", 
                    "padding": "0",
                    "height": f"{container_height}px",
                    "position": "relative",
                    "z-index": 0,
                    "display": "block"
                }
                
                # Safety check for necessary attributes
                if not hasattr(self, 'horizon_app') or not hasattr(self, 'df'):
                    logger.error("Required attributes missing")
                    return default_html, container_style
                    
                # Create a safe copy of the dataframe to avoid modification issues
                try:
                    df = self.df.copy() if hasattr(self, 'df') else pd.DataFrame()
                except Exception as df_error:
                    logger.error(f"Error copying dataframe: {str(df_error)}")
                    return f"<p>Error accessing data: {str(df_error)}</p>", container_style
                
                # Apply filters if dataframe is not empty
                if not df.empty:
                    # Filter by project
                    if selected_project and selected_project != 'ALL':
                        try:
                            df = df[df['project_id'].astype(str) == str(selected_project)]
                        except Exception as e:
                            logger.error(f"Error filtering by project: {str(e)}")
                    
                    # Filter by subproject
                    if selected_subproject and selected_subproject != 'ALL' and 'subproject' in df.columns:
                        try:
                            df = df[df['subproject'].astype(str) == str(selected_subproject)]
                        except Exception as e:
                            logger.error(f"Error filtering by subproject: {str(e)}")
                
                # Check if filtered data is empty
                if df.empty:
                    return "<p>No data available for the selected filters</p>", container_style
                
                # Generate plot HTML
                try:
                    if hasattr(self.horizon_app, '_get_horizon_plot_html'):
                        html_content = self.horizon_app._get_horizon_plot_html(
                            df=df,
                            taxa_count=taxa_count or 20
                        )
                        
                        if not html_content:
                            return "<p>No content generated by plot function</p>", container_style
                            
                        return html_content, container_style
                    else:
                        return "<p>Horizon plot function not available</p>", container_style
                except Exception as plot_error:
                    logger.error(f"Error generating plot: {str(plot_error)}")
                    logger.error(traceback.format_exc())
                    return f"<p>Error generating plot: {str(plot_error)}</p>", container_style
                
            except Exception as e:
                import traceback
                logger.error(f"Unhandled error in update_horizon_plot: {str(e)}")
                logger.error(traceback.format_exc())
                return f"<p>Application error: {str(e)}</p>", {
                    "width": "100%", 
                    "margin": "0", 
                    "padding": "0",
                    "height": "700px",
                    "position": "relative",
                    "z-index": 0,
                    "display": "block"
                }
        
        # MAG callbacks
        @self.app.callback(
            [Output('mag-content', 'children', allow_duplicate=True),
             Output('mag-genome-warning', 'children'),
             Output('mag-genome-warning', 'style'),
             Output('mag-info-section', 'style')],
            [Input('mag-table-selected-rows', 'data')],
            [State('mag-table-data', 'children')],
            prevent_initial_call=True
        )
        def update_mag_content(selected_rows, table_data_json):
            """Update the MAG information content when a MAG is selected."""
            try:
                if not selected_rows or len(selected_rows) == 0:
                    # No MAG selected, show empty content
                    return (
                        "",
                        "Please select a MAG from the table to view detailed information.",
                        {"display": "block"},
                        {"display": "none"}
                    )
                
                # Parse the JSON table data
                table_data = json.loads(table_data_json)
                
                selected_row = table_data[selected_rows[0]]
                mag_name = selected_row.get('name')
                
                if not mag_name:
                    return (
                        html.Div("Selected MAG has no name"),
                        "Invalid MAG selection",
                        {"display": "block", "padding": "10px", "backgroundColor": "#f8d7da",  "borderRadius": "5px"},
                        {"display": "none"}
                    )
                
                logger.info(f"Updating MAG content for: {mag_name}")
                
                # Create simplified info section without accessing MAG data directly
                # This avoids large payload issues
                info_section = html.Div([
                    html.H4(f"MAG: {mag_name}", className="mb-3"),
                    
                    # Basic info card - use data already in table_data
                    html.Div([
                        html.H5("Basic Information", className="card-header"),
                        html.Div([
                            html.Table([
                                html.Tr([html.Th("Property"), html.Th("Value")], 
                                        style={"backgroundColor": "#f5f5f5"}),
                                html.Tr([html.Td("Taxonomy"), html.Td(selected_row.get('taxonomy', 'Unknown'))]),
                                html.Tr([html.Td("Completeness"), 
                                        html.Td(f"{selected_row.get('completeness', 0) * 100:.1f}%" if selected_row.get('completeness', 0) <= 1 else 
                                              f"{selected_row.get('completeness', 0):.1f}%")]),
                                html.Tr([html.Td("Contamination"), 
                                        html.Td(f"{selected_row.get('contamination', 0) * 100:.1f}%" if selected_row.get('contamination', 0) <= 1 else 
                                              f"{selected_row.get('contamination', 0):.1f}%")]),
                                html.Tr([html.Td("Quality"), 
                                        html.Td(selected_row.get('quality', 'Unknown').capitalize())]),
                                html.Tr([html.Td("Size"), html.Td(f"{selected_row.get('size_mb', 0):.2f} Mb")]),
                                html.Tr([html.Td("GC Content"), 
                                        html.Td(f"{selected_row.get('gc_content', 0) * 100:.1f}%" if selected_row.get('gc_content', 0) <= 1 else 
                                              f"{selected_row.get('gc_content', 0):.1f}%")]),
                                html.Tr([html.Td("Contigs"), html.Td(f"{selected_row.get('num_contigs', 0):,}")]),
                                html.Tr([html.Td("N50"), html.Td(f"{selected_row.get('n50', 0):,} bp")]),
                                html.Tr([html.Td("Longest Contig"), html.Td(f"{selected_row.get('longest_contig', 0):,} bp")])
                            ], className="table table-bordered", 
                              style={"width": "100%", "marginBottom": "0"})
                        ], className="card-body")
                    ], className="card mb-4"),
                    
                    # Instructions for user
                    html.Div([
                        html.H5("Genome Features", className="card-header"),
                    html.Div([
                            html.P("Select this MAG in the table to view its genome visualization and sequence data."),
                            html.P("Use the sequence selector dropdown to view specific contigs."),
                            html.P("Use the search box to find specific genes or functions.")
                        ], className="card-body")
                    ], className="card")
                ])
                
                # Check if sequence data might be available
                has_sequences = selected_row.get('num_contigs', 0) > 0
                
                # Set appropriate warning message
                if has_sequences:
                    return (
                        info_section,
                        "",
                        {"display": "none"},
                        {"display": "block"}
                    )
                else:
                    return (
                        info_section,
                        "No sequence data available for this MAG",
                        {"display": "block", "padding": "10px", "backgroundColor": "#fff3cd", "color": "#856404", "borderRadius": "5px"},
                        {"display": "block"}
                    )
                
            except Exception as e:
                logger.error(f"Error updating MAG content: {str(e)}", exc_info=True)
                return (
                    html.Div([
                        html.H5("Error Loading MAG Information", className="alert-heading"),
                        html.P(f"An error occurred: {str(e)}"),
                        html.Hr(),
                        html.P("Please try selecting a different MAG.")
                    ]),
                    "Error loading MAG data",
                    {"display": "block", "padding": "10px", "backgroundColor": "#f8d7da",  "borderRadius": "5px"},
                    {"display": "none"}
                )

                
        # Add a new callback to update the debug info
        @self.app.callback(
            Output('mags-debug-info', 'children'),
            [Input('mag-table', 'selected_rows'),
             Input('mag-table', 'data'),
             Input('sequence-selector', 'value')],
            prevent_initial_call=False
        )
        def update_debug_info(selected_rows, table_data, selected_sequence):
            """Update the debug information panel with current state."""
            try:
                now = datetime.now()
                timestamp = now.strftime('%Y-%m-%d %H:%M:%S')
                
                # Get user ID
                user_id = self.user_id
                
                # Check if MAGs app is initialized
                mags_app_initialized = hasattr(self, 'mags_app') and self.mags_app is not None
                
                # Get MAG directory
                mag_dir = getattr(self.mags_app, 'mag_dir', 'Unknown') if mags_app_initialized else 'Unknown'
                
                # Check if a MAG is selected
                mag_selected = selected_rows is not None and len(selected_rows) > 0 and table_data is not None
                
                # Get selected MAG info
                selected_mag_id = None
                selected_mag_name = None
                selected_mag_taxonomy = None
                
                if mag_selected and selected_rows[0] < len(table_data):
                    try:
                        selected_row = table_data[selected_rows[0]]
                        selected_mag_id = selected_row.get('bin_id', 'Unknown')
                        selected_mag_name = selected_row.get('name', selected_mag_id)
                        selected_mag_taxonomy = selected_row.get('taxonomy', 'Unknown')
                    except (IndexError, KeyError, TypeError) as e:
                        logger.error(f"Error getting MAG details: {str(e)}")
                        selected_mag_id = "Error"
                        selected_mag_name = f"Error: {str(e)}"
                        selected_mag_taxonomy = "Unknown"
                
                # Check NGCircos status
                ngcircos_available = "Unknown"
                try:
                    ngcircos_status = dcc.Store(id="ngcircos-status", data="checking")
                    self.app.clientside_callback(
                        """
                        function(data) {
                            return window.NGCircos ? "Available" : "Not Available";
                        }
                        """,
                        Output("ngcircos-status", "data", allow_duplicate=True),
                        [Input("ngcircos-status", "data")],
                        prevent_initial_call=True
                    )
                    ngcircos_available = "Checking..."
                except Exception as e:
                    logger.error(f"Error checking NGCircos: {str(e)}")
                    ngcircos_available = f"Error checking: {str(e)}"
                
                # Create debug info text
                debug_info = f"""Time: {timestamp}
                User ID: {user_id}
                MAGs app initialized: {mags_app_initialized}
                MAG directory: {mag_dir}
                
                MAG selection active: {'Yes' if mag_selected else 'No'}
                Selected MAG: {selected_mag_id if mag_selected else 'None'}
                Selected MAG name: {selected_mag_name if mag_selected else 'None'}
                Selected MAG taxonomy: {selected_mag_taxonomy if mag_selected else 'None'}
                Selected sequence: {selected_sequence if selected_sequence else 'None'}
                
                NGCircos status: {ngcircos_available}
                
                Note: MAGs are loaded from the database, not from the MAG directory.
                The MAG directory is only used for temporary files.
                
                If no MAGs are displayed, please check:
                1. That the user ID is correct
                2. That there are MAGs in the database for this user
                3. Check the browser console for errors
                
                If sequence viewer is not showing:
                1. Check that a MAG is selected in the table
                2. Check that the MAG has sequence data
                3. Open browser console (F12) and check for errors
                """
                
                return html.Pre(
                    debug_info,
                    style={
                        'whiteSpace': 'pre-wrap',
                        'fontSize': '12px',
                        'fontFamily': 'monospace',
                        'backgroundColor': '#f8f9fa',
                        'padding': '10px',
                        'border': '1px solid #ddd',
                        'borderRadius': '5px',
                        'maxHeight': '300px',
                        'overflow': 'auto'
                    }
                )

            except Exception as e:
                logger.error(f"Error updating debug info: {str(e)}")
                logger.error(traceback.format_exc())
            
                return html.Pre(
                    f"Error updating debug info: {str(e)}\n\n{traceback.format_exc()}",
                    style={
                        'whiteSpace': 'pre-wrap',
                        'fontSize': '12px',
                        'fontFamily': 'monospace',
                        'color': 'red',
                        'backgroundColor': '#f8f9fa',
                        'padding': '10px',
                        'border': '1px solid #ddd',
                        'borderRadius': '5px',
                        'maxHeight': '300px',
                        'overflow': 'auto'
                    }
                )

        @self.app.callback(
            Output('average-read-length-plot', 'figure'),
            [Input('sample-dropdown-qc', 'value')],
            prevent_initial_call=True
        )
        def update_qc_plots(selected_sample):
            if not selected_sample:
                empty_fig = go.Figure()
                empty_fig.update_layout(
                    title="No sample selected",
                    annotations=[{
                        'text': "Please select a sample to display QC plots",
                        'xref': 'paper',
                        'yref': 'paper',
                        'showarrow': False,
                        'font': {'size': 14}
                    }]
                )
                return empty_fig

            try:
                # Get QC plots with improved styling
                qc_plots = self.qc_app.update_qc_plots(selected_sample)
                distribution_plots = self.qc_app.create_quality_distribution_plots(selected_sample)
                
                # Create subplots
                fig = make_subplots(
                    rows=2, cols=2,
                    subplot_titles=(
                        "Read Length Distribution",
                        "GC Content Distribution",
                        "Cumulative Bases",
                        "Quality Score Distribution"
                    )
                )
                
                # Add traces to appropriate subplots
                for i, plot in enumerate(qc_plots + distribution_plots):
                    row = (i // 2) + 1
                    col = (i % 2) + 1
                    
                    for trace in plot.data:
                        # Fix for GC content plot (normalize and smooth)
                        if i == 1:  # GC Content plot
                            if hasattr(trace, 'y') and len(trace.y) > 0:
                                # Normalize
                                trace.y = trace.y / np.sum(trace.y)
                                # Smooth using moving average
                                trace.y = pd.Series(trace.y).rolling(window=5, center=True).mean()
                        
                        # Fix for cumulative bases plot
                        elif i == 2:  # Cumulative Bases plot
                            if hasattr(trace, 'y') and len(trace.y) > 0:
                                # Ensure monotonically increasing
                                trace.y = np.cumsum(np.abs(trace.y))
                        
                        fig.add_trace(trace, row=row, col=col)
                
                # Check if this is likely 16S data (average read length around 1500 bp)
                avg_read_length = self.qc_app.get_average_read_length(selected_sample)
                if 1300 <= avg_read_length <= 1700:  # Typical 16S range
                    # Add 16S reference line to read length plot
                    fig.add_hline(
                        y=1500,
                        line=dict(color="red", width=2, dash="dash"),
                        annotation_text="Typical 16S length (1500 bp)",
                        annotation_position="top right",
                        row=1, col=1
                    )
                
                # Update layout with better styling
                fig.update_layout(
                    template="plotly_white",
                    showlegend=True,
                    legend=dict(
                        yanchor="bottom",
                        y=1.02,
                        xanchor="right",
                        x=1
                    ),
                    margin=dict(l=50, r=50, t=100, b=50),
                    font=dict(family="Arial", size=12),
                    height=800,
                    title=dict(
                        text=f"QC Metrics for Sample: {selected_sample}",
                        x=0.5,
                        xanchor="center"
                    )
                )
                
                # Update axes for better readability
                fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')
                fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')
                
                # Specific updates for each subplot
                fig.update_xaxes(title_text="Read Length (bp)", row=1, col=1)
                fig.update_yaxes(title_text="Count", row=1, col=1)
                
                fig.update_xaxes(title_text="GC Content (%)", row=1, col=2)
                fig.update_yaxes(title_text="Frequency", row=1, col=2)
                
                fig.update_xaxes(title_text="Read Length (bp)", row=2, col=1)
                fig.update_yaxes(title_text="Cumulative Bases", row=2, col=1)
                
                fig.update_xaxes(title_text="Quality Score", row=2, col=2)
                fig.update_yaxes(title_text="Count", row=2, col=2)
                
                return fig
                
            except Exception as e:
                print(f"Error in update_qc_plots: {str(e)}")
                empty_fig = go.Figure()
                empty_fig.update_layout(
                    title="Error loading QC plots",
                    annotations=[{
                        'text': str(e),
                        'xref': 'paper',
                        'yref': 'paper',
                        'showarrow': False,
                        'font': {'size': 14}
                    }]
                )
                return empty_fig

        @self.app.callback(
            [Output('read-length-dist', 'figure'),
             Output('quality-dist', 'figure'),
             Output('gc-dist', 'figure'),
             Output('base-composition', 'figure')],
            [Input('sample-selector', 'value')]
        )
        def update_qc_sample_plots(selected_sample):
            """Update QC plots for the selected sample."""
            if not hasattr(self, 'user_id'):
                return [go.Figure() for _ in range(4)]

            try:
                # Use the cached QC instance
                qc_instance = self.qc_app
                
                # If no sample is selected, use the first available sample
                if not selected_sample and hasattr(qc_instance, 'stats_df') and not qc_instance.stats_df.empty:
                    selected_sample = qc_instance.stats_df['sample_name'].iloc[0]
                
                if not selected_sample:
                    return [go.Figure() for _ in range(4)]

                # Let QC class handle the plot creation with caching
                return qc_instance.update_qc_plots(selected_sample)
                
            except Exception as e:
                logger.error(f"Error in update_qc_sample_plots: {str(e)}")
                return [go.Figure() for _ in range(4)]

        @self.app.callback(
            Output("analysis-output", "children"),
            Input("analyze-button", "n_clicks"),
            [State("sample-select-value-tax", "value")],
            prevent_initial_call=True
        )
        def analyze_taxonomy_button_clicked(n_clicks, selected_samples):
            if n_clicks is None:
                raise PreventUpdate
            
            try:
                # Filter data for selected samples
                df_to_analyze = self.taxonomy_app.df_sorted.copy()
                if selected_samples:
                    df_to_analyze = df_to_analyze[df_to_analyze['sample_id'].isin(selected_samples)]
                
                if df_to_analyze.empty:
                    return "No data found for selected samples"
                
                # Generate analysis
                analysis_json = self.taxonomy_app.analyze_taxonomy_data(df_to_analyze)
                
                return analysis_json
                
            except Exception as e:
                error_msg = {
                    "error": str(e),
                    "traceback": traceback.format_exc(),
                    "debug_info": {
                        "selected_samples": selected_samples,
                        "df_shape": self.taxonomy_app.df_sorted.shape if hasattr(self.taxonomy_app, 'df_sorted') else None
                    }
                }
                return json.dumps(error_msg)

        @self.app.callback(
            Output('date-range-picker', 'disabled'),
            [Input('use-date-range', 'checked')]
        )
        def toggle_date_picker(use_date_range):
            try:
                return not use_date_range
            except Exception as e:
                logger.error(f"Error toggling date picker: {str(e)}")
                logger.error(traceback.format_exc())
                return True  # Default to enabled if there's an error

        @self.app.callback(
            Output("download-pdf", "data"),
            Input("download-pdf-button", "n_clicks"),
            [State('project-select-horizon', 'value'),
             State('subproject-select-horizon', 'value'),
             State('taxa-count', 'value')],
            prevent_initial_call=True
        )
        def download_horizon_pdf(n_clicks, selected_project, selected_subproject, taxa_count):
            if not n_clicks:
                raise PreventUpdate

            try:
                print("\n=== Debug: download_horizon_pdf ===")
                print(f"Generating PNG for project: {selected_project}")
                
                # Get filtered data
                df = self.horizon_app.df.copy()
                if selected_project and selected_project != 'ALL':
                    df = df[df['project_id'].astype(str).str.strip() == str(selected_project).strip()]
                if selected_subproject and selected_subproject != 'ALL':
                    df = df[df['subproject_id'].astype(str).str.strip() == str(selected_subproject).strip()]
                
                
                if len(df) == 0:
                    raise PreventUpdate
                
                # Generate plot HTML
                try:
                    if hasattr(self.horizon_app, '_get_horizon_plot_html'):
                        html_content = self.horizon_app._get_horizon_plot_html(
                            df=df,
                            taxa_count=taxa_count or 20
                        )
                        
                        if not html_content:
                            return "<p>No content generated by plot function</p>", container_style
                            
                        return html_content, container_style
                    else:
                        return "<p>Horizon plot function not available</p>", container_style
                except Exception as plot_error:
                    logger.error(f"Error generating plot: {str(plot_error)}")
                    logger.error(traceback.format_exc())
                    return f"<p>Error generating plot: {str(plot_error)}</p>", container_style
                
            except Exception as e:
                print(f"Error in download_horizon_pdf: {str(e)}")
                traceback.print_exc()
                raise PreventUpdate

        
        
        

        
        @self.app.callback(
            Output('circos-plot-iframe', 'srcDoc'),
            [Input('mag-table', 'selected_rows')],
            [State('mag-table', 'data')],
            prevent_initial_call=True
        )
        def update_circos_plot_callback(selected_rows, table_data):
            """Update the circos plot for the selected MAG."""
            start_time = time.time()
            logger.info(f"Circos plot callback triggered with selected_rows: {selected_rows}")
            
            try:
                if not selected_rows or not table_data or len(selected_rows) == 0:
                    return html.Div("Select a MAG from the table to view its genome visualization", 
                                   style={'padding': '20px', 'textAlign': 'center'})
                
                selected_row = table_data[selected_rows[0]]
                mag_name = selected_row.get('name', 'unknown')
                bin_id = selected_row.get('bin_id', 'unknown')
                
                # Log detailed information about the selected MAG
                logger.info(f"Creating circos plot for MAG: {mag_name} (bin_id: {bin_id})")
                
                # Generate the circos plot HTML content from the mags app, passing the selected MAG ID
                try:
                    circos_html = self.mags_app._create_test_circos_plot(mag_id=bin_id)
                    return circos_html
                except Exception as e:
                    logger.error(f"Error generating circos plot HTML: {str(e)}", exc_info=True)
                    # Fallback to basic plot without MAG ID
                    circos_html = self.mags_app._create_test_circos_plot(mag_id=None)
                
                # Create the iframe with appropriate settings for HTML content
                
                
                
                # Alternative approach - also provide a direct HTML component for browsers where iframe might not work properly
                direct_html_div = html.Div([
                    html.Div(id="direct-circos-container", style={'width': '100%', 'height': '700px'}),
                    # Use dcc.Markdown to inject raw HTML safely
                    dcc.Markdown(circos_html, dangerously_allow_html=True, 
                                 style={'display': 'none'}, id="raw-circos-html"),
                    # JavaScript to move content into the container div
                    html.Script('''
                    document.addEventListener("DOMContentLoaded", function() {
                        try {
                            const rawHtml = document.getElementById("raw-circos-html").textContent;
                            const container = document.getElementById("direct-circos-container");
                            if (container) {
                                container.innerHTML = rawHtml;
                                console.log("Injected direct HTML for circos plot");
                            } else {
                                console.error("Container for direct HTML not found");
                            }
                        } catch(e) {
                            console.error("Error injecting direct HTML:", e);
                        }
                    });
                    ''')
                ], style={'display': 'none'})  # Keep hidden initially as a fallback
                
                # Create container with header and plot
                circos_container = html.Div([
                    html.H4(f"Genome Visualization for {mag_name}", 
                           style={'marginBottom': '10px', 'paddingLeft': '10px'}),
                    html.Div([
                        html_iframe,
                        direct_html_div  # Include the alternative approach as a fallback
                    ], style={
                        'width': '100%',
                        'marginBottom': '20px'
                    }),
                    html.Div([
                        html.P("Note: The plot shows contigs (chromosomes) and gene annotations from the MAG's GFF file if available.", 
                              style={'color': '#666', 'fontSize': '12px', 'fontStyle': 'italic'})
                    ], style={'padding': '10px'})
                ])
                
                elapsed_time = time.time() - start_time
                logger.info(f"Successfully created circos plot for {mag_name} in {elapsed_time:.2f} seconds")
                return circos_container
                
            except Exception as e:
                logger.error(f"Error in circos plot callback: {str(e)}")
                logger.error(traceback.format_exc())
                
                # Return a simple fallback plot with error information
                try:
                    # Attempt to create a simple default plot
                    default_html = """
                        <html>
                        <body style="text-align: center; font-family: Arial, sans-serif;">
                            <div style="margin: 50px auto; max-width: 500px; padding: 20px; border: 1px solid #ddd; border-radius: 5px; background-color: #f9f9f9;">
                                <h3>Genome Visualization</h3>
                                <p>There was an error generating the genome visualization.</p>
                                <p style="font-size: 12px; color: #666;">Error details have been logged for the administrator.</p>
                                <div style="width: 400px; height: 400px; margin: 20px auto; border: 1px solid #ccc; border-radius: 50%; background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);">
                                    <div style="padding-top: 180px;">Genome Visualization Placeholder</div>
                                </div>
                            </div>
                        </body>
                        </html>
                    """
                    
                    fallback_iframe = html.Iframe(
                        srcDoc=default_html,
                        style={
                            'width': '100%',
                            'height': '600px',
                            'border': '1px solid #ddd',
                            'borderRadius': '5px',
                            'backgroundColor': '#ffffff',
                        },
                        sandbox='allow-scripts'
                    )
                    
                    return html.Div([
                        html.H4("Genome Visualization (Error Occurred)", 
                               style={'marginBottom': '10px', 'paddingLeft': '10px', 'color': '#c62828'}),
                        html.Div([
                            fallback_iframe
                        ], style={
                            'width': '100%',
                            'marginBottom': '20px'
                        })
                    ])
                    
                except Exception:
                    # Last resort fallback
                    return html.Div([
                        html.H4("Error Creating Genome Visualization", style={'color': 'red'}),
                        html.P(f"An error occurred: {str(e)}"),
                        html.P("Try refreshing the page or selecting a different MAG.")
                    ], style={'padding': '20px'})

        
        # Remove the duplicated callback that serves no purpose
        # self.app.clientside_callback(
        #     """
        #     function(trigger) {
        #         // Simple function that just returns empty content
        #         // We keep this to ensure the scripts get loaded when needed
        #         console.log("Secondary NGCircos callback triggered");
        #         return [];
        #     }
        #     """,
        #     Output('global-error-container', 'children', allow_duplicate=True),
        #     Input('circos-plot', 'children'),
        #     prevent_initial_call=True
        # )

        # Add missing components to suppress errors
        self.app.layout.children[0].children.append(html.Div([
            html.Button("Search", id="gene-search-button", n_clicks=0, style={"display": "none"}),
            html.Button("Copy FASTA", id="copy-fasta-btn", n_clicks=0, style={"display": "none"}),
            dcc.Input(id="gene-search", type="text", style={"display": "none"}),
            html.Div("", id="sequence-data", style={"display": "none"}),
            html.Div("", id="fasta-data", style={"display": "none"}),
        ], id="suppress-missing-components", style={"display": "none"}))
        
        # Modify the gene search callback to use PreventUpdate to avoid unnecessary processing
        @self.app.callback(
            Output('gene-search-results', 'children'),
            [Input('gene-search-button', 'n_clicks'),
             Input('mag-table', 'selected_rows')],
            [State('gene-search', 'value'),
             State('mag-table', 'data')],
            prevent_initial_call=True
        )
        def display_search_results(n_clicks, selected_rows, search_term, table_data):
            """Display search results for gene search."""
            if not n_clicks or not search_term or not selected_rows or not table_data:
                return html.Div("Enter a search term and select a MAG to search for genes.", 
                               style={'fontStyle': 'italic', 'color': '#666'})
                
            try:
                # Get selected MAG data
                selected_mag = table_data[selected_rows[0]]
                mag_name = selected_mag.get("name") or selected_mag.get("Name")
                bin_id = selected_mag.get("bin_id")
                
                if not mag_name or not bin_id:
                    return html.Div("Could not identify selected MAG.", 
                                   style={'color': 'red'})
                
                # Get GFF data for the MAG
                try:
                    from users.models import MAG
                    mag = MAG.objects.get(name=bin_id, user_id=self.user_id)
                    gff_file = mag.gff_file
                    
                    if not gff_file:
                        return html.Div("No annotation data available for this MAG.", 
                                       style={'color': 'orange'})
                    
                    # Search for the term in GFF data
                    search_term = search_term.lower()
                    matched_features = []
                    
                    for line in gff_file.split('\n'):
                        if line.startswith('#') or not line.strip():
                            continue
                            
                        if search_term in line.lower():
                            matched_features.append(line)
                    
                    # Format and return results
                    if not matched_features:
                        return html.Div(f"No genes found matching '{search_term}'",
                                       style={'color': '#666'})
                    
                    result_items = [html.H5(f"Found {len(matched_features)} matches:")]
                    
                    for i, feature in enumerate(matched_features[:10]):  # Limit to first 10
                        parts = feature.split('\t')
                        if len(parts) >= 9:
                            feature_id = parts[8].split(';')[0].replace('ID=', '')
                            result_items.append(html.Div([
                                html.Strong(f"{i+1}. {feature_id}"),
                                html.Span(f" - {parts[2]} at {parts[0]}:{parts[3]}-{parts[4]}")
                            ]))
                    
                    if len(matched_features) > 10:
                        result_items.append(html.Div(f"...and {len(matched_features) - 10} more matches",
                                                   style={'fontStyle': 'italic'}))
                    
                    return html.Div(result_items)
                
                except MAG.DoesNotExist:
                    return html.Div(f"MAG '{bin_id}' not found in database.", 
                                   style={'color': 'red'})
                except Exception as e:
                    return html.Div(f"Error searching genes: {str(e)}", 
                                   style={'color': 'red'})
                
            except Exception as e:
                import traceback
                logger.error(f"Error in gene search: {str(e)}")
                logger.error(traceback.format_exc())
                return html.Div(f"Error: {str(e)}", style={'color': 'red'})

        # Set default sequence when options change
        @self.app.callback(
            Output('sequence-selector', 'value'),
            [Input('sequence-selector', 'options')],
            prevent_initial_call=True,allow_duplicate=True
        )
        def set_default_sequence(options):
            """Set the default sequence to the first option when options change."""
            try:
                if not options or len(options) == 0:
                    logger.warning("No sequence options available")
                    return None
                logger.info(f"Setting default sequence to: {options[0]['value']}")
                return options[0]['value']
            except Exception as e:
                logger.error(f"Error setting default sequence: {str(e)}")
                return None

        # Add navigation callback
        @self.app.callback(
            [Output("app-subtitle", "children"),
             Output("page-content", "children"),
             Output("header-title", "children")],
            [Input("url", "pathname")],
            prevent_initial_call=False
        )
        def update_page_content(pathname):
            if pathname is None:
                return "Welcome", html.Div("Please select an app from the menu."), ""
                
            if pathname in self._apps:
                app_info = self._apps[pathname]
                try:
                    app_instance = app_info['instance']
                    # Add debug logging
                    logger.info(f"Loading app {app_info['name']} with instance {app_instance}")
                    
                    if hasattr(app_instance, 'layout'):
                        app_layout = app_instance.layout
                    elif hasattr(app_instance, 'app'):
                        app_layout = app_instance.app.layout
                    else:
                        logger.error(f"App {app_info['name']} has no layout")
                        return "Error", html.Div(f"Error: App {app_info['name']} has no layout"), ""
                    
                    logger.info(f"Successfully loaded app: {app_info['name']}")
                    return app_info['name'], app_layout, ""
                except Exception as e:
                    logger.error(f"Error loading app {app_info['name']}: {str(e)}")
                    return "Error", html.Div(f"Error loading app: {str(e)}"), ""
            else:
                return "Welcome", html.Div("Please select an app from the menu."), ""

        # Correlations app callbacks
        @self.app.callback(
            Output('download-metadata-csv', 'data'),
            [Input('download-csv-button', 'n_clicks')],
            prevent_initial_call=True
        )
        def download_metadata_template(n_clicks):
            if n_clicks is None:
                raise PreventUpdate
            try:
                logger.info("Downloading metadata template")
                return self.correlations_app.handle_template_download()
            except Exception as e:
                logger.error(f"Error downloading metadata template: {str(e)}")
                logger.error(traceback.format_exc())
                raise PreventUpdate
            
        @self.app.callback(
            Output('download-corr-csv', 'data'),
            [Input('btn-download-corr', 'n_clicks')],
            prevent_initial_call=True
        )
        def download_correlations(n_clicks):
            if n_clicks is None:
                raise PreventUpdate
            try:
                logger.info("Downloading correlations data")
                return self.correlations_app.handle_correlations_download()
            except Exception as e:
                logger.error(f"Error downloading correlations data: {str(e)}")
                logger.error(traceback.format_exc())
                raise PreventUpdate
            
        @self.app.callback(
            Output('notification-output', 'children', allow_duplicate=True),
            [Input('upload-data', 'contents')],
            [State('upload-data', 'filename')],
            prevent_initial_call=True
        )
        def handle_upload(contents, filename):
            if contents is None:
                raise PreventUpdate
            try:
                logger.info(f"Handling upload for {filename}")
                return self.correlations_app.handle_file_upload(contents, filename)
            except Exception as e:
                logger.error(f"Error handling file upload: {str(e)}")
                logger.error(traceback.format_exc())
                return dmc.Alert(
                    id="upload-error",
                    title="Error",
                    children=f"Error uploading file: {str(e)}",
                    color="red",
                    withCloseButton=True,
                    variant="filled",
                    style={"marginBottom": "15px"}
                )

        @self.app.callback(
            Output('download-svg-diversity', 'data'),
            [Input('btn-download-svg-correlations', 'n_clicks')],
            [State('heatmap-dendrogram', 'figure')],
            prevent_initial_call=True
        )
        def download_svg(n_clicks, figure):
            if n_clicks is None or figure is None:
                raise PreventUpdate
                
            try:
                logger.info("Downloading SVG for correlations heatmap")
                import plotly.io as pio
                svg_figure = pio.to_image(figure, format='svg', width=1200, height=800)
                filename = f"correlations_heatmap_{datetime.now().strftime('%Y%m%d_%H%M%S')}.svg"
                return dict(content=svg_figure, filename=filename)
            except Exception as e:
                logger.error(f"Error generating SVG: {str(e)}")
                raise PreventUpdate

        @self.app.callback(
            [Output("graph1", "figure"),
             Output("download-plot", "data")],
            [Input("legend-items-input", "value"),
             Input("export-svg-button", "n_clicks"),
             Input("plot-type-dropdown-tax", "value"),
             Input("tax-rank-dropdown-tax", "value"),
             Input("use-date-checkbox-tax", "checked"),
             Input("sample-select-value-tax", "value"),
             
             Input("legend-toggle-checkbox-tax", "checked")],
            [State("width-input", "value"),
             State("height-input", "value"),
             State("graph1", "figure")],
            prevent_initial_call=False
        )
        def update_and_export_plot(n_items, n_clicks, plot_type, tax_rank, use_date, selected_samples, show_legend, width, height, figure):
            """Update the plot based on user selections and export if requested."""
            try:
                ctx = dash.callback_context
                triggered_id = ctx.triggered[0]['prop_id'].split('.')[0] if ctx.triggered else None
                
                logger.info(f"Plot update triggered by: {triggered_id} with {len(selected_samples) if selected_samples else 0} samples selected")
                
                # Initialize download_data to no_update
                download_data = no_update
                
                
                # Prevent update if critical inputs are missing
                if not plot_type or not tax_rank:
                    logger.warning("Missing plot type or taxonomic rank")
                    return no_update, no_update
                
                # Handle case when no samples are selected
                if not selected_samples or len(selected_samples) == 0:
                    logger.info("No samples selected, showing empty plot with message")
                    empty_fig = go.Figure()
                    empty_fig.add_annotation(
                        text="No samples selected - please select samples to display data",
                        xref="paper", yref="paper",
                        x=0.5, y=0.5, showarrow=False,
                        font=dict(size=20)
                    )
                    return empty_fig, no_update
                
                # Filter dataframe based on selected samples
                df_filtered = self.taxonomy_app.df[self.taxonomy_app.df['sample_id'].isin(selected_samples)]
                
                # Debug: Print dataframe info to browser console
                
                # If triggered by sample selection, always create a new plot
                if triggered_id == "sample-select-value-tax":
                    logger.info("Sample selection changed, creating completely new plot")
                    
                    # DEBUG: Log what dataframe is being used
                    logger.info(f"Using filtered dataframe with shape: {df_filtered.shape}")
                    logger.info(f"Sample IDs in filtered data: {df_filtered['sample_id'].unique().tolist()}")
                    
                    # Create new plot based on selected type
                    if plot_type == 'stackedbar':
                        # NOTE: This one uses self.df_sorted instead of df_filtered - this might be the bug
                        logger.info(f"Stacked bar plot using self.df_sorted with shape: {self.df_sorted.shape}")
                        logger.info(f"Sample IDs in self.df_sorted: {self.df_sorted['sample_id'].unique().tolist()}")
                        
                        # Use filtered data instead of self.df_sorted
                        fig = self.taxonomy_app.plot_stacked_bar(df_filtered, use_date, tax_rank)
                    elif plot_type == 'groupedbar':
                        fig = self.taxonomy_app.plot_grouped_bar(df_filtered, use_date, tax_rank)
                    elif plot_type == 'area':
                        fig = self.taxonomy_app.plot_area(df_filtered, use_date, tax_rank)
                    elif plot_type == 'line':
                        fig = self.taxonomy_app.plot_line(df_filtered, use_date, tax_rank)
                    elif plot_type == 'heatmap':
                        fig = self.taxonomy_app.plot_heatmap(df_filtered, use_date, tax_rank)
                    elif plot_type == 'scatter':
                        fig = self.taxonomy_app.plot_scatter(df_filtered, use_date, tax_rank)
                    elif plot_type == 'pie' and len(selected_samples) == 1:
                        # For pie chart, we need exactly one sample
                        fig = self.taxonomy_app.plot_pie(df_filtered, selected_samples[0], tax_rank)
                    else:
                        # Default to stacked bar if invalid type
                        fig = self.taxonomy_app.plot_stacked_bar(df_filtered, use_date, tax_rank)
                
                # If triggered by export button
                elif triggered_id == "export-svg-button" and n_clicks:
                    logger.info("Export button clicked, preparing SVG file")
                    try:
                        current_time = datetime.now().strftime("%Y%m%d-%H%M%S")
                        filename = f"taxonomy_plot_{current_time}.svg"
                        
                        fig = figure  # Use current figure for export
                        if width and height:
                            fig['layout']['width'] = width
                            fig['layout']['height'] = height
                        
                        # Convert the figure to SVG format
                        svg_data = pio.to_image(fig, format='svg')
                        
                        # Set up the download with the SVG data
                        download_data = dict(content=svg_data.decode('utf-8'), filename=filename)
                        return figure, download_data
                    except Exception as e:
                        logger.error(f"Error exporting plot: {str(e)}", exc_info=True)
                        return no_update, no_update
                
                # For all other triggers, create a new plot
                else:
                    # DEBUG: Log what dataframe is being used
                    logger.info(f"Using filtered dataframe with shape: {df_filtered.shape}")
                    logger.info(f"Sample IDs in filtered data: {df_filtered['sample_id'].unique().tolist()}")
                    
                    # Create new plot based on selected type
                    if plot_type == 'stackedbar':
                        # NOTE: This one uses self.df_sorted instead of df_filtered - this might be the bug
                        logger.info(f"Stacked bar plot using self.df_sorted with shape: {self.df_sorted.shape}")
                        logger.info(f"Sample IDs in self.df_sorted: {self.df_sorted['sample_id'].unique().tolist()}")
                        
                        # Use filtered data instead of self.df_sorted
                        fig = self.taxonomy_app.plot_stacked_bar(df_filtered, use_date, tax_rank)
                    elif plot_type == 'groupedbar':
                        fig = self.taxonomy_app.plot_grouped_bar(df_filtered, use_date, tax_rank)
                    elif plot_type == 'area':
                        fig = self.taxonomy_app.plot_area(df_filtered, use_date, tax_rank)
                    elif plot_type == 'line':
                        fig = self.taxonomy_app.plot_line(df_filtered, use_date, tax_rank)
                    elif plot_type == 'heatmap':
                        fig = self.taxonomy_app.plot_heatmap(df_filtered, use_date, tax_rank)
                    elif plot_type == 'scatter':
                        fig = self.taxonomy_app.plot_scatter(df_filtered, use_date, tax_rank)
                    elif plot_type == 'pie' and len(selected_samples) == 1:
                        # For pie chart, we need exactly one sample
                        fig = self.taxonomy_app.plot_pie(df_filtered, selected_samples[0], tax_rank)
                    else:
                        # Default to stacked bar if invalid type
                        fig = self.taxonomy_app.plot_stacked_bar(df_filtered, use_date, tax_rank)
                
                
                if show_legend is not None:
                    fig.update_layout(showlegend=show_legend)
                
                # Update legend visibility based on abundance if n_items is set
                if n_items:
                    logger.info(f"Updating legend visibility for top {n_items} items")
                    trace_abundances = []
                    for i, trace in enumerate(fig.data):
                        try:
                            total = sum(trace.y) if hasattr(trace, 'y') and trace.y is not None else 0
                            trace_abundances.append((total, i))
                        except Exception as e:
                            logger.error(f"Error calculating abundance for trace {i}: {str(e)}")
                            trace_abundances.append((0, i))
                    
                    # Sort by abundance and get indices of traces to hide in legend
                    # trace_abundances.sort(key=lambda x: x[0], reverse=True)
                    # hidden_indices = {idx for _, idx in trace_abundances[n_items:]}
                    
                    # Update showlegend for each trace
                    # for i in range(len(fig.data)):
                    #     if i < len(fig.data) and hasattr(fig.data[i], 'showlegend'):
                    #         fig.data[i].showlegend = i not in hidden_indices
                
                logger.info("Successfully created plot")
                return fig, no_update
                
            except Exception as e:
                logger.error(f"Error in update_and_export_plot: {str(e)}", exc_info=True)
                empty_fig = go.Figure()
                empty_fig.add_annotation(
                    text=f"An error occurred: {str(e)}",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5, showarrow=False,
                    font=dict(size=20)
                )
                return empty_fig, no_update


        
        @self.app.callback(
            [Output("sample-select-value-tax", "value"),
             Output("project-dropdown-tax", "value"),
             Output("subproject-dropdown-tax", "value"),
             Output("date-dropdown-tax", "value")],
            [Input("project-dropdown-tax", "value"),
             Input("subproject-dropdown-tax", "value"),
             Input("date-dropdown-tax", "value"),
             Input("refresh-taxonomy-button", "n_clicks")],
            prevent_initial_call=True
        )
        def update_taxonomy_sample_selection(project, subproject, date, n_clicks):
            """Update sample selection based on project, subproject, and date filters."""
            try:
                ctx = dash.callback_context
                if not ctx.triggered:
                    raise PreventUpdate
                    
                # Log which input triggered the callback
                trigger_id = ctx.triggered[0]['prop_id'].split('.')[0] if ctx.triggered else None
                logger.info(f"Sample selection triggered by: {trigger_id}")
                    
                # If refresh button was clicked, clear the taxonomy app cache and reload data
                if trigger_id == "refresh-taxonomy-button":
                    logger.info("Refresh button clicked, reloading taxonomy data")
                    self.taxonomy_app._clear_cache()
                    # Get the fresh dataframe
                    df = self.taxonomy_app.df.copy()
                    logger.info(f"Reloaded taxonomy data with {len(df)} records")
                    logger.info(f"Unique sample IDs after reload: {df['sample_id'].unique().tolist()}")
                else:
                    df = self.taxonomy_app.df.copy()
                
                # Apply filters if they exist and aren't 'ALL'
                if project and project != 'ALL':
                    df = df[df['project_id'].astype(str) == str(project)]
                if subproject and subproject != 'ALL':
                    df = df[df['subproject'].astype(str) == str(subproject)]
                if date and date != 'ALL':
                    df = df[df['date'].astype(str) == str(date)]
                    
                # Get unique sample IDs after filtering
                filtered_samples = df['sample_id'].unique().tolist()
                
                # Log the result for debugging
                logger.info(f"Filtered samples: {len(filtered_samples)} samples found")
                logger.info(f"Filtered sample IDs: {filtered_samples}")
                
                # Return the filtered sample list and preserve the filter values
                return filtered_samples, project, subproject, date
                
            except Exception as e:
                logger.error(f"Error in update_taxonomy_sample_selection: {str(e)}")
                logger.error(traceback.format_exc())
                return [], None, None, None

        
        # Add diversity app callbacks
        @self.app.callback(
            [Output('sample-select-value-diversity', 'options'),
             Output('project-dropdown-diversity', 'options'),
             Output('subproject-dropdown-diversity', 'options'),
             Output('date-dropdown-diversity', 'options')],
            [Input('url', 'pathname')]
        )
        def update_diversity_dropdowns(pathname):
            try:
                logger.info(f"update_diversity_dropdowns called with pathname: {pathname}")
                
                # Return empty lists if not on diversity page or no diversity app
                if not pathname or '/diversity' not in pathname:
                    return [], [], [], []

                if not hasattr(self, 'diversity_app') or self.diversity_app is None:
                    logger.warning("Diversity app not initialized")
                    return [], [], [], []
                
                # Initialize with default empty options
                sample_options = []
                project_options = [{'label': 'All Projects', 'value': 'ALL'}]
                subproject_options = [{'label': 'All Subprojects', 'value': 'ALL'}]
                date_options = [{'label': 'All Dates', 'value': 'ALL'}]
                
                # Safely get unique sample IDs
                if hasattr(self.diversity_app, 'unique_sample_ids') and self.diversity_app.unique_sample_ids is not None:
                    try:
                        sample_options = [{'label': str(s), 'value': str(s)} for s in self.diversity_app.unique_sample_ids]
                        logger.info(f"Found {len(sample_options)} sample options")
                    except Exception as e:
                        logger.error(f"Error processing sample IDs: {str(e)}")
                
                # Safely get unique projects
                if hasattr(self.diversity_app, 'unique_projects_ids') and self.diversity_app.unique_projects_ids is not None:
                    try:
                        project_options.extend([
                            {'label': str(p), 'value': str(p)} 
                            for p in self.diversity_app.unique_projects_ids if p is not None
                        ])
                        logger.info(f"Found {len(project_options)-1} project options")
                    except Exception as e:
                        logger.error(f"Error processing project IDs: {str(e)}")
                
                # Safely get unique subprojects
                if hasattr(self.diversity_app, 'unique_subprojects') and self.diversity_app.unique_subprojects is not None:
                    try:
                        subproject_options.extend([
                            {'label': str(s), 'value': str(s)} 
                            for s in self.diversity_app.unique_subprojects if s is not None
                        ])
                        logger.info(f"Found {len(subproject_options)-1} subproject options")
                    except Exception as e:
                        logger.error(f"Error processing subprojects: {str(e)}")
                
                # Safely get unique dates
                if hasattr(self.diversity_app, 'unique_dates') and self.diversity_app.unique_dates is not None:
                    try:
                        date_options.extend([
                            {'label': str(d), 'value': str(d)} 
                            for d in self.diversity_app.unique_dates if d is not None
                        ])
                        logger.info(f"Found {len(date_options)-1} date options")
                    except Exception as e:
                        logger.error(f"Error processing dates: {str(e)}")

                logger.info("Successfully prepared all dropdown options")
                return sample_options, project_options, subproject_options, date_options

            except Exception as e:
                logger.error(f"Error in update_diversity_dropdowns: {str(e)}")
                logger.error(traceback.format_exc())
                return [], [], [], []

        
        
        # Add MAG-specific callbacks
        @self.app.callback(
            [Output('circos-plot', 'children'),
             Output('mag-sequence-info', 'children'),
             Output('mag-gene-content', 'children')],
            [Input('mag-table-selected-rows', 'data')],
            [State('mag-table-data', 'children')],
            prevent_initial_call=True
        )
        def update_mag_details(selected_rows, table_data_json):
            """Update MAG visualizations when a row is selected."""
            try:
                if not selected_rows or len(selected_rows) == 0:
                    # No MAG selected
                    return dash.no_update, dash.no_update, dash.no_update
                
                # Parse the JSON data from the hidden div
                table_data = json.loads(table_data_json)
                selected_mag = table_data[selected_rows[0]]
                logger.info(f"Selected MAG: {selected_mag.get('name', 'unknown')}")
                
                # Create the visualizations for the selected MAG using mags_app methods
                if hasattr(self, 'mags_app') and self.mags_app:
                    # Generate circos plot HTML
                    circos_plot_html = self.mags_app._create_test_circos_plot(selected_mag.get('id'))
                    circos_plot = html.Div([
                        html.Iframe(
                            srcDoc=circos_plot_html,
                            style={
                                'width': '100%',
                                'height': '800px',
                                'border': 'none'
                            },
                            sandbox='allow-scripts allow-same-origin'
                        )
                    ])
                    
                    # Get basic MAG info for display
                    mag_name = selected_mag.get('name', 'Unknown')
                    taxonomy = selected_mag.get('taxonomy', 'Unknown')
                    completeness = selected_mag.get('completeness', 0)
                    contamination = selected_mag.get('contamination', 0)
                    size_mb = selected_mag.get('size_mb', 0)
                    gc_content = selected_mag.get('gc_content', 0)
                    n_contigs = selected_mag.get('num_contigs', 0)
                    quality = selected_mag.get('quality', 'unknown')
                    
                    # Format values for display
                    completeness_display = f"{completeness * 100:.1f}%" if completeness <= 1 else f"{completeness:.1f}%"
                    contamination_display = f"{contamination * 100:.1f}%" if contamination <= 1 else f"{contamination:.1f}%"
                    gc_content_display = f"{gc_content * 100:.1f}%" if gc_content <= 1 else f"{gc_content:.1f}%"
                    
                    # Create sequence info panel
                    sequence_info = html.Div([
                        html.H4(mag_name, className="mb-3"),
                        dmc.Grid([
                            dmc.GridCol([
                                html.Div([
                                    html.Strong("Taxonomy: "),
                                    html.Span(taxonomy)
                                ], className="mb-2"),
                                html.Div([
                                    html.Strong("Quality: "),
                                    html.Span(quality.capitalize())
                                ], className="mb-2"),
                                html.Div([
                                    html.Strong("Completeness: "),
                                    html.Span(completeness_display)
                                ], className="mb-2"),
                            ], span=6),
                            dmc.GridCol([
                                html.Div([
                                    html.Strong("Contamination: "),
                                    html.Span(contamination_display)
                                ], className="mb-2"),
                                html.Div([
                                    html.Strong("Size: "),
                                    html.Span(f"{size_mb:.2f} Mb")
                                ], className="mb-2"),
                                html.Div([
                                    html.Strong("Contigs: "),
                                    html.Span(f"{n_contigs:,}")
                                ], className="mb-2"),
                                html.Div([
                                    html.Strong("GC Content: "),
                                    html.Span(gc_content_display)
                                ], className="mb-2"),
                            ], span=6),
                        ]),
                    ], style={'padding': '15px', 'backgroundColor': '#f8f9fa', 'borderRadius': '5px', 'marginBottom': '20px'})
                    
                    # Create gene content figure 
                    gene_content = self.mags_app._create_gene_content_plot(selected_mag.get('bin_id'))
                    gene_content_div = html.Div([
                        html.H5("Gene Content Overview"),
                        dcc.Graph(figure=gene_content)
                    ])
                    
                    return circos_plot, sequence_info, gene_content_div
                else:
                    # MAGs app not initialized
                    error_message = html.Div("MAGs app not initialized", style={
                        'textAlign': 'center',
                        'padding': '20px',
                        'backgroundColor': '#f8d7da',
                        'color': '#721c24',
                        'borderRadius': '5px'
                    })
                    return error_message, error_message, error_message
            except Exception as e:
                logger.error(f"Error in update_mag_details: {str(e)}", exc_info=True)
                error_div = html.Div(f"Error: {str(e)}", style={
                    'textAlign': 'center',
                    'padding': '20px',
                    'backgroundColor': '#f8d7da',
                    'color': '#721c24',
                    'borderRadius': '5px'
                })
                return error_div, error_div, error_div

        @self.app.callback(
            Output('sequence-viewer-container', 'children'),
            [Input('sequence-selector', 'value')],
            [State('mag-table', 'selected_rows'),
             State('mag-table', 'data')],
            prevent_initial_call=True
        )
        def update_sequence_viewer(sequence_id, selected_rows, table_data):
            """Update the sequence viewer with the selected sequence."""
            try:
                from users.models import MAG
                
                if not sequence_id:
                    return html.Div(
                        "Please select a sequence from the dropdown to view its details.",
                        style={
                            'padding': '20px',
                            'textAlign': 'center',
                            'backgroundColor': '#f8f9fa',
                            'borderRadius': '5px'
                        }
                    )
                
                # Get the selected MAG
                if not selected_rows or len(selected_rows) == 0 or not table_data or len(table_data) <= selected_rows[0]:
                    return html.Div("Please select a MAG from the table first.")
                
                selected_mag = table_data[selected_rows[0]].get('bin_id')
                if not selected_mag:
                    return html.Div("Selected MAG has no bin_id")
                
                logger.info(f"Getting sequence {sequence_id} from MAG {selected_mag}")
                
                # Get sequence data from database directly
                try:
                    mag = MAG.objects.get(name=selected_mag, user_id=self.user_id)
                    if not mag.fasta_file:
                        return html.Div("No FASTA data available for this MAG.")
                    
                    # Find the specific sequence by parsing the FASTA file
                    sequence_details = None
                    current_header = None
                    current_seq = ""
                    
                    for line in mag.fasta_file.split('\n'):
                        line = line.strip()
                        if not line:
                            continue
                        
                        if line.startswith('>'):
                            if current_header and current_seq:
                                seq_id = current_header.split()[0]
                                if seq_id == sequence_id:
                                    # Calculate GC content
                                    gc_count = current_seq.upper().count('G') + current_seq.upper().count('C')
                                    gc_content = (gc_count / len(current_seq) * 100) if current_seq else 0
                                    
                                    sequence_details = {
                                        'id': seq_id,
                                        'header': current_header,
                                        'sequence': current_seq,
                                        'length': len(current_seq),
                                        'gc_content': gc_content,
                                        'parent_mag': selected_mag,
                                        'taxonomy': mag.taxonomy or 'Unknown'
                                    }
                                    break
                            current_header = line[1:]
                            current_seq = ""
                        else:
                            current_seq += line
                    
                    # Check for the last sequence
                    if not sequence_details and current_header and current_seq:
                        seq_id = current_header.split()[0]
                        if seq_id == sequence_id:
                            gc_count = current_seq.upper().count('G') + current_seq.upper().count('C')
                            gc_content = (gc_count / len(current_seq) * 100) if current_seq else 0
                            
                            sequence_details = {
                                'id': seq_id,
                                'header': current_header,
                                'sequence': current_seq,
                                'length': len(current_seq),
                                'gc_content': gc_content,
                                'parent_mag': selected_mag,
                                'taxonomy': mag.taxonomy or 'Unknown'
                            }
                    
                    if not sequence_details:
                        return html.Div(f"Sequence {sequence_id} not found in MAG {selected_mag}.")
                    
                    # Create formatted sequence view
                    return self._create_sequence_viewer_content(sequence_details)
                    
                except MAG.DoesNotExist:
                    return html.Div(f"MAG {selected_mag} not found in database.")
                except Exception as e:
                    logger.error(f"Error retrieving sequence data: {str(e)}", exc_info=True)
                    return html.Div(f"Error retrieving sequence data: {str(e)}")
                
            except Exception as e:
                logger.error(f"Error processing sequence data: {str(e)}", exc_info=True)
                return html.Div(
                    f"Error processing sequence data: {str(e)}",
                    style={
                        'padding': '20px',
                        'textAlign': 'center',
                        'backgroundColor': '#f8f9fa',
                        'borderRadius': '5px',
                        'color': 'red'
                    }
                )

        
        # MAG table filter callback
        @self.app.callback(
            [Output('mag-table-data', 'children')],
            [Input('quality-filter', 'value')],
            prevent_initial_call=False
        )
        def filter_mag_table_data(quality_value):
            """Filter the MAG table data based on quality selection and update the hidden data store."""
            try:
                logger.info(f"Filtering MAGs by quality: {quality_value}")
                
                # Check if MAGs app is initialized
                if not hasattr(self, 'mags_app') or self.mags_app is None:
                    logger.error("MAGs app not initialized")
                    return [dash.no_update]
                
                # Get all table data
                table_data = self.mags_app._create_table_data()
                
                if not table_data:
                    logger.warning("No MAG data available for filtering")
                    return [dash.no_update]
                
                # Filter by quality if specified
                if quality_value and quality_value.lower() != 'all':
                    filtered_data = [
                        row for row in table_data 
                        if row.get('quality', '').lower() == quality_value.lower()
                    ]
                    logger.info(f"Filtered to {len(filtered_data)} MAGs with quality {quality_value}")
                else:
                    # No filter or 'all' selected
                    filtered_data = table_data
                    logger.info(f"Showing all {len(filtered_data)} MAGs")
                
                # Limit the size of table data to avoid SchemaLengthValidationError
                # Truncate large text fields and limit fasta_file content
                for row in filtered_data:
                    # Remove fasta_file content completely - it's not needed for table display
                    if 'fasta_file' in row:
                        # Just store a placeholder indicating the file exists
                        row['fasta_file'] = 'FASTA_DATA_AVAILABLE' if row['fasta_file'] else ''
                    
                    # Truncate taxonomy to reasonable size
                    if 'taxonomy' in row and isinstance(row['taxonomy'], str) and len(row['taxonomy']) > 200:
                        row['taxonomy'] = row['taxonomy'][:197] + '...'
                
                # Return truncated data as JSON string (filtered_data must be serializable)
                return [json.dumps(filtered_data)]
                
            except Exception as e:
                logger.error(f"Error filtering MAG table: {str(e)}", exc_info=True)
                return [json.dumps([])]

        # Add MAG button selection callback
        @self.app.callback(
            Output('mag-table-selected-rows', 'data'),
            [Input({'type': 'mag-select-btn', 'index': ALL}, 'n_clicks')],
            prevent_initial_call=True
        )
        def update_selected_rows_from_button(n_clicks_list):
            """Update the selected rows when a select button is clicked."""
            try:
                # Get the triggered input
                ctx = dash.callback_context
                if not ctx.triggered:
                    return []
                
                # Get the ID of the button that was clicked
                button_id = ctx.triggered[0]['prop_id'].split('.')[0]
                if not button_id:
                    return []
                    
                # Extract the index from the button ID
                try:
                    button_index = json.loads(button_id)['index']
                    logger.info(f"MAG selection: Button at index {button_index} was clicked")
                    return [button_index]
                except:
                    logger.error(f"Could not parse button ID: {button_id}")
                    return []
                
            except Exception as e:
                logger.error(f"Error in update_selected_rows_from_button: {str(e)}", exc_info=True)
                return []

        # Add callback to update MAG table when data changes
        @self.app.callback(
            Output('mag-table', 'data'),
            [Input('mag-table-data', 'children')],
            prevent_initial_call=False
        )
        def update_mag_table_from_data(table_data_json):
            """Update the MAG table with the filtered data."""
            try:
                if not table_data_json:
                    logger.warning("No table data available")
                    return []
                
                # Parse the data - table_data_json is now a JSON string 
                try:
                    table_data = json.loads(table_data_json)
                    logger.info(f"Successfully parsed table data with {len(table_data)} rows")
                except Exception as e:
                    logger.error(f"Error parsing table data JSON: {str(e)}")
                    return []
                
                if not table_data:
                    logger.warning("Empty table data after parsing JSON")
                    return []
                
                # Return the data directly - DataTable expects a list of dicts
                return table_data
                
            except Exception as e:
                logger.error(f"Error updating MAG table: {str(e)}", exc_info=True)
                return []

        @self.app.callback(
            [Output('alpha_diversities_plot1', 'figure'),
             Output('alpha_diversities_plot2', 'figure'),
             Output('diversity_time_series', 'figure'),
             Output('beta_diversity_heatmap', 'figure'),
             Output('pcoa_plot_container', 'figure'),
             Output('rarefaction_plot', 'figure')],
             
            [Input('sample-select-value-diversity', 'value'),
             Input('diversity-metric-dropdown', 'value'),
             Input('project-dropdown-diversity', 'value'),
             Input('subproject-dropdown-diversity', 'value'),
             Input('date-dropdown-diversity', 'value'),
             Input('color-by-dropdown', 'value'),
             Input('taxonomic-rank-dropdown', 'value'),
             Input('kmeans-input', 'value'),
             Input('arrow-toggle', 'value'),
             Input('pcoa-3d-toggle', 'value')],
            prevent_initial_call=False
        )
        def update_diversity_plots(selected_samples, selected_metrics, selected_project,
                                 selected_subproject, selected_date, color_by,
                                 tax_rank, kmeans_clusters, show_arrows, show_3d):
            try:
                logger.info("\n=== Starting Diversity Plot Updates ===")
                logger.info(f"Selected samples: {selected_samples}")
                logger.info(f"Selected metrics: {selected_metrics}")
                
                # Create empty figure for error cases
                def create_empty_figure(message="No data available"):
                    fig = go.Figure()
                    fig.add_annotation(
                        text=message,
                        xref="paper",
                        yref="paper",
                        x=0.5,
                        y=0.5,
                        showarrow=False
                    )
                    fig.update_layout(
                        plot_bgcolor='white',
                        paper_bgcolor='white'
                    )
                    return fig
                
                # Check if diversity app is initialized
                if not hasattr(self, 'diversity_app') or self.diversity_app is None:
                    error_msg = "Diversity app not initialized"
                    logger.error(error_msg)
                    empty_fig = create_empty_figure(error_msg)
                    return [empty_fig] * 6 + [[], []]

                # Update taxonomic rank with error handling
                try:
                    if tax_rank:
                        self.diversity_app.update_taxonomic_rank(tax_rank)
                    else:
                        self.diversity_app.update_taxonomic_rank('species')
                except Exception as e:
                    logger.error(f"Error updating taxonomic rank: {str(e)}")
                    logger.error(traceback.format_exc())
                
                # Get all samples if none selected
                if not selected_samples:
                    if hasattr(self.diversity_app, 'unique_sample_ids') and self.diversity_app.unique_sample_ids is not None:
                        selected_samples = list(self.diversity_app.unique_sample_ids)
                    else:
                        selected_samples = []
                elif isinstance(selected_samples, str):
                    selected_samples = [selected_samples]
                
                # Make sure selected_samples is a list
                if not isinstance(selected_samples, list):
                    selected_samples = list(selected_samples) if hasattr(selected_samples, '__iter__') else [selected_samples]
                    
                # Filter samples based on project/subproject/date
                filtered_samples = selected_samples.copy()
                
                if selected_project and selected_project != 'ALL' and hasattr(self.diversity_app, 'sample_to_project_dict'):
                    filtered_samples = [s for s in filtered_samples 
                                     if self.diversity_app.sample_to_project_dict.get(s) == selected_project]
                
                if selected_subproject and selected_subproject != 'ALL' and hasattr(self.diversity_app, 'sample_to_subproject_dict'):
                    filtered_samples = [s for s in filtered_samples 
                                     if self.diversity_app.sample_to_subproject_dict.get(s) == selected_subproject]
                
                if selected_date and selected_date != 'ALL' and hasattr(self.diversity_app, 'sample_to_date_dict'):
                    filtered_samples = [s for s in filtered_samples 
                                     if str(self.diversity_app.sample_to_date_dict.get(s)) == str(selected_date)]

                if not filtered_samples:
                    logger.warning("No samples match the selected filters")
                    empty_fig = create_empty_figure("No samples match the selected filters")
                    return [empty_fig] * 6 + [[], []]

                # Set default values for parameters
                if not selected_metrics:
                    selected_metrics = ['shannon', 'simpson']
                elif isinstance(selected_metrics, str):
                    selected_metrics = [selected_metrics]
                    
                # Convert types and set defaults
                try:
                    kmeans_clusters = int(kmeans_clusters) if kmeans_clusters else 3
                except (ValueError, TypeError):
                    kmeans_clusters = 3
                    
                show_arrows = bool(show_arrows)
                show_3d = bool(show_3d)
                
                # Safely get color_by parameter
                if not color_by:
                    color_by = 'project'

                # Calculate summary statistics
                try:
                    stats = self.diversity_app.calculate_summary_statistics()
                    if not stats:
                        stats = {}
                except Exception as e:
                    logger.error(f"Error calculating statistics: {str(e)}")
                    stats = {}
                
                # Create summary text components
                diversity_summary = dmc.List([
                    dmc.ListItem([
                        dmc.Text("Shannon Diversity (H′):"),
                        dmc.Text(
                            f"Mean: {stats.get('shannon', {}).get('mean', 0):.2f}" if stats.get('shannon', {}).get('mean') is not None else "No data",
                            
                        )
                    ]),
                    dmc.ListItem([
                        dmc.Text("Simpson Diversity (D):"),
                        dmc.Text(
                            f"Mean: {stats.get('simpson', {}).get('mean', 0):.2f}" if stats.get('simpson', {}).get('mean') is not None else "No data",
                            
                        )
                    ]),
                    dmc.ListItem([
                        dmc.Text("Species Richness:"),
                        dmc.Text(
                            f"Total: {stats.get('total_species', 0)}, Mean per sample: {stats.get('mean_species_per_sample', 0):.1f}" if stats.get('total_species') is not None else "No data",
                            
                        )
                    ])
                ])
                
                sample_summary = dmc.List([
                    dmc.ListItem([
                        dmc.Text("Total Samples:"),
                        dmc.Text(
                            str(stats.get('total_samples', 0)) if stats.get('total_samples') is not None else "No data",
                            
                        )
                    ]),
                    dmc.ListItem([
                        dmc.Text("Projects:"),
                        dmc.Text(
                            str(stats.get('total_projects', 0)) if stats.get('total_projects') is not None else "No data",
                            
                        )
                    ])
                ])

                # Use a direct approach for alpha diversity plots if the main function fails
                try:
                    # Try to get plots from diversity app
                    plots = self.diversity_app.update_diversity_plots(
                        selected_samples=filtered_samples,
                        diversity_metric=selected_metrics,
                        color_by=color_by,
                        k_value=kmeans_clusters,
                        project=selected_project,
                        subproject=selected_subproject,
                        date=selected_date,
                        is_3d=show_3d,
                        show_arrows=show_arrows,
                        tax_rank=tax_rank
                    )
                    
                    # Verify that we got valid plots back
                    if not plots or len(plots) < 6:
                        logger.warning(f"Invalid or missing plots returned from diversity_app: Got {len(plots) if plots else 0} plots, expected 6")
                        plots = [create_empty_figure("Error generating plot")] * 6
                    
                    # Validate each plot
                    for i, plot in enumerate(plots):
                        if not isinstance(plot, go.Figure):
                            logger.warning(f"Plot {i} is not a valid figure, replacing with empty figure")
                            plots[i] = create_empty_figure(f"Error with plot {i}")
                    
                    return plots
                    
                except Exception as e:
                    logger.error(f"Error in update_diversity_plots using app method: {str(e)}")
                    logger.error(traceback.format_exc())
                    
                    # Fallback to creating basic plots
                    try:
                        # Create alpha diversity plots directly
                        if hasattr(self.diversity_app, "create_alpha_diversity_plots"):
                            alpha_fig1, alpha_fig2 = self.diversity_app.create_alpha_diversity_plots(filtered_samples)
                        else:
                            alpha_fig1 = create_empty_figure("Alpha diversity plot unavailable")
                            alpha_fig2 = create_empty_figure("Alpha diversity plot unavailable")
                            
                        # Create time series plot directly
                        if hasattr(self.diversity_app, "create_diversity_time_series"):
                            time_series_fig = self.diversity_app.create_diversity_time_series(filtered_samples)
                        else:
                            time_series_fig = create_empty_figure("Time series plot unavailable")
                            
                        # Create placeholder plots for other visualizations
                        beta_fig = create_empty_figure("Beta diversity plot unavailable")
                        pcoa_fig = create_empty_figure("PCoA plot unavailable")
                        rare_fig = create_empty_figure("Rarefaction plot unavailable")
                        
                        return [alpha_fig1, alpha_fig2, time_series_fig, beta_fig, pcoa_fig, rare_fig, diversity_summary, sample_summary]
                    except Exception as e2:
                        logger.error(f"Error in fallback plot creation: {str(e2)}")
                        logger.error(traceback.format_exc())
                        
                        # Ultimate fallback - empty figures
                        empty_figs = [create_empty_figure(f"Error: {str(e)}")] * 6
                        return empty_figs + [diversity_summary, sample_summary]

            except Exception as e:
                logger.error(f"Unexpected error in update_diversity_plots: {str(e)}")
                logger.error(traceback.format_exc())
                
                empty_fig = go.Figure()
                empty_fig.add_annotation(
                    text=f"Error: {str(e)}",
                    xref="paper",
                    yref="paper",
                    x=0.5,
                    y=0.5,
                    showarrow=False
                )
                
                return [empty_fig] * 6 + [[], []]

        # Other callback registrations can go here if needed
        pass
        
        # Add callback for downloading diversity data as CSV
        @self.app.callback(
            Output('download-diversity-csv', 'data'),
            Input('download-diversity-csv-btn', 'n_clicks'),
            prevent_initial_call=True
        )
        def download_diversity_csv(n_clicks):
            """Generate and download CSV file with diversity metrics and statistics."""
            try:
                if not n_clicks:
                    raise PreventUpdate
                
                if not hasattr(self, 'diversity_app') or self.diversity_app is None:
                    logger.error("Diversity app not initialized")
                    raise PreventUpdate
                
                # Get the diversity data
                logger.info("Preparing diversity data for CSV export")
                data_dict = self.diversity_app.prepare_diversity_csv_data()
                
                if not data_dict or not isinstance(data_dict, dict):
                    logger.error("No valid diversity data available for export")
                    raise PreventUpdate
                
                # Extract the DataFrames
                diversity_df = data_dict.get('diversity_data', pd.DataFrame())
                stats_df = data_dict.get('summary_stats', pd.DataFrame())
                
                if diversity_df.empty:
                    logger.warning("No diversity data available for export")
                    raise PreventUpdate
                
                # Create a CSV string with two sections
                import io
                
                buffer = io.StringIO()
                
                # Write the diversity data section
                buffer.write("# SAMPLE DIVERSITY METRICS\n")
                diversity_df.to_csv(buffer, index=False)
                
                # Add a separator between sections
                buffer.write("\n\n# SUMMARY STATISTICS\n")
                
                # Write the summary statistics if available
                if not stats_df.empty:
                    stats_df.to_csv(buffer, index=False)
                else:
                    buffer.write("No summary statistics available\n")
                
                # Get current date for filename
                from datetime import datetime
                current_date = datetime.now().strftime("%Y-%m-%d")
                
                # Return the CSV data
                return dcc.send_string(
                    buffer.getvalue(),
                    filename=f"diversity_metrics_{current_date}.csv"
                )
                
            except Exception as e:
                logger.error(f"Error generating diversity CSV: {str(e)}")
                logger.error(traceback.format_exc())
                raise PreventUpdate
                
        # Other callbacks go here
        
            # Add clipboard callbacks
        @self.app.callback(
            [Output('copy-alert', 'style'),
            Output('copy-alert', 'children')],
            [Input('copy-sequence-btn', 'n_clicks')],
            [State('sequence-data', 'children')],
            prevent_initial_call=True
        )
        def copy_sequence_to_clipboard(seq_clicks, sequence):
            """Handle sequence copying to clipboard."""
            try:
                if not dash.callback_context.triggered:
                    raise PreventUpdate
                    
                button_id = dash.callback_context.triggered[0]['prop_id'].split('.')[0]
                
                if button_id == 'copy-sequence-btn' and seq_clicks:
                    return {'display': 'block', 'marginTop': '10px'}, "Sequence copied to clipboard!"
                    
                return {'display': 'none'}, ""
            except Exception as e:
                logger.error(f"Error copying sequence to clipboard: {str(e)}")
                logger.error(traceback.format_exc())
                return {'display': 'block', 'marginTop': '10px', 'color': 'red'}, f"Error copying sequence: {str(e)}"

        @self.app.callback(
            Output('copy-alert', 'style', allow_duplicate=True),
            [Input('copy-alert', 'children')],
            prevent_initial_call=True
        )
        def hide_copy_alert(alert_text):
            """Hide the copy alert after a delay."""
            try:
                if not alert_text:
                    raise PreventUpdate
                    
                # Hide after 3 seconds
                time.sleep(3)
                return {'display': 'none'}
            except Exception as e:
                logger.error(f"Error hiding copy alert: {str(e)}")
                logger.error(traceback.format_exc())
                return {'display': 'none'}
    
    def _get_app_icon(self, app_name: str) -> str:
        """Get the icon name for a given app."""
        # For custom image apps, use a special prefix to indicate custom icons
        custom_icons = {
            'Taxonomy': 'custom:microorganism.png',
            'MAGs': 'custom:circ-chromosome.png'
        }
        
        if app_name in custom_icons:
            return custom_icons[app_name]
            
        # Standard iconify icons
        icon_map = {
            'Horizon': 'mdi:chart-timeline',
            'Diversity': 'tabler:chart-scatter-3d',
            'Correlations': 'carbon:qq-plot',
            'QC': 'material-symbols:fact-check-outline'
        }
        return icon_map.get(app_name, 'mdi:application')

    def _get_cache_key(self, prefix, *args, **kwargs):
        """Generate a cache key from arguments"""
        key_dict = {
            'args': args,
            'kwargs': kwargs,
            'user_id': self.user_id
        }
        key_str = json.dumps(key_dict, sort_keys=True)
        return f"{prefix}:{hashlib.md5(key_str.encode()).hexdigest()}"

    @lru_cache(maxsize=128)
    def _get_cached_data(self):
        """Get cached data with performance optimizations"""
        start = time.time()
        try:
            records = NanoporeRecord.objects.filter(user_id=self.user_id).values()
            df = pd.DataFrame.from_records(records)
            duration = time.time() - start
            logger.info(f"Data fetching and processing took {duration:.2f} seconds")
            return df
        except Exception as e:
            logger.error(f"Error fetching data: {str(e)}")
            return pd.DataFrame()

    def _optimize_dataframe(self, df):
        """Optimize DataFrame memory usage"""
        start = time.time()
        
        # Convert date columns to datetime
        date_cols = [col for col in df.columns if 'date' in col.lower()]
        for col in date_cols:
            df[col] = pd.to_datetime(df[col], errors='ignore')
        
        # Optimize numeric columns
        for col in df.select_dtypes(include=['float64']).columns:
            df[col] = pd.to_numeric(df[col], downcast='float')
        
        for col in df.select_dtypes(include=['int64']).columns:
            df[col] = pd.to_numeric(df[col], downcast='integer')
        
        # Categorize string columns with low cardinality
        for col in df.select_dtypes(include=['object']).columns:
            unique_ratio = df[col].nunique() / len(df)
            if unique_ratio < 0.5:  # If less than 50% unique values
                df[col] = df[col].astype('category')
        
        duration = time.time() - start
        memory_usage = df.memory_usage(deep=True).sum() / 1024**2
        logger.info(f"DataFrame optimization took {duration:.2f} seconds. Memory usage: {memory_usage:.2f} MB")
        
        return df

    def _precalculate_aggregations(self):
        """Pre-calculate common aggregations to improve performance"""
        if len(self.df) == 0:
            return
        
        # Calculate common aggregations
        self.df_sorted = self.df.sort_values(['project_id', 'subproject', 'sample_id'])
        self.unique_samples = self.df['sample_id'].unique()
        self.unique_projects = self.df['project_id'].unique()
        
        # Ensure dates are properly converted to datetime before finding min/max
        try:
            # Convert to datetime if not already
            if 'date' in self.df.columns:
                if not pd.api.types.is_datetime64_dtype(self.df['date']):
                    date_series = pd.to_datetime(self.df['date'], errors='coerce')
                else:
                    date_series = self.df['date']
                    
                # Get unique dates as regular datetime objects, not categorical
                self.unique_dates = date_series.dropna().unique()
                
                # Get min and max dates safely
                min_date = None
                max_date = None
                if len(self.unique_dates) > 0:
                    min_date = self.unique_dates.min()
                    max_date = self.unique_dates.max()
            else:
                self.unique_dates = np.array([])
                min_date = None
                max_date = None
        except Exception as e:
            logger.error(f"Error processing dates: {str(e)}")
            self.unique_dates = np.array([])
            min_date = None
            max_date = None
        
        # Store in cache
        self._cache['aggregations'] = {
            'sample_counts': self.df.groupby('sample_id', observed=True).size().to_dict(),
            'project_stats': self.df.groupby('project_id', observed=True).agg({
                'sample_id': 'nunique',
                'taxonomy': 'count'
            }).to_dict(),
            'date_range': {
                'min': min_date,
                'max': max_date
            }
        }

  
        
        # Add callback to update gene content plot
    
    def _get_app_icon(self, app_name: str) -> str:
        """Get the icon name for a given app."""
        # For custom image apps, use a special prefix to indicate custom icons
        custom_icons = {
            'Taxonomy': 'custom:microorganism.png',
            'MAGs': 'custom:circ-chromosome.png'
        }
        
        if app_name in custom_icons:
            return custom_icons[app_name]
            
        # Standard iconify icons
        icon_map = {
            'Horizon': 'mdi:chart-timeline',
            'Diversity': 'tabler:chart-scatter-3d',
            'Correlations': 'carbon:qq-plot',
            'QC': 'material-symbols:fact-check-outline'
        }
        return icon_map.get(app_name, 'mdi:application')

    def _get_cache_key(self, prefix, *args, **kwargs):
        """Generate a cache key from arguments"""
        key_dict = {
            'args': args,
            'kwargs': kwargs,
            'user_id': self.user_id
        }
        key_str = json.dumps(key_dict, sort_keys=True)
        return f"{prefix}:{hashlib.md5(key_str.encode()).hexdigest()}"

    @lru_cache(maxsize=128)
    def _get_cached_data(self):
        """Get cached data with performance optimizations"""
        start = time.time()
        try:
            records = NanoporeRecord.objects.filter(user_id=self.user_id).values()
            df = pd.DataFrame.from_records(records)
            duration = time.time() - start
            logger.info(f"Data fetching and processing took {duration:.2f} seconds")
            return df
        except Exception as e:
            logger.error(f"Error fetching data: {str(e)}")
            return pd.DataFrame()

    def _optimize_dataframe(self, df):
        """Optimize DataFrame memory usage"""
        start = time.time()
        
        # Convert date columns to datetime
        date_cols = [col for col in df.columns if 'date' in col.lower()]
        for col in date_cols:
            df[col] = pd.to_datetime(df[col], errors='ignore')
        
        # Optimize numeric columns
        for col in df.select_dtypes(include=['float64']).columns:
            df[col] = pd.to_numeric(df[col], downcast='float')
        
        for col in df.select_dtypes(include=['int64']).columns:
            df[col] = pd.to_numeric(df[col], downcast='integer')
        
        # Categorize string columns with low cardinality
        for col in df.select_dtypes(include=['object']).columns:
            unique_ratio = df[col].nunique() / len(df)
            if unique_ratio < 0.5:  # If less than 50% unique values
                df[col] = df[col].astype('category')
        
        duration = time.time() - start
        memory_usage = df.memory_usage(deep=True).sum() / 1024**2
        logger.info(f"DataFrame optimization took {duration:.2f} seconds. Memory usage: {memory_usage:.2f} MB")
        
        return df

    

    def _get_sequence_details(self, sequence_id):
        """Get details for a specific sequence."""
        try:
            if not sequence_id or not hasattr(self, 'mags_app'):
                return None
                
            # Get MAG data with fasta content
            mag_data = self.mags_app.get_mag_data(include_fasta=True)
            
            if not mag_data:
                logger.warning("No MAG data available")
                return None
            
            # Check all MAGs for the sequence
            for mag_id, mag_info in mag_data.items():
                if not isinstance(mag_info, dict):
                    continue
                    
                fasta_data = mag_info.get('fasta', [])
                
                # Parse fasta file if needed
                if not fasta_data and 'fasta_file' in mag_info:
                    fasta_data = self.mags_app.parse_fasta_for_viewer(mag_info['fasta_file'])
                    
                # Look for the sequence in this MAG
                for seq in fasta_data:
                    if seq.get('id') == sequence_id:
                        # Calculate GC content if not present
                        if 'gc_content' not in seq:
                            sequence = seq.get('sequence', '')
                            gc_count = sequence.upper().count('G') + sequence.upper().count('C')
                            total_len = len(sequence)
                            if total_len > 0:
                                seq['gc_content'] = (gc_count / total_len) * 100
                            else:
                                seq['gc_content'] = 0
                        
                        # Add parent MAG info for context
                        seq['parent_mag'] = mag_info.get('bin_id', mag_id)
                        seq['taxonomy'] = mag_info.get('taxonomy', 'Unknown')
                        
                        # Try to get annotation data for this sequence
                        try:
                            from users.models import MAG
                            mag_obj = MAG.objects.get(name=mag_info.get('bin_id', mag_id), user_id=self.user_id)
                            if hasattr(mag_obj, 'gff_file') and mag_obj.gff_file:
                                seq['has_annotations'] = True
                        except Exception as e:
                            logger.error(f"Error checking for annotations: {str(e)}")
                            seq['has_annotations'] = False
                            
                        return seq
            
            # If sequence not found, try as a fallback mechanism
            if sequence_id in ['complete_sequence', 'sequence_1', 'emergency_sequence']:
                logger.info(f"Using fallback mechanism for sequence {sequence_id}")
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
                            'gc_content': self._calculate_gc_content_for_sequence(cleaned_sequence),
                            'parent_mag': mag.get('bin_id', mag_id),
                            'taxonomy': mag.get('taxonomy', 'Unknown')
                        }
            
            return None
            
        except Exception as e:
            logger.error(f"Error getting sequence details: {str(e)}", exc_info=True)
            return None
            
    def _calculate_gc_content_for_sequence(self, sequence):
        """Calculate GC content for a DNA sequence as a percentage."""
        if not sequence:
            return 0.0
            
        sequence = sequence.upper()
        gc_count = sequence.count('G') + sequence.count('C')
        total_bases = len(sequence)
        
        return ((gc_count / total_bases) * 100)

        # Removed NGCircos clientside callback since we now use pyCirclize for visualization
        # PyCirclize doesn't require JavaScript dependencies as it generates static images server-side

    def _create_sequence_viewer_content(self, sequence_details):
        """Create content for the sequence viewer with enhanced features."""
        logger = logging.getLogger(__name__)
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
