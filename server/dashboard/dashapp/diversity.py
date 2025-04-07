"""
The Diversity app provides comprehensive analysis of microbial diversity in your samples.

Features:
- Alpha Diversity: Measures diversity within individual samples using multiple metrics
  (Shannon, Simpson, and Observed Features)
- Beta Diversity: Compares diversity between samples using distance metrics and heatmaps
- PCoA Analysis: Principal Coordinate Analysis to visualize sample relationships
- Rarefaction Curves: Evaluate sampling depth and species richness
- Interactive Filtering: Filter by project, subproject, date, and taxonomic rank
- Dynamic Visualization: Customize plots with different coloring schemes and clustering options

This app helps you understand:
- How diverse are your individual samples?
- How do samples compare to each other?
- Are there patterns in microbial composition across projects or time?
- Is your sampling depth sufficient?
"""

from dash import dcc, html, callback_context, Input, Output, State, no_update, dash_table
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import skbio
from django_plotly_dash import DjangoDash
from skbio import DistanceMatrix
from skbio.stats import ordination
from users.models import NanoporeRecord, SequencingStatistics
from dash_iconify import DashIconify
from sklearn.cluster import KMeans
import traceback
import time
import hashlib
import logging
from functools import lru_cache
import json
from typing import Dict, List, Tuple, Any
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy import stats

logger = logging.getLogger(__name__)
# Suppress numpy warnings
# warnings.filterwarnings('ignore', category=RuntimeWarning)
# warnings.filterwarnings('ignore', category=FutureWarning)
# warnings.filterwarnings('ignore', category=DeprecationWarning)


def distance_matrix_to_dataframe(distance_matrix):
    """Convert distance matrix to pandas DataFrame with proper type handling"""
    matrix_df = pd.DataFrame(distance_matrix.data.astype(float),
                            index=distance_matrix.ids,
                            columns=distance_matrix.ids)
    return matrix_df


def calculate_normalized_counts(self):
    counts_df = pd.DataFrame.from_records(NanoporeRecord.objects.filter(user_id=self.user_id).values())
    stats_df = pd.DataFrame.from_records(SequencingStatistics.objects.filter(user_id=self.user_id).values())

    merged_df = pd.merge(counts_df, stats_df, left_on='sample_id', right_on="sample_name")
    merged_df['normalized_count'] = (merged_df['count'] / merged_df['total_bases_sequenced']) * 10000000
    normalized_counts_df = merged_df[['sample_id', 'taxonomy', 'abundance', 'count', 'normalized_count']]
    return normalized_counts_df


def rarefy_sample(sample, depth):
    if sum(sample) >= depth:
        return np.random.choice(np.where(sample > 0)[0], size=depth, replace=False, p=sample / sum(sample))
    else:
        raise ValueError("Sample size is less than the specified rarefaction depth.")
    return sample


class Diversity:
    def __init__(self, user_id):
        
        
        
        app_name = 'diversity'
        
        # Initialize the Dash app with appropriate external stylesheets
        external_stylesheets = [dbc.themes.BOOTSTRAP]
        
        try:
            # Initialize the Dash app
            self.app = DjangoDash(
                app_name,
                add_bootstrap_links=False,
                external_stylesheets=external_stylesheets,
                serve_locally=True  # Use local resources for better stability
            )
            
            # Initialize the app layout
            self.app.layout = self.create_layout()
            
            logger.info("MAGs app initialized successfully")
        except Exception as e:
            logger.error(f"Error initializing MAGs app: {str(e)}", exc_info=True)
            # Create a minimal app to avoid breaking the whole dashboard
            self.app = DjangoDash(app_name, suppress_callback_exceptions=True)
            self.app.layout = html.Div("Error loading MAGs app")
        
        # Note: All callbacks for this app are now defined in index.py
        # This ensures no duplicate callbacks and better centralization
        
        self.user_id = user_id
        self._cache = {}
        self._cache_timeout = 3600  # 1 hour cache timeout
        self._last_data_update = None
        self._performance_logs = []
        self.get_data()
        
        self.records = NanoporeRecord.objects.filter(user_id=self.user_id)
        if not self.records.exists():
            print(f"No records found for user_id: {self.user_id}")
            
        self.plot_ids = ['alpha_diversities_plot1', 'alpha_diversities_plot2', 'beta_diversity_heatmap',
                         'pcoa_plot_container', 'rarefaction_plot', 'diversity_time_series']
        
        self.df = pd.DataFrame()
        if self.records.exists():
            self.df = self.get_data()
            
            if not self.df.empty:
                self._init_data()
                
        self._init_layout()

    def _get_cache_key(self, prefix: str, *args, **kwargs) -> str:
        """Generate a unique cache key."""
        key_dict = {
            'args': args,
            'kwargs': kwargs,
            'user_id': self.user_id,
            'last_update': self._last_data_update
        }
        key_str = json.dumps(key_dict, sort_keys=True)
        return f"{prefix}:{hashlib.md5(key_str.encode()).hexdigest()}"
        
    def _get_cached_result(self, key: str) -> Tuple[bool, Any]:
        """Get a cached result if it exists and is valid."""
        if key in self._cache:
            cached_time, result = self._cache[key]
            if time.time() - cached_time < self._cache_timeout:
                return True, result
        return False, None
        
    def _set_cached_result(self, key: str, result: Any):
        """Store a result in the cache."""
        self._cache[key] = (time.time(), result)
        
        if len(self._cache) > 100:  
            oldest_key = min(self._cache.keys(), key=lambda k: self._cache[k][0])
            del self._cache[oldest_key]
            
    def _init_data(self):
        """Initialize all data-related attributes"""
        try:
            print("Initializing data attributes...")
            self.unique_sample_ids = self.records.values('sample_id').distinct()
            self.unique_sample_ids = [item['sample_id'] for item in self.unique_sample_ids]
            self.unique_samples = self.records.values('sample_id').distinct()
            self.unique_projects_ids = self.records.values('project_id').distinct()
            self.unique_projects_ids = [item['project_id'] for item in self.unique_projects_ids]

            self.unique_subprojects = self.records.values('subproject').distinct()
            self.unique_subprojects = [item['subproject'] for item in self.unique_subprojects]
            
            # Get unique taxonomic ranks
            self.unique_species = self.df['taxonomy'].unique()
            self.unique_genera = self.df['tax_genus'].unique()
            self.unique_families = self.df['tax_family'].unique()
            self.unique_classes = self.df['tax_class'].unique()
            self.unique_orders = self.df['tax_order'].unique()
            self.unique_phyla = self.df['tax_phylum'].unique()

            print("Creating diversity matrices...")
            # Store the full dataframe for different taxonomic ranks
            self.df_by_rank = {
                'species': self.df.pivot_table(
                    index='sample_id',
                    columns='taxonomy',
                    values='abundance',
                    fill_value=0
                ),
                'genus': self.df.pivot_table(
                    index='sample_id',
                    columns='tax_genus',
                    values='abundance',
                    aggfunc='sum',
                    fill_value=0
                ),
                'family': self.df.pivot_table(
                    index='sample_id',
                    columns='tax_family',
                    values='abundance',
                    aggfunc='sum',
                    fill_value=0
                ),
                'class': self.df.pivot_table(
                    index='sample_id',
                    columns='tax_class',
                    values='abundance',
                    aggfunc='sum',
                    fill_value=0
                ),
                'order': self.df.pivot_table(
                    index='sample_id',
                    columns='tax_order',
                    values='abundance',
                    aggfunc='sum',
                    fill_value=0
                ),
                 'phylum': self.df.pivot_table(
                    index='sample_id',
                    columns='tax_phylum',
                    values='abundance',
                    aggfunc='sum',
                    fill_value=0
                )
            }
            
            self.df_full_for_diversity = self.df_by_rank['species']  # Default to species level
            
            print("Calculating alpha diversity...")
            self.calculate_alpha_diversity(use_normalized_counts=False)
            
            print("Calculating beta diversity...")
            try:
                self.beta_diversity_matrix = self.calculate_beta_diversity()
            except Exception as e:
                print(f"Error calculating beta diversity: {e}")
                print(traceback.format_exc())

            print("Getting unique dates...")
            self.unique_dates = self.records.values('date').distinct()
            self.unique_dates = [item['date'] for item in self.unique_dates]
            
            # Create mapping dictionaries for sample metadata
            print("Creating sample metadata mappings...")
            sample_metadata = self.records.values('sample_id', 'project_id', 'subproject', 'date')
            self.sample_to_project_dict = {item['sample_id']: item['project_id'] for item in sample_metadata}
            self.sample_to_subproject_dict = {item['sample_id']: item['subproject'] for item in sample_metadata}
            self.sample_to_date_dict = {item['sample_id']: item['date'] for item in sample_metadata}
            
            print("Data initialization complete")
            
        except Exception as e:
            print(f"Error in _init_data: {str(e)}")
            print(traceback.format_exc())
            
    def _benchmark(func):
        """Decorator to benchmark function execution time"""
        def wrapper(self, *args, **kwargs):
            start_time = time.time()
            result = func(self, *args, **kwargs)
            
            # Log performance data
            perf_data = {
                'function': func.__name__,
                'execution_time': time.time() - start_time,
                'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                'args': str(args),
                'kwargs': str(kwargs)
            }
            self._performance_logs.append(perf_data)
            
            # Keep only last 100 logs
            if len(self._performance_logs) > 100:
                self._performance_logs = self._performance_logs[-100:]
            
            print(f"Performance: {func.__name__} took {time.time() - start_time:.3f} seconds")
            return result
        return wrapper

    def _clear_cache(self):
        """Clear the cache when data is updated"""
        self._cache = {}
        self._last_data_update = time.time()

    @_benchmark
    def create_alpha_diversity_plots(self, selected_samples):
        """Create alpha diversity plots with improved hover information."""
        if not selected_samples:
            # Return empty plots with messages
            empty_line_fig = go.Figure()
            empty_line_fig.add_annotation(
                text="No samples selected",
                xref="paper", yref="paper",
                x=0.5, y=0.5,
                showarrow=False
            )
            empty_line_fig.update_layout(
                title='Alpha Diversity Measures Over Samples',
                plot_bgcolor='white',
                paper_bgcolor='white'
            )
            
            empty_box_fig = go.Figure()
            empty_box_fig.add_annotation(
                text="No samples selected",
                xref="paper", yref="paper",
                x=0.5, y=0.5,
                showarrow=False
            )
            empty_box_fig.update_layout(
                title='Distribution of Alpha Diversity Measures',
                plot_bgcolor='white',
                paper_bgcolor='white'
            )
            
            return empty_line_fig, empty_box_fig

        try:
            # Calculate diversity metrics for selected samples
            diversity_df = self.calculate_alpha_diversity()
            if diversity_df is None or diversity_df.empty:
                print("No diversity data available")
                empty_fig1 = go.Figure()
                empty_fig1.add_annotation(
                    text="No diversity data available",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5,
                    showarrow=False
                )
                empty_fig1.update_layout(
                    title='Alpha Diversity Measures Over Samples',
                    plot_bgcolor='white',
                    paper_bgcolor='white'
                )
                
                empty_fig2 = go.Figure()
                empty_fig2.add_annotation(
                    text="No diversity data available",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5,
                    showarrow=False
                )
                empty_fig2.update_layout(
                    title='Distribution of Alpha Diversity Measures',
                    plot_bgcolor='white',
                    paper_bgcolor='white'
                )
                
                return empty_fig1, empty_fig2
            
            # Filter for selected samples - handle case where some samples aren't in the data
            available_samples = [s for s in selected_samples if s in diversity_df.index]
            
            if not available_samples:
                print("None of the selected samples have diversity data")
                empty_fig1 = go.Figure()
                empty_fig1.add_annotation(
                    text="Selected samples have no diversity data",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5,
                    showarrow=False
                )
                empty_fig1.update_layout(
                    title='Alpha Diversity Measures Over Samples',
                    plot_bgcolor='white',
                    paper_bgcolor='white'
                )
                
                empty_fig2 = go.Figure()
                empty_fig2.add_annotation(
                    text="Selected samples have no diversity data",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5,
                    showarrow=False
                )
                empty_fig2.update_layout(
                    title='Distribution of Alpha Diversity Measures',
                    plot_bgcolor='white',
                    paper_bgcolor='white'
                )
                
                return empty_fig1, empty_fig2
            
            diversity_df = diversity_df.loc[available_samples]
            
            # Create a simpler line plot
            line_fig = go.Figure()
            colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
            
            for i, metric in enumerate(diversity_df.columns):
                # Use modulo operation to cycle through colors if we have more metrics than colors
                color_idx = i % len(colors)
                
                line_fig.add_trace(
                    go.Scatter(
                        x=diversity_df.index,
                        y=diversity_df[metric],
                        name=metric.capitalize(),
                        mode='lines+markers',
                        marker=dict(size=8),
                        line=dict(color=colors[color_idx], width=2),
                        hovertemplate='<b>%{x}</b><br>' + f'{metric.capitalize()}: ' + '%{y:.3f}<extra></extra>'
                    )
                )
            
            line_fig.update_layout(
                title='Alpha Diversity Measures Over Samples',
                xaxis_title='Sample',
                yaxis_title='Diversity Value',
                showlegend=True,
                plot_bgcolor='white',
                paper_bgcolor='white',
                height=400,
                xaxis=dict(tickangle=45)
            )
            
            # Create a simpler box plot
            box_fig = go.Figure()
            
            for i, metric in enumerate(diversity_df.columns):
                color_idx = i % len(colors)
                
                box_fig.add_trace(
                    go.Box(
                        y=diversity_df[metric],
                        name=metric.capitalize(),
                        boxpoints='all',
                        jitter=0.3,
                        pointpos=-1.8,
                        marker_color=colors[color_idx],
                        hovertemplate='<b>' + f'{metric.capitalize()}' + '</b><br>Value: %{y:.3f}<extra></extra>'
                    )
                )
            
            if len(diversity_df.columns) == 0:
                box_fig.add_annotation(
                    text="No valid metrics to display",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5,
                    showarrow=False
                )
            
            box_fig.update_layout(
                title='Distribution of Alpha Diversity Measures',
                showlegend=False,
                plot_bgcolor='white',
                paper_bgcolor='white',
                height=400,
                yaxis_title='Diversity Value',
                xaxis_title='Diversity Measure'
            )

            return line_fig, box_fig
            
        except Exception as e:
            print(f"Error in create_alpha_diversity_plots: {str(e)}")
            print(traceback.format_exc())
            
            # Return empty figures with error message
            error_fig1 = go.Figure()
            error_fig1.add_annotation(
                text=f"Error creating plot: {str(e)}",
                xref="paper", yref="paper",
                x=0.5, y=0.5,
                showarrow=False
            )
            error_fig1.update_layout(
                title='Alpha Diversity Measures Over Samples',
                plot_bgcolor='white',
                paper_bgcolor='white'
            )
            
            error_fig2 = go.Figure()
            error_fig2.add_annotation(
                text=f"Error creating plot: {str(e)}",
                xref="paper", yref="paper",
                x=0.5, y=0.5,
                showarrow=False
            )
            error_fig2.update_layout(
                title='Distribution of Alpha Diversity Measures',
                plot_bgcolor='white',
                paper_bgcolor='white'
            )
            
            return error_fig1, error_fig2

    def create_summary_statistics_layout(self, summary_stats):
        """Create the layout for summary statistics."""
        if not summary_stats:
            return dmc.Paper(
                children=[
                    dmc.Title("Summary Statistics", order=2),
                    dmc.Text("No data available", color="dimmed")
                ],
                p="lg",
                radius="md",
                withBorder=True
            )

        # Calculate trend data if possible
        trend_data = self.calculate_diversity_trends()

        # Prepare the metric information
        metrics_data = []
        
        # Format Shannon Diversity
        if 'shannon' in summary_stats and summary_stats['shannon']['mean'] is not None:
            shannon_row = {
                "Metric": "Shannon Diversity (H′)",
                "Mean": f"{summary_stats['shannon']['mean']:.2f}",
                "Std": f"{summary_stats['shannon']['std']:.2f}",
                "Min": f"{summary_stats['shannon']['min']:.2f}",
                "Max": f"{summary_stats['shannon']['max']:.2f}",
                "Trend": self._get_trend_icon(trend_data.get('shannon', {}).get('trend', 0))
            }
            metrics_data.append(shannon_row)
            
        # Format Simpson Diversity
        if 'simpson' in summary_stats and summary_stats['simpson']['mean'] is not None:
            simpson_row = {
                "Metric": "Simpson Diversity (D)",
                "Mean": f"{summary_stats['simpson']['mean']:.2f}",
                "Std": f"{summary_stats['simpson']['std']:.2f}",
                "Min": f"{summary_stats['simpson']['min']:.2f}",
                "Max": f"{summary_stats['simpson']['max']:.2f}",
                "Trend": self._get_trend_icon(trend_data.get('simpson', {}).get('trend', 0))
            }
            metrics_data.append(simpson_row)
            
        # Format Chao1 Richness
        if 'chao1' in summary_stats and summary_stats['chao1']['mean'] is not None:
            chao1_row = {
                "Metric": "Chao1 Richness",
                "Mean": f"{summary_stats['chao1']['mean']:.2f}",
                "Std": f"{summary_stats['chao1']['std']:.2f}",
                "Min": f"{summary_stats['chao1']['min']:.2f}",
                "Max": f"{summary_stats['chao1']['max']:.2f}",
                "Trend": self._get_trend_icon(trend_data.get('chao1', {}).get('trend', 0))
            }
            metrics_data.append(chao1_row)
            
        # Format ACE Diversity
        if 'ace' in summary_stats and summary_stats['ace']['mean'] is not None:
            ace_row = {
                "Metric": "ACE Diversity",
                "Mean": f"{summary_stats['ace']['mean']:.2f}",
                "Std": f"{summary_stats['ace']['std']:.2f}",
                "Min": f"{summary_stats['ace']['min']:.2f}",
                "Max": f"{summary_stats['ace']['max']:.2f}",
                "Trend": self._get_trend_icon(trend_data.get('ace', {}).get('trend', 0))
            }
            metrics_data.append(ace_row)
            
        # Format Faith's PD
        if 'faith_pd' in summary_stats and summary_stats['faith_pd']['mean'] is not None:
            faith_pd_row = {
                "Metric": "Faith's PD",
                "Mean": f"{summary_stats['faith_pd']['mean']:.2f}",
                "Std": f"{summary_stats['faith_pd']['std']:.2f}",
                "Min": f"{summary_stats['faith_pd']['min']:.2f}",
                "Max": f"{summary_stats['faith_pd']['max']:.2f}",
                "Trend": self._get_trend_icon(trend_data.get('faith_pd', {}).get('trend', 0))
            }
            metrics_data.append(faith_pd_row)
            
        # Format Species Richness
        if ('total_species' in summary_stats and summary_stats['total_species'] is not None and
            'mean_species_per_sample' in summary_stats and summary_stats['mean_species_per_sample'] is not None):
            species_row = {
                "Metric": "Species Richness",
                "Mean": f"{summary_stats['mean_species_per_sample']:.1f}",
                "Std": f"{summary_stats['std_species_per_sample']:.1f}" if 'std_species_per_sample' in summary_stats else "N/A",
                "Min": f"{summary_stats['min_species_per_sample']:.0f}" if 'min_species_per_sample' in summary_stats else "N/A",
                "Max": f"{summary_stats['max_species_per_sample']:.0f}" if 'max_species_per_sample' in summary_stats else "N/A",
                "Trend": self._get_trend_icon(trend_data.get('species_richness', {}).get('trend', 0))
            }
            metrics_data.append(species_row)
            
        # Create the metrics table
        metrics_table = dash_table.DataTable(
            data=metrics_data,
            columns=[
                {"name": "Metric", "id": "Metric"},
                {"name": "Mean", "id": "Mean"},
                {"name": "Std Dev", "id": "Std"},
                {"name": "Min", "id": "Min"},
                {"name": "Max", "id": "Max"},
                {"name": "Trend", "id": "Trend", "presentation": "markdown"}
            ],
            style_cell={
                'textAlign': 'left',
                'padding': '8px',
                'font-family': 'Arial',
                'fontSize': '14px'
            },
            style_header={
                'backgroundColor': '#f8f9fa',
                'fontWeight': 'bold',
                'border': '1px solid #ddd'
            },
            style_data={
                'border': '1px solid #ddd'
            },
            style_table={
                'overflowX': 'auto',
                'width': '100%'
            },
            markdown_options={"html": True}
        )
        
        # Create sample statistics card
        sample_stats_items = []
        
        if 'total_samples' in summary_stats and summary_stats['total_samples'] is not None:
            sample_stats_items.append(dmc.Group([
                dmc.Text("Total Samples:"),
                dmc.Text(str(summary_stats['total_samples']))
            ]))
            
        if 'total_projects' in summary_stats and summary_stats['total_projects'] is not None:
            sample_stats_items.append(dmc.Group([
                dmc.Text("Projects:"),
                dmc.Text(str(summary_stats['total_projects']))
            ]))
            
        if 'date_range' in summary_stats and summary_stats['date_range'] is not None:
            sample_stats_items.append(dmc.Group([
                dmc.Text("Date Range:"),
                dmc.Text(f"{summary_stats['date_range']} days")
            ]))
                
        # Add trend insights
        trend_insights = self._generate_trend_insights(trend_data)
        
        return dmc.Paper(
            children=[
                dmc.Title("Diversity Statistics Overview", order=2),
                
                # Metrics Table
                dmc.Text("Diversity Metrics", size="lg"),
                html.Div(metrics_table, style={"marginBottom": "20px"}),
                
                # Sample Statistics
                dmc.Grid(
                    [
                        dmc.Text("Sample Statistics", size="lg"),
                        dmc.Stack(sample_stats_items),
                        
                        dmc.Text("Trend Analysis", size="lg"),
                        dmc.Text(trend_insights)
                    ],
                ),
            ],
            p="lg", 
            radius="md",
            withBorder=True
        )

    def _get_trend_icon(self, trend_value):
        """Return a markdown formatted trend indicator based on the trend value."""
        if trend_value > 0.1:
            return "<span style='color: #28a745;'>▲ <b>Increasing</b></span>"
        elif trend_value < -0.1:
            return "<span style='color: #dc3545;'>▼ <b>Decreasing</b></span>"
        else:
            return "<span style='color: #6c757d;'>■ <b>Stable</b></span>"
    
    def _generate_trend_insights(self, trend_data):
        """Generate natural language insights based on trend data."""
        if not trend_data:
            return "Not enough data to analyze trends over time."
            
        insights = []
        
        # Shannon trends
        if 'shannon' in trend_data:
            shannon_trend = trend_data['shannon']['trend']
            if shannon_trend > 0.1:
                insights.append("Shannon diversity is increasing over time, suggesting growing microbial diversity and evenness.")
            elif shannon_trend < -0.1:
                insights.append("Shannon diversity is decreasing, indicating potential reduction in community complexity.")
                
        # Simpson trends
        if 'simpson' in trend_data:
            simpson_trend = trend_data['simpson']['trend']
            if simpson_trend > 0.1:
                insights.append("Simpson diversity index is rising, showing improved evenness in species distribution.")
            elif simpson_trend < -0.1:
                insights.append("Declining Simpson diversity suggests increasing dominance by fewer species.")
                
        # Species richness trends
        if 'species_richness' in trend_data:
            richness_trend = trend_data['species_richness']['trend']
            if richness_trend > 0.1:
                insights.append("Species richness is increasing, with more unique species being observed over time.")
            elif richness_trend < -0.1:
                insights.append("Species richness is declining, with fewer unique species detected in recent samples.")
        
        if not insights:
            insights.append("Diversity metrics appear relatively stable over the observed time period.")
            
        return " ".join(insights)
        
    def calculate_diversity_trends(self):
        """Calculate trends in diversity metrics over time."""
        try:
            # Calculate alpha diversity
            diversity_df = self.calculate_alpha_diversity()
            if diversity_df is None or diversity_df.empty:
                return {}
                
            # Add date information
            dates = {}
            for sample in diversity_df.index:
                if sample in self.sample_to_date_dict:
                    dates[sample] = self.sample_to_date_dict[sample]
            
            if not dates:
                return {}
                
            diversity_df['date'] = pd.Series(dates)
            
            # Sort by date
            diversity_df = diversity_df.sort_values('date')
            
            # Need at least 3 data points for meaningful trend
            if len(diversity_df) < 3:
                return {}
                
            # Calculate trends for each metric
            trends = {}
            
            for metric in ['shannon', 'simpson', 'chao1', 'ace', 'faith_pd']:
                if metric in diversity_df.columns:
                    # Get numeric x values for trend calculation
                    x_numeric = np.arange(len(diversity_df))
                    
                    if len(diversity_df[metric].dropna()) > 2:
                        # Calculate simple linear regression for trend
                        try:
                            slope, intercept, r_value, p_value, std_err = stats.linregress(
                                x_numeric, 
                                diversity_df[metric].values
                            )
                            
                            # Normalize trend by mean value
                            mean_value = diversity_df[metric].mean()
                            if mean_value != 0:
                                normalized_slope = slope / mean_value
                            else:
                                normalized_slope = slope
                                
                            trends[metric] = {
                                'trend': normalized_slope * len(diversity_df),  # Scale by dataset size for more intuitive values
                                'p_value': p_value,
                                'start': diversity_df[metric].iloc[0],
                                'end': diversity_df[metric].iloc[-1],
                                'change_pct': ((diversity_df[metric].iloc[-1] - diversity_df[metric].iloc[0]) / 
                                           max(abs(diversity_df[metric].iloc[0]), 0.001)) * 100
                            }
                        except Exception as e:
                            print(f"Error calculating trend for {metric}: {str(e)}")
            
            # Calculate species richness trend
            try:
                species_per_sample = diversity_df.apply(lambda row: sum(row.drop('date') > 0), axis=1)
                x_numeric = np.arange(len(species_per_sample))
                
                if len(species_per_sample) > 2:
                    slope, intercept, r_value, p_value, std_err = stats.linregress(
                        x_numeric, 
                        species_per_sample.values
                    )
                    
                    mean_value = species_per_sample.mean()
                    if mean_value != 0:
                        normalized_slope = slope / mean_value
                    else:
                        normalized_slope = slope
                        
                    trends['species_richness'] = {
                        'trend': normalized_slope * len(species_per_sample),
                        'p_value': p_value,
                        'start': species_per_sample.iloc[0],
                        'end': species_per_sample.iloc[-1],
                        'change_pct': ((species_per_sample.iloc[-1] - species_per_sample.iloc[0]) / 
                                    max(abs(species_per_sample.iloc[0]), 0.001)) * 100
                    }
            except Exception as e:
                print(f"Error calculating species richness trend: {str(e)}")
            
            return trends
            
        except Exception as e:
            print(f"Error in calculate_diversity_trends: {str(e)}")
            return {}

    def _init_layout(self):
        # Initialize dropdowns only if we have data
        sample_options = []
        project_options = [{'value': 'ALL', 'label': 'All Projects'}]
        subproject_options = [{'value': 'ALL', 'label': 'All Subprojects'}]
        date_options = [{'value': 'ALL', 'label': 'All Dates'}]
        
        if hasattr(self, 'unique_sample_ids'):
            sample_options = [{'value': sample, 'label': sample} for sample in self.unique_sample_ids]
        
        if hasattr(self, 'unique_projects_ids'):
            project_options.extend([{'value': project_id, 'label': project_id} for project_id in self.unique_projects_ids])
        
        if hasattr(self, 'unique_subprojects'):
            subproject_options.extend([{'value': subproject, 'label': subproject} for subproject in self.unique_subprojects])
            
        if hasattr(self, 'unique_dates'):
            date_options.extend([{'value': str(date), 'label': str(date)} for date in self.unique_dates])

        # Calculate summary statistics
        summary_stats = self.calculate_summary_statistics() if hasattr(self, 'df_full_for_diversity') else None

        # Create components
        diversity_metrics_dropdown = dcc.Dropdown(
            id='diversity-metric-dropdown',  # Fixed ID
            options=[
                {'label': 'Shannon Index (H′)', 'value': 'shannon'},
                {'label': 'Simpson Index (D)', 'value': 'simpson'},
                {'label': 'Chao1', 'value': 'chao1'},
                {'label': 'ACE', 'value': 'ace'},
                {'label': "Faith's PD", 'value': 'faith_pd'}
            ],
            value=['shannon', 'simpson'],
            multi=True,
            placeholder="Select diversity metrics",
            style={"width": "100%", "marginBottom": "10px"}
        )

        # Add download button
        download_button = dmc.Button(
            "Download Diversity Data",
            id="download-diversity-csv-btn",
            leftIcon=[DashIconify(icon="mdi:download")],
            color="blue",
            style={"marginBottom": "15px"}
        )
        
        # Add download component
        diversity_download = dcc.Download(id="download-diversity-csv")

        sample_select_dropdown = dcc.Dropdown(
            id='sample-select-value-diversity',  # Fixed ID
            options=sample_options,
            value=[opt['value'] for opt in sample_options] if sample_options else [],
            multi=True,
            placeholder="Select samples to display",
            style={"width": "100%"}
        )

        project_dropdown = dcc.Dropdown(
            id='project-dropdown-diversity',
            options=project_options,
            value='ALL',
            placeholder="Select Project",
            style={"width": "200px"}
        )

        subproject_dropdown = dcc.Dropdown(
            id='subproject-dropdown-diversity',
            options=subproject_options,
            value='ALL',
            placeholder="Select Subproject",
            style={"width": "200px"}
        )

        date_dropdown = dcc.Dropdown(
            id='date-dropdown-diversity',
            options=date_options,
            value='ALL',
            placeholder="Select Date",
            style={"width": "200px"}
        )

        taxonomic_rank_dropdown = dcc.Dropdown(
            id='taxonomic-rank-dropdown',
            options=[
                {'value': 'species', 'label': 'Species'},
                {'value': 'genus', 'label': 'Genus'},
                {'value': 'family', 'label': 'Family'},
                {'value': 'class', 'label': 'Class'},
                {'value': 'order', 'label': 'Order'},
                {'value': 'phylum', 'label': 'Phylum'}
            ],
            value='species',
            placeholder="Select Taxonomic Rank",
            style={"width": "200px"}
        )

        color_by_dropdown = dcc.Dropdown(
            id='color-by-dropdown',
            options=[
                {'label': 'Project', 'value': 'project'},
                {'label': 'Subproject', 'value': 'subproject'},
                {'label': 'Date', 'value': 'date'},
                {'label': 'K-means', 'value': 'kmeans'}
            ],
            value='project',
            placeholder="Color by",
            style={"width": "200px"}
        )

        kmeans_input = dmc.NumberInput(
            id='kmeans-input',
            label="Number of clusters",
            value=3,
            min=2,
            max=10,
            style={"width": "150px"}
        )

        # PCoA plot controls
        pcoa_controls = dmc.Group([
            dmc.Switch(
                id="arrow-toggle",
                label="Show Arrows",
                size="sm",
                checked=False
            ),
            dmc.Switch(
                id="pcoa-3d-toggle",
                label="3D Plot",
                size="sm",
                checked=False
            )
        ])

        # Create empty figures for initial load
        empty_fig = go.Figure()
        empty_fig.update_layout(
            plot_bgcolor='white',
            paper_bgcolor='white'
        )

        # Create summary statistics section
        summary_section = self.create_summary_statistics_layout(summary_stats)

        layout = dmc.Container(
            children=[
                dmc.Stack([
                    dmc.Title("Diversity Analysis", order=2),
                    
                    # Add introduction alert
                    dmc.Alert(
                        title="About Diversity Analysis",
                        children=[
                            "Explore different aspects of microbial diversity in your samples:",
                            html.Br(),
                            html.Ul([
                                html.Li("Alpha Diversity: Measures diversity within individual samples"),
                                html.Li("Beta Diversity: Shows relationships between samples"),
                                html.Li("Diversity Over Time: Track changes in diversity metrics across different dates"),
                                html.Li("Rarefaction Curves: Evaluate sampling depth sufficiency"),
                                html.Li("Use the controls above to filter samples and customize visualizations")
                            ])
                        ],
                        color="blue",
                        variant="light",
                        style={"marginBottom": "20px"}
                    ),
                    
                    # Add summary statistics section
                    summary_section,
                    
                    dmc.Space(h=20),
                    
                    dmc.Grid(
                        children=[
                            dmc.Group([
                                dmc.Stack([
                                    dmc.Text("Sample Selection"),
                                    sample_select_dropdown
                                ]),
                                dmc.Stack([
                                    dmc.Text("Diversity Metrics"),
                                    diversity_metrics_dropdown
                                ])
                            ], align="flex-start", style={"width": "100%"}),
                            
                            dmc.Space(h=10),
                            dmc.Group([
                                download_button,
                                diversity_download
                            ]),
                            
                            dmc.Space(h=20),
                            dmc.Group([   
                                dmc.Stack([
                                    dmc.Text("Filters"),
                                    dmc.Group([
                                        project_dropdown,
                                        subproject_dropdown,
                                        date_dropdown
                                    ])
                                ]),
                                dmc.Stack([
                                    dmc.Text("Visualization"),
                                    dmc.Group([
                                        taxonomic_rank_dropdown,
                                        color_by_dropdown,
                                        kmeans_input
                                    ])
                                ])
                            ], align="flex-start", style={"width": "100%"})
                        ],
                    ),
                    
                    dmc.Space(h=20),
                    
                    dmc.Paper([
                        dmc.Title("Alpha Diversity", order=3),
                        dmc.Grid(
                            children=[
                                dcc.Graph(id='alpha_diversities_plot1', figure=empty_fig),
                                dcc.Graph(id='alpha_diversities_plot2', figure=empty_fig)
                            ],
                        )
                    ], p="md", withBorder=True),
                    
                    dmc.Space(h=20),
                    
                    # Add new Diversity Over Time section
                    dmc.Paper([
                        dmc.Title("Diversity Over Time", order=3),
                        dcc.Graph(id='diversity_time_series', figure=empty_fig)
                    ], p="md", withBorder=True),
                    
                    dmc.Space(h=20),
                    
                    dmc.Paper([
                        dmc.Title("Beta Diversity", order=3),
                        dcc.Graph(id='beta_diversity_heatmap', figure=empty_fig)
                    ], p="md", withBorder=True),
                    
                    dmc.Space(h=20),
                    
                    dmc.Paper([
                        dmc.Title("PCoA Analysis", order=3),
                        pcoa_controls,
                        dcc.Graph(id='pcoa_plot_container', figure=empty_fig)
                    ], p="md", withBorder=True),
                    
                    dmc.Space(h=20),
                    
                    dmc.Paper([
                        dmc.Title("Rarefaction Curves", order=3),
                        dcc.Graph(id='rarefaction_plot', figure=empty_fig)
                    ], p="md", withBorder=True)
                ])
            ],
            fluid=True,
            style={'width': '100%', 'padding': '20px'}
        )

        self.app.layout = layout
        self.register_callbacks()

    def register_callbacks(self):
        """Register all callbacks for the Dash application."""
        # Remove callback registration from here since it's handled in index.py
        pass

    @_benchmark
    def update_diversity_plots(self, selected_samples, diversity_metric, color_by, k_value, project, subproject, date, is_3d, show_arrows, tax_rank="species"):
        print(f"\n=== Starting Diversity Plot Updates ===")
        print(f"Selected samples: {selected_samples}")
        print(f"Taxonomic rank: {tax_rank}")
        
        try:
            # Update the data frame based on taxonomic rank
            self.update_taxonomic_rank(tax_rank)
            
            # Use only selected samples for all plots
            if not selected_samples:
                selected_samples = self.unique_sample_ids
                
            # Create alpha diversity plots
            alpha_fig1, alpha_fig2 = self.create_alpha_diversity_plots(selected_samples)
            
            # Create diversity time series plot
            time_series_fig = self.create_diversity_time_series(selected_samples)
            
            # Create beta diversity heatmap
            beta_fig = self.create_beta_diversity_heatmap(selected_samples)
            
            # Create PCoA plot
            pcoa_fig = self.create_pcoa_figure(selected_samples, color_by, k_value, is_3d, show_arrows)
            
            # Create rarefaction plot
            rarefaction_data = self.calculate_rarefaction_curves(selected_samples)
            rarefaction_fig = self.create_rarefaction_plot(rarefaction_data)
            
            # Convert plots to a list before returning
            plots = [alpha_fig1, alpha_fig2, time_series_fig, beta_fig, pcoa_fig, rarefaction_fig]
            return plots
            
        except Exception as e:
            print(f"[ERROR] Error in update_diversity_plots: {str(e)}")
            print(traceback.format_exc())
            empty_fig = go.Figure()
            empty_fig.add_annotation(
                text=f"Error: {str(e)}",
                xref="paper", yref="paper",
                x=0.5, y=0.5,
                showarrow=False
            )
            return [empty_fig for _ in range(6)]

    def _create_empty_figure(self, error_message="No data available"):
        """Create an empty figure with an error message."""
        return go.Figure().add_annotation(
            text=error_message,
            xref="paper",
            yref="paper",
            x=0.5,
            y=0.5,
            showarrow=False
        )

    def update_taxonomic_rank(self, rank):
        """Update the data frame used for diversity calculations based on the selected taxonomic rank."""
        try:
            print(f"\nUpdating taxonomic rank to {rank}")
            
            # First, let's look at the raw data
            print("\nRaw data sample:")
            print(self.df[['sample_id', 'abundance', f'tax_{rank}' if rank != 'taxonomy' else rank]].head())
            
            # Sum abundances for each taxa within each sample
            grouped = self.df.groupby(['sample_id', f'tax_{rank}' if rank != 'taxonomy' else rank], observed=True)['abundance'].sum().reset_index()
            print("\nGrouped data sample:")
            print(grouped.head())
            
            # Create pivot table with samples as rows and taxa as columns
            pivot_data = pd.pivot_table(
                grouped,
                values='abundance',
                index='sample_id',
                columns=f'tax_{rank}' if rank != 'taxonomy' else rank,
                fill_value=0,
                aggfunc='sum',  # Explicitly specify sum aggregation
                observed=True   # Add observed=True to fix deprecation warning
            )
            
            print(f"\nPivot table shape: {pivot_data.shape}")
            print("Sample of pivot data:")
            print(pivot_data.head())
            print("\nColumn names:")
            print(pivot_data.columns[:5])  # Show first 5 column names
            print(f"\nSum of values: {pivot_data.values.sum()}")
            print("Row sums:")
            print(pivot_data.sum(axis=1).head())
            
            self.df_full_for_diversity = pivot_data
            print("Data update complete")
            
        except Exception as e:
            print(f"Error in update_taxonomic_rank: {str(e)}")
            traceback.print_exc()

    def calculate_alpha_diversity(self, use_normalized_counts=False):
        """Calculate alpha diversity metrics with proper error handling and data validation."""
        try:
            print("\nStarting alpha diversity calculation...")
            if self.df_full_for_diversity is None or self.df_full_for_diversity.empty:
                print("No data available for alpha diversity calculation")
                return pd.DataFrame()

            counts = self.df_full_for_diversity.values
            results = {
                'shannon': [],
                'simpson': [],
                'chao1': [],
                'ace': [],
                'faith_pd': []
            }

            for i, sample_id in enumerate(self.df_full_for_diversity.index):
                sample_counts = counts[i]
                total_counts = np.sum(sample_counts)
                
                if total_counts > 0:
                    proportions = sample_counts / total_counts
                    nonzero_counts = sample_counts[sample_counts > 0]
                    
                    if len(proportions) > 0:
                        # Shannon diversity
                        shannon = -np.sum(proportions[proportions > 0] * np.log(proportions[proportions > 0]))
                        
                        # Simpson diversity
                        simpson = 1 - np.sum(proportions ** 2)
                        
                        # Chao1 estimator
                        singletons = np.sum(sample_counts == 1)
                        doubletons = np.sum(sample_counts == 2)
                        observed_species = np.sum(sample_counts > 0)
                        if doubletons > 0:
                            chao1 = observed_species + (singletons * (singletons - 1)) / (2 * (doubletons + 1))
                        else:
                            chao1 = observed_species + (singletons * (singletons - 1)) / 2
                        
                        # ACE estimator (improved)
                        rare_threshold = 10
                        rare_counts = nonzero_counts[nonzero_counts <= rare_threshold]
                        abundant_counts = nonzero_counts[nonzero_counts > rare_threshold]
                        
                        if len(rare_counts) > 0:
                            s_rare = len(rare_counts)
                            n_rare = np.sum(rare_counts)
                            f1 = np.sum(rare_counts == 1)
                            
                            # Calculate ACE components
                            c_ace = 1 - (f1 / n_rare)
                            if c_ace > 0:
                                # Calculate gamma component
                                f_rare = np.array([np.sum(rare_counts == i) for i in range(1, rare_threshold + 1)])
                                s_rare_squared = s_rare * (s_rare - 1)
                                if s_rare_squared > 0:
                                    gamma_ace = (s_rare * np.sum(f_rare * (np.arange(1, len(f_rare) + 1)) * 
                                              (np.arange(1, len(f_rare) + 1) - 1))) / (c_ace * n_rare * (n_rare - 1))
                                    gamma_ace = max(0, gamma_ace - 1)
                                    
                                    # Final ACE calculation
                                    ace = len(abundant_counts) + s_rare/c_ace + (f1/c_ace) * gamma_ace
                                else:
                                    ace = len(nonzero_counts)
                            else:
                                ace = len(nonzero_counts)
                        else:
                            ace = len(nonzero_counts)
                        
                        # Faith's PD (improved approximation)
                        # Using branch lengths approximated by abundance differences
                        sorted_props = np.sort(proportions[proportions > 0])
                        branch_lengths = np.diff(np.concatenate([[0], sorted_props]))
                        faith_pd = np.sum(branch_lengths) * observed_species
                        
                        results['shannon'].append(shannon)
                        results['simpson'].append(simpson)
                        results['chao1'].append(chao1)
                        results['ace'].append(ace)
                        results['faith_pd'].append(faith_pd)
                    else:
                        results['shannon'].append(0)
                        results['simpson'].append(0)
                        results['chao1'].append(0)
                        results['ace'].append(0)
                        results['faith_pd'].append(0)
                else:
                    results['shannon'].append(0)
                    results['simpson'].append(0)
                    results['chao1'].append(0)
                    results['ace'].append(0)
                    results['faith_pd'].append(0)

            diversity_df = pd.DataFrame(
                results,
                index=self.df_full_for_diversity.index
            )
            
            return diversity_df

        except Exception as e:
            print(f"Error in alpha diversity calculation: {str(e)}")
            traceback.print_exc()
            return pd.DataFrame()

    def calculate_beta_diversity(self):
        """Calculate beta diversity with proper error handling."""
        try:
            print("Calculating beta diversity...")
            if self.df_full_for_diversity is None or self.df_full_for_diversity.empty:
                print("No data available for beta diversity calculation")
                return None

            # Convert DataFrame to numpy array for beta diversity calculation
            abundance_data = self.df_full_for_diversity.values.astype(float)
            sample_ids = list(self.df_full_for_diversity.index)

            # Calculate beta diversity using Bray-Curtis distance
            beta_div = skbio.diversity.beta_diversity('braycurtis', abundance_data, sample_ids)
            
            # Convert to numpy array and ensure float type
            distance_matrix = beta_div.data.astype(float)
            
            # Ensure the matrix is symmetric and contains valid distances
            if not np.allclose(distance_matrix, distance_matrix.T):
                print("Warning: Distance matrix is not symmetric, fixing...")
                distance_matrix = (distance_matrix + distance_matrix.T) / 2
            
            if not np.all(np.isfinite(distance_matrix)):
                print("Warning: Distance matrix contains non-finite values")
                distance_matrix[~np.isfinite(distance_matrix)] = 1.0
            
            # Create a new DistanceMatrix object with the cleaned data
            beta_div = DistanceMatrix(distance_matrix, ids=sample_ids)
            
            return beta_div

        except Exception as e:
            print(f"Error in calculate_beta_diversity: {str(e)}")
            print(traceback.format_exc())
            return None

    def get_data(self):
        records = NanoporeRecord.objects.filter(user_id=self.user_id)
        self.records = records
        if not records.exists():
            return pd.DataFrame()
        df = pd.DataFrame.from_records(records.values())
        df_sorted = df.sort_values(by=['project_id', 'subproject', 'sample_id'])
        return df_sorted

    def calculate_normalized_counts(self):
        counts_df = pd.DataFrame.from_records(NanoporeRecord.objects.filter(user_id=self.user_id).values())
        stats_df = pd.DataFrame.from_records(SequencingStatistics.objects.filter(user_id=self.user_id).values())

        try:
            merged_df = pd.merge(counts_df, stats_df, left_on='sample_id', right_on="sample_name")
        except KeyError as e:
            print(f"KeyError occurred: {e}")
            merged_df = counts_df
        try:
            merged_df['normalized_count'] = (merged_df['count'] / merged_df['total_bases_sequenced']) * 10000000
        except KeyError as e:
            print(
                f"KeyError occurred. Check if total_bases_sequenced is in merged_df setting normalized count to not normalized count for now: {e}")
            merged_df['normalized_count'] = (merged_df['count'])

        normalized_counts_df = merged_df[['sample_id', 'taxonomy', 'abundance', 'count', 'normalized_count']]
        return normalized_counts_df

    def create_beta_diversity_heatmap(self, selected_samples):
        try:
            if self.beta_diversity_matrix is None:
                print("No beta diversity matrix available")
                return go.Figure()

            df_clean = distance_matrix_to_dataframe(self.beta_diversity_matrix)
            
            # Check if selected samples exist in the matrix
            valid_samples = [s for s in selected_samples if s in df_clean.index and s in df_clean.columns]
            if not valid_samples:
                print("No valid samples found in beta diversity matrix")
                return go.Figure()
                
            df_clean = df_clean.loc[valid_samples, valid_samples]
            
            # Validate data before plotting
            if df_clean.empty or df_clean.isnull().values.all():
                print("Beta diversity matrix contains invalid data")
                return go.Figure()

            # Ensure numeric type and handle potential bool values
            df_clean = df_clean.astype(float)
            
            # Check for any remaining NaN values
            if df_clean.isnull().any().any():
                print("Beta diversity matrix contains NaN values")
                return go.Figure()
            
            # Create annotated heatmap
            beta_fig = go.Figure(data=go.Heatmap(
                z=df_clean.values,
                x=df_clean.columns,
                y=df_clean.index,
                colorscale='RdBu_r',
                zmin=0,
                zmax=1,
                colorbar=dict(title='Bray-Curtis Distance')
            ))

            beta_fig.update_layout(
                title='Beta Diversity Heatmap',
                xaxis_title='Sample ID',
                yaxis_title='Sample ID',
                height=1000,
                width=1200,
                margin=dict(l=100, r=100, b=100, t=100)
            )

            return beta_fig

        except Exception as e:
            print(f"Error in create_beta_diversity_heatmap: {str(e)}")
            print(traceback.format_exc())
            return go.Figure()

    def calculate_rarefaction_curves(self, selected_samples, max_depth=None, num_points=20):
        """Calculate rarefaction curves for selected samples."""
        try:
            print("\nStarting rarefaction calculation...")
            # Get abundance data for selected samples
            abundance_data = self.df_full_for_diversity.loc[selected_samples]
            
            if abundance_data.empty:
                print("No abundance data available")
                return None
            
            print(f"Abundance data shape: {abundance_data.shape}")
            print("Sample of abundance data:")
            print(abundance_data.head())
            
            # Get row sums (total abundance per sample)
            row_sums = abundance_data.sum(axis=1)
            print("\nTotal abundances per sample:")
            print(row_sums.head())
            
            # Determine max depth if not provided
            if max_depth is None:
                max_depth = row_sums.min()  # Use minimum total abundance
                if max_depth < 100:  # If very low counts, use the maximum instead
                    max_depth = row_sums.max()
            
            print(f"Using max depth: {max_depth}")
            
            # Create depth points (use more points at lower depths where changes are more dramatic)
            depths = np.concatenate([
                np.linspace(10, min(1000, max_depth/2), num_points//2, dtype=int),
                np.linspace(min(1000, max_depth/2), max_depth, num_points//2, dtype=int)
            ])
            depths = np.unique(depths)  # Remove any duplicates
            
            # Calculate rarefaction curves
            rarefaction_data = pd.DataFrame(index=depths)
            
            for sample in selected_samples:
                try:
                    print(f"\nProcessing sample: {sample}")
                    sample_abundances = abundance_data.loc[sample]
                    total_abundance = sample_abundances.sum()
                    
                    if total_abundance == 0:
                        print(f"Skipping sample {sample} - no abundance")
                        continue
                    
                    # Get non-zero species and their abundances
                    nonzero_species = sample_abundances[sample_abundances > 0]
                    print(f"Number of species with non-zero abundance: {len(nonzero_species)}")
                    print(f"Total abundance: {total_abundance:.2f}")
                    
                    # Calculate probabilities
                    probs = nonzero_species / total_abundance
                    
                    # Initialize list for this sample's rarefaction curve
                    sample_richness = []
                    
                    # Calculate richness at each depth
                    for depth in depths:
                        if depth > total_abundance:
                            # If depth is greater than total abundance, use the total number of species
                            richness = len(nonzero_species)
                        else:
                            # Perform multiple rarefactions and take the mean
                            richness_values = []
                            n_iterations = 10  # Number of random subsamples to average over
                            
                            for _ in range(n_iterations):
                                # Randomly subsample the community at this depth
                                subsample = np.random.multinomial(depth, probs)
                                # Count number of non-zero species
                                richness = np.sum(subsample > 0)
                                richness_values.append(richness)
                            
                            # Take the mean richness at this depth
                            richness = np.mean(richness_values)
                        
                        sample_richness.append(richness)
                        if len(sample_richness) <= 5:  # Print first few points
                            print(f"Depth {depth}: richness = {richness:.2f}")
                    
                    # Add this sample's curve to the dataframe
                    rarefaction_data[sample] = sample_richness
                    print(f"Completed rarefaction for sample {sample}")
                    
                except Exception as e:
                    print(f"Error processing sample {sample}: {str(e)}")
                    continue
        
            print("\nRarefaction calculation complete")
            print("Final rarefaction data shape:", rarefaction_data.shape)
            return rarefaction_data
        
        except Exception as e:
            print(f"Error in rarefaction calculation: {str(e)}")
            traceback.print_exc()
            return None

    def create_rarefaction_plot(self, rarefaction_data):
        """Create a rarefaction curve plot."""
        try:
            if rarefaction_data is None or rarefaction_data.empty:
                return self._create_empty_figure("No rarefaction data available")
                
            # Create figure
            fig = go.Figure()
            
            # Add a line for each sample
            for sample in rarefaction_data.columns:
                fig.add_trace(
                    go.Scatter(
                        x=rarefaction_data.index,
                        y=rarefaction_data[sample],
                        name=sample,
                        mode='lines',
                        hovertemplate='Depth: %{x}<br>Species: %{y:.0f}<extra></extra>'
                    )
                )
                
            # Update layout
            fig.update_layout(
                title='Rarefaction Curves',
                xaxis_title='Sequencing Depth',
                yaxis_title='Number of Observed Species',
                showlegend=True,
                plot_bgcolor='white',
                paper_bgcolor='white',
                height=500,
                # Add a note about methodology
                annotations=[
                    dict(
                        text='Each curve shows species accumulation with increased sampling depth',
                        xref='paper',
                        yref='paper',
                        x=0,
                        y=1.1,
                        showarrow=False,
                        font=dict(size=10)
                    )
                ]
            )
            
            return fig
            
        except Exception as e:
            print(f"Error in rarefaction plotting: {str(e)}")
            return self._create_empty_figure(str(e))

    def calculate_summary_statistics(self):
        """Calculate summary statistics for the diversity data."""
        if self.df_full_for_diversity is None or self.df_full_for_diversity.empty:
            return {
                'shannon': {'mean': None, 'std': None, 'min': None, 'max': None, 'median': None},
                'simpson': {'mean': None, 'std': None, 'min': None, 'max': None, 'median': None},
                'chao1': {'mean': None, 'std': None, 'min': None, 'max': None, 'median': None},
                'ace': {'mean': None, 'std': None, 'min': None, 'max': None, 'median': None},
                'faith_pd': {'mean': None, 'std': None, 'min': None, 'max': None, 'median': None},
                'total_species': None,
                'mean_species_per_sample': None,
                'std_species_per_sample': None,
                'min_species_per_sample': None,
                'max_species_per_sample': None,
                'total_samples': None,
                'total_projects': None,
                'date_range': None
            }

        try:
            # Calculate diversity indices
            diversity_metrics = self.calculate_alpha_diversity()
            
            # Calculate comprehensive statistics for each diversity metric
            stats = {}
            for metric in ['shannon', 'simpson', 'chao1', 'ace', 'faith_pd']:
                if metric in diversity_metrics:
                    metric_data = diversity_metrics[metric]
                    stats[metric] = {
                        'mean': float(metric_data.mean()),
                        'std': float(metric_data.std()),
                        'min': float(metric_data.min()),
                        'max': float(metric_data.max()),
                        'median': float(metric_data.median())
                    }
                else:
                    stats[metric] = {
                        'mean': None, 'std': None, 'min': None, 'max': None, 'median': None
                    }
            
            # Calculate species counts with enhanced statistics
            species_counts = self.df_full_for_diversity.apply(lambda row: sum(1 for x in row if x > 0))
            
            # Add species statistics
            if not species_counts.empty:
                stats.update({
                    'total_species': int(species_counts.sum()),
                    'mean_species_per_sample': float(species_counts.mean()),
                    'std_species_per_sample': float(species_counts.std()),
                    'min_species_per_sample': float(species_counts.min()),
                    'max_species_per_sample': float(species_counts.max()),
                    'total_samples': len(self.df_full_for_diversity),
                    'total_projects': len(set(self.sample_to_project_dict.values()))
                })
            else:
                stats.update({
                    'total_species': None,
                    'mean_species_per_sample': None,
                    'std_species_per_sample': None,
                    'min_species_per_sample': None,
                    'max_species_per_sample': None,
                    'total_samples': 0,
                    'total_projects': 0
                })
            
            # Calculate date range if dates are available
            if self.sample_to_date_dict and len(self.sample_to_date_dict) > 0:
                dates = [d for d in self.sample_to_date_dict.values() if d is not None]
                if dates:
                    stats['date_range'] = (max(dates) - min(dates)).days
                else:
                    stats['date_range'] = None
            else:
                stats['date_range'] = None

            return stats
            
        except Exception as e:
            print(f"Error calculating summary statistics: {str(e)}")
            print(traceback.format_exc())
            return {
                'shannon': {'mean': None, 'std': None, 'min': None, 'max': None, 'median': None},
                'simpson': {'mean': None, 'std': None, 'min': None, 'max': None, 'median': None},
                'chao1': {'mean': None, 'std': None, 'min': None, 'max': None, 'median': None},
                'ace': {'mean': None, 'std': None, 'min': None, 'max': None, 'median': None},
                'faith_pd': {'mean': None, 'std': None, 'min': None, 'max': None, 'median': None},
                'total_species': None,
                'mean_species_per_sample': None,
                'std_species_per_sample': None,
                'min_species_per_sample': None,
                'max_species_per_sample': None,
                'total_samples': None,
                'total_projects': None,
                'date_range': None
            }

    @_benchmark
    def create_diversity_time_series(self, selected_samples):
        """Create a time series plot showing how diversity metrics change over time."""
        try:
            # Calculate diversity metrics for selected samples
            diversity_df = self.calculate_alpha_diversity()
            if diversity_df is None or diversity_df.empty:
                empty_fig = go.Figure()
                empty_fig.add_annotation(
                    text="No diversity data available for time series",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5,
                    showarrow=False
                )
                empty_fig.update_layout(
                    title='Diversity Metrics Over Time',
                    plot_bgcolor='white',
                    paper_bgcolor='white'
                )
                return empty_fig
            
            # Filter for selected samples
            available_samples = [s for s in selected_samples if s in diversity_df.index]
            if not available_samples:
                empty_fig = go.Figure()
                empty_fig.add_annotation(
                    text="Selected samples have no diversity data for time series",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5,
                    showarrow=False
                )
                empty_fig.update_layout(
                    title='Diversity Metrics Over Time',
                    plot_bgcolor='white',
                    paper_bgcolor='white'
                )
                return empty_fig
            
            diversity_df = diversity_df.loc[available_samples]
            
            # Add date information to the diversity DataFrame
            dates = {}
            for sample in diversity_df.index:
                if sample in self.sample_to_date_dict:
                    dates[sample] = self.sample_to_date_dict[sample]
            
            if not dates:
                empty_fig = go.Figure()
                empty_fig.add_annotation(
                    text="No date information available for samples",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5,
                    showarrow=False
                )
                empty_fig.update_layout(
                    title='Diversity Metrics Over Time',
                    plot_bgcolor='white',
                    paper_bgcolor='white'
                )
                return empty_fig
            
            # Add date column and project column
            diversity_df['date'] = pd.Series(dates)
            diversity_df['project'] = pd.Series({sample: self.sample_to_project_dict.get(sample, 'Unknown') 
                                              for sample in diversity_df.index})
            
            # Sort by date
            diversity_df = diversity_df.sort_values('date')
            
            # Create figure
            fig = go.Figure()
            
            # Add traces for Shannon and Simpson diversity
            if 'shannon' in diversity_df.columns:
                fig.add_trace(go.Scatter(
                    x=diversity_df['date'],
                    y=diversity_df['shannon'],
                    mode='lines+markers',
                    name='Shannon Diversity',
                    line=dict(color='#1f77b4', width=2),
                    marker=dict(size=8),
                    hovertemplate='<b>%{text}</b><br>Date: %{x|%Y-%m-%d}<br>Shannon: %{y:.3f}<extra></extra>',
                    text=diversity_df.index
                ))
            
            if 'simpson' in diversity_df.columns:
                fig.add_trace(go.Scatter(
                    x=diversity_df['date'],
                    y=diversity_df['simpson'],
                    mode='lines+markers',
                    name='Simpson Diversity',
                    line=dict(color='#ff7f0e', width=2, dash='dash'),
                    marker=dict(size=8),
                    hovertemplate='<b>%{text}</b><br>Date: %{x|%Y-%m-%d}<br>Simpson: %{y:.3f}<extra></extra>',
                    text=diversity_df.index
                ))
            
            # Add observed species richness if available
            try:
                # Try to calculate observed species (excluding date and project columns)
                if 'date' in diversity_df.columns and 'project' in diversity_df.columns:
                    cols_to_exclude = ['date', 'project']
                elif 'date' in diversity_df.columns:
                    cols_to_exclude = ['date']
                elif 'project' in diversity_df.columns:
                    cols_to_exclude = ['project']
                else:
                    cols_to_exclude = []
                
                # Calculate observed species by excluding metadata columns
                observed_species = diversity_df.loc[:, ~diversity_df.columns.isin(cols_to_exclude)].gt(0).sum(axis=1)
                diversity_df['observed_species'] = observed_species
                
                fig.add_trace(go.Scatter(
                    x=diversity_df['date'],
                    y=diversity_df['observed_species'],
                    mode='lines+markers',
                    name='Observed Species',
                    line=dict(color='#2ca02c', width=2),
                    marker=dict(size=8),
                    yaxis='y2',
                    hovertemplate='<b>%{text}</b><br>Date: %{x|%Y-%m-%d}<br>Observed Species: %{y}<extra></extra>',
                    text=diversity_df.index
                ))
            except Exception as e:
                print(f"Error calculating observed species: {str(e)}")
            
            # Update layout
            fig.update_layout(
                title='Diversity Metrics Over Time',
                xaxis=dict(
                    title='Date',
                    showgrid=True,
                    gridwidth=1,
                    gridcolor='lightgray',
                    type='date'
                ),
                yaxis=dict(
                    title='Diversity Index',
                    showgrid=True,
                    gridwidth=1,
                    gridcolor='lightgray',
                    side='left'
                ),
                yaxis2=dict(
                    title='Observed Species',
                    showgrid=False,
                    side='right',
                    overlaying='y'
                ),
                plot_bgcolor='white',
                paper_bgcolor='white',
                hovermode='closest',
                height=500,
                legend=dict(
                    orientation='h',
                    yanchor='bottom',
                    y=1.02,
                    xanchor='right',
                    x=1
                ),
                margin=dict(l=80, r=80, t=100, b=80)
            )
            
            # Add range slider
            fig.update_layout(
                xaxis=dict(
                    rangeselector=dict(
                        buttons=list([
                            dict(count=7, label="1w", step="day", stepmode="backward"),
                            dict(count=1, label="1m", step="month", stepmode="backward"),
                            dict(count=3, label="3m", step="month", stepmode="backward"),
                            dict(count=6, label="6m", step="month", stepmode="backward"),
                            dict(count=1, label="1y", step="year", stepmode="backward"),
                            dict(step="all")
                        ])
                    ),
                    rangeslider=dict(visible=True),
                    type="date"
                )
            )
            
            # Add annotations explaining the metrics
            fig.add_annotation(
                text="Shannon and Simpson diversity indices measure the richness and evenness of species",
                xref="paper", yref="paper",
                x=0, y=1.08,
                showarrow=False,
                font=dict(size=10, color="gray"),
                align="left"
            )
            
            # Add trend lines using LOWESS smoothing if we have enough points
            if len(diversity_df) >= 5:
                try:
                    if 'shannon' in diversity_df.columns:
                        # Get numeric x values for LOWESS
                        x_numeric = np.array([(d - diversity_df['date'].min()).total_seconds() for d in diversity_df['date']])
                        lowess_data = lowess(diversity_df['shannon'].values, x_numeric, frac=0.5)
                        
                        # Map back to dates for plotting
                        trend_dates = [diversity_df['date'].min() + pd.Timedelta(seconds=x) for x in lowess_data[:, 0]]
                        
                        fig.add_trace(go.Scatter(
                            x=trend_dates,
                            y=lowess_data[:, 1],
                            mode='lines',
                            name='Shannon Trend',
                            line=dict(color='#1f77b4', width=3, dash='dot'),
                            opacity=0.7,
                            showlegend=False
                        ))
                except Exception as e:
                    print(f"Error creating trend line: {str(e)}")
            
            return fig
            
        except Exception as e:
            print(f"Error in create_diversity_time_series: {str(e)}")
            print(traceback.format_exc())
            
            # Return empty figure with error message
            error_fig = go.Figure()
            error_fig.add_annotation(
                text=f"Error creating time series plot: {str(e)}",
                xref="paper", yref="paper",
                x=0.5, y=0.5,
                showarrow=False
            )
            error_fig.update_layout(
                title='Diversity Metrics Over Time',
                plot_bgcolor='white',
                paper_bgcolor='white'
            )
            
            return error_fig

    @_benchmark
    def create_pcoa_figure(self, selected_samples, color_by, k_value, is_3d=False, show_arrows=True):
        """Create PCoA plot using Bray-Curtis distance."""
        try:
            # Get abundance data for selected samples
            if not hasattr(self, 'df_full_for_diversity') or self.df_full_for_diversity is None:
                return self._create_empty_figure("No abundance data available for PCoA analysis")
                
            # Filter to selected samples
            available_samples = [s for s in selected_samples if s in self.df_full_for_diversity.index]
            if not available_samples:
                return self._create_empty_figure("Selected samples have no abundance data for PCoA analysis")
                
            abundance_data = self.df_full_for_diversity.loc[available_samples].values

            # Calculate Bray-Curtis distance
            try:
                distance_matrix = skbio.diversity.beta_diversity(
                    metric='braycurtis',
                    counts=abundance_data,
                    ids=available_samples
                )
            except Exception as e:
                print(f"Error calculating beta diversity: {str(e)}")
                return self._create_empty_figure(f"Error in beta diversity calculation: {str(e)}")

            # Perform PCoA
            try:
                pcoa_results = skbio.stats.ordination.pcoa(distance_matrix)
            except Exception as e:
                print(f"Error performing PCoA: {str(e)}")
                return self._create_empty_figure(f"Error in PCoA calculation: {str(e)}")
            
            # Convert to DataFrame
            df_pcoa = pd.DataFrame(
                data=pcoa_results.samples,
                columns=['PC1', 'PC2', 'PC3'],
                index=available_samples
            )
            
            # Add metadata for coloring
            df_pcoa['sample_id'] = df_pcoa.index
            
            if color_by == 'kmeans':
                # Use K-means clustering
                kmeans = KMeans(n_clusters=k_value, random_state=42)
                df_pcoa['cluster'] = kmeans.fit_predict(pcoa_results.samples[:, :3])
                color_column = 'cluster'
                df_pcoa['cluster'] = df_pcoa['cluster'].astype(str)  # Convert to string for consistent coloring
            elif color_by == 'project' and hasattr(self, 'sample_to_project_dict'):
                # Color by project
                df_pcoa['project'] = df_pcoa.index.map(lambda x: self.sample_to_project_dict.get(x, 'Unknown'))
                color_column = 'project'
            elif color_by == 'subproject' and hasattr(self, 'sample_to_subproject_dict'):
                # Color by subproject
                df_pcoa['subproject'] = df_pcoa.index.map(lambda x: self.sample_to_subproject_dict.get(x, 'Unknown'))
                color_column = 'subproject'
            elif color_by == 'date' and hasattr(self, 'sample_to_date_dict'):
                # Color by date
                df_pcoa['date'] = df_pcoa.index.map(lambda x: str(self.sample_to_date_dict.get(x, 'Unknown')))
                color_column = 'date'
            else:
                # Default coloring
                df_pcoa['group'] = 'All Samples'
                color_column = 'group'
                
            # Create the plot
            if is_3d:
                fig = px.scatter_3d(
                    df_pcoa,
                    x='PC1',
                    y='PC2',
                    z='PC3',
                    color=color_column,
                    title='PCoA Plot (3D)',
                    labels={
                        'PC1': f'PC1 ({pcoa_results.proportion_explained[0]:.2%})',
                        'PC2': f'PC2 ({pcoa_results.proportion_explained[1]:.2%})',
                        'PC3': f'PC3 ({pcoa_results.proportion_explained[2]:.2%})'
                    },
                    hover_name='sample_id'
                )
                
                # Update marker appearance for 3D plot
                fig.update_traces(
                    marker=dict(
                        size=8,
                        line=dict(width=1, color='DarkSlateGrey')
                    )
                )
            else:
                fig = px.scatter(
                    df_pcoa,
                    x='PC1',
                    y='PC2',
                    color=color_column,
                    title='PCoA Plot',
                    labels={
                        'PC1': f'PC1 ({pcoa_results.proportion_explained[0]:.2%})',
                        'PC2': f'PC2 ({pcoa_results.proportion_explained[1]:.2%})'
                    },
                    hover_name='sample_id'
                )
                
                # Update marker appearance for 2D plot
                fig.update_traces(
                    marker=dict(
                        size=10,
                        line=dict(width=1, color='DarkSlateGrey')
                    )
                )
                
                # Add arrows if requested (only for 2D plot)
                if show_arrows:
                    try:
                        # Calculate loadings (feature contributions to PCs)
                        features = self.df_full_for_diversity.columns
                        
                        # Only use top features for cleaner visualization
                        n_arrows = min(20, len(features))
                        
                        # Get the feature importance
                        feature_importance = np.abs(abundance_data).sum(axis=0)
                        top_features_idx = np.argsort(feature_importance)[-n_arrows:]
                        
                        # Get coordinates for arrows (simplified approach)
                        x_arrows = pcoa_results.samples[:, 0].mean()
                        y_arrows = pcoa_results.samples[:, 1].mean()
                        
                        # Add arrows for top features
                        for idx in top_features_idx:
                            feature_name = features[idx]
                            # Use correlation as direction with error handling
                            try:
                                # Check if we have enough data points for correlation
                                if len(abundance_data[:, idx]) < 2 or len(pcoa_results.samples[:, 0]) < 2:
                                    continue
                                    
                                # Check for NaN values in the data
                                if np.isnan(abundance_data[:, idx]).any() or np.isnan(pcoa_results.samples).any():
                                    continue
                                    
                                # Calculate correlations with error handling
                                try:
                                    corr_x = np.corrcoef(abundance_data[:, idx], pcoa_results.samples[:, 0])[0, 1]
                                    corr_y = np.corrcoef(abundance_data[:, idx], pcoa_results.samples[:, 1])[0, 1]
                                except Exception as e:
                                    print(f"Error calculating correlation for feature {feature_name}: {str(e)}")
                                    continue
                                
                                # Skip if correlation is too low or NaN
                                if np.isnan(corr_x) or np.isnan(corr_y) or (abs(corr_x) < 0.3 and abs(corr_y) < 0.3):
                                    continue
                                    
                                # Scale arrow length
                                arrow_length = np.sqrt(corr_x**2 + corr_y**2)
                                if arrow_length == 0:
                                    continue
                                    
                                # Normalize and scale
                                scale_factor = 0.2
                                dx = corr_x * scale_factor
                                dy = corr_y * scale_factor
                                
                                # Add the arrow
                                fig.add_annotation(
                                    x=x_arrows,
                                    y=y_arrows,
                                    ax=x_arrows + dx,
                                    ay=y_arrows + dy,
                                    arrowhead=2,
                                    arrowsize=1,
                                    arrowwidth=1,
                                    arrowcolor="rgba(0,0,0,0.5)",
                                    showarrow=True,
                                    text=feature_name,
                                    font=dict(size=10),
                                    hoverinfo="name",
                                )
                            except Exception as e:
                                # Silently skip features that cause errors
                                continue
                    except Exception as e:
                        print(f"Error adding arrows to PCoA plot: {str(e)}")
                
            # Update layout
            fig.update_layout(
                plot_bgcolor='white',
                paper_bgcolor='white',
                showlegend=True,
                legend=dict(
                    title=dict(text=color_by.capitalize()),
                    orientation='h',
                    yanchor='bottom',
                    y=1.02,
                    xanchor='right',
                    x=1
                ),
                margin=dict(l=50, r=50, t=80, b=50)
            )
            
            # Add explanatory annotation
            fig.add_annotation(
                text="Principal Coordinate Analysis (PCoA) visualizes sample relationships based on their taxonomic composition",
                xref="paper", yref="paper",
                x=0, y=1.07,
                showarrow=False,
                font=dict(size=10, color="gray"),
                align="left"
            )
            
            return fig
            
        except Exception as e:
            print(f"Error in create_pcoa_figure: {str(e)}")
            print(traceback.format_exc())
            return self._create_empty_figure(f"Error creating PCoA plot: {str(e)}")

    def prepare_diversity_csv_data(self):
        """Prepare diversity data for CSV export, including all metrics and statistics."""
        try:
            # Calculate alpha diversity for all samples
            diversity_df = self.calculate_alpha_diversity()
            
            if diversity_df is None or diversity_df.empty:
                return pd.DataFrame()
            
            # Add sample metadata
            diversity_df['sample_id'] = diversity_df.index
            
            # Add project information if available
            if hasattr(self, 'sample_to_project_dict'):
                diversity_df['project'] = diversity_df.index.map(
                    lambda x: self.sample_to_project_dict.get(x, 'Unknown')
                )
            
            # Add subproject information if available
            if hasattr(self, 'sample_to_subproject_dict'):
                diversity_df['subproject'] = diversity_df.index.map(
                    lambda x: self.sample_to_subproject_dict.get(x, 'Unknown')
                )
            
            # Add date information if available
            if hasattr(self, 'sample_to_date_dict'):
                diversity_df['date'] = diversity_df.index.map(
                    lambda x: self.sample_to_date_dict.get(x, 'Unknown')
                )
            
            # Calculate observed species for each sample
            try:
                if self.df_full_for_diversity is not None:
                    observed_species = self.df_full_for_diversity.apply(lambda row: (row > 0).sum(), axis=1)
                    diversity_df['observed_species'] = diversity_df.index.map(
                        lambda x: observed_species.get(x, 0) if x in observed_species.index else 0
                    )
            except Exception as e:
                print(f"Error calculating observed species: {str(e)}")
            
            # Calculate summary statistics
            stats = self.calculate_summary_statistics()
            
            # Create a separate DataFrame for global statistics
            if stats:
                stats_rows = []
                
                # Add Shannon statistics
                if 'shannon' in stats:
                    shannon_stats = stats['shannon']
                    stats_rows.append({
                        'statistic_type': 'Shannon Diversity',
                        'mean': shannon_stats.get('mean'),
                        'std': shannon_stats.get('std'),
                        'min': shannon_stats.get('min'),
                        'max': shannon_stats.get('max'),
                        'median': shannon_stats.get('median')
                    })
                
                # Add Simpson statistics
                if 'simpson' in stats:
                    simpson_stats = stats['simpson']
                    stats_rows.append({
                        'statistic_type': 'Simpson Diversity',
                        'mean': simpson_stats.get('mean'),
                        'std': simpson_stats.get('std'),
                        'min': simpson_stats.get('min'),
                        'max': simpson_stats.get('max'),
                        'median': simpson_stats.get('median')
                    })
                
                # Add Chao1 statistics
                if 'chao1' in stats:
                    chao1_stats = stats['chao1']
                    stats_rows.append({
                        'statistic_type': 'Chao1 Richness',
                        'mean': chao1_stats.get('mean'),
                        'std': chao1_stats.get('std'),
                        'min': chao1_stats.get('min'),
                        'max': chao1_stats.get('max'),
                        'median': chao1_stats.get('median')
                    })
                
                # Add species richness statistics
                if 'mean_species_per_sample' in stats:
                    stats_rows.append({
                        'statistic_type': 'Species Richness',
                        'mean': stats.get('mean_species_per_sample'),
                        'std': stats.get('std_species_per_sample'),
                        'min': stats.get('min_species_per_sample'),
                        'max': stats.get('max_species_per_sample'),
                        'median': None
                    })
                
                # Add total samples and projects
                stats_rows.append({
                    'statistic_type': 'Sample Count',
                    'value': stats.get('total_samples')
                })
                
                stats_rows.append({
                    'statistic_type': 'Project Count',
                    'value': stats.get('total_projects')
                })
                
                # Create stats DataFrame
                stats_df = pd.DataFrame(stats_rows)
            else:
                stats_df = pd.DataFrame()
            
            # Return both DataFrames
            return {
                'diversity_data': diversity_df,
                'summary_stats': stats_df
            }
            
        except Exception as e:
            print(f"Error preparing diversity CSV data: {str(e)}")
            import traceback
            traceback.print_exc()
            return {'diversity_data': pd.DataFrame(), 'summary_stats': pd.DataFrame()}
