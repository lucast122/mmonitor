import json
import random
import logging
import traceback

import dash
from dash import html, dcc  # Updated Dash imports
import dash_mantine_components as dmc
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
from django_plotly_dash import DjangoDash
import numpy as np
from functools import lru_cache
import time
import hashlib
from typing import Dict, List, Any, Optional, Union, Tuple, cast
from users.models import SequencingStatistics

# Initialize logger
logger = logging.getLogger(__name__)

class QC:
    def __init__(self, user_id):
        self.user_id = user_id
        self.app = DjangoDash('QC', add_bootstrap_links=True)
        
        # Get initial data
        self.stats_df, self.stats = self.get_data()
        
        # Initialize layout after data is loaded
        self._init_layout()

    def get_data(self):
        """Get data and update cache timestamp."""
        logger.info(f"Fetching QC data for user_id: {self.user_id}")
        # Only fetch the fields we actually need
        stats = SequencingStatistics.objects.filter(user_id=self.user_id).only(
            'sample_name', 'q20_score', 'q30_score', 'mean_read_length',
            'mean_quality_score', 'mean_gc_content', 'read_lengths',
            'avg_qualities', 'gc_contents_per_sequence'
        )
        
        if not stats.exists():
            logger.warning(f"No QC statistics found for user_id: {self.user_id}")
            return pd.DataFrame(), stats
        
        # Create a smaller DataFrame with only needed columns
        df = pd.DataFrame.from_records(
            stats.values('sample_name', 'q20_score', 'q30_score', 
                       'mean_read_length', 'mean_quality_score', 'mean_gc_content')
        )
        
        return df, stats

    def _init_layout(self):
        if self.stats:
            self.unique_sample_ids = self.stats.values('sample_name').distinct()
            self.unique_sample_ids = [item['sample_name'] for item in self.unique_sample_ids]
        else:
            self.unique_sample_ids = []

        qc_layout = self._create_qc_layout() if not self.stats_df.empty else self._create_empty_layout()
        self.app.layout = dmc.Container(
            qc_layout, 
            fluid=True, 
            size="100%", 
            style={
                "width": "100%",
                "paddingLeft": "0",
                "paddingRight": "0",
                "margin": 0
            }
        )

        # Register callbacks
        self._init_callbacks()

    def _init_callbacks(self):
        """Initialize callbacks for the QC app."""
        pass  # Callbacks are handled in the Index class

    def _create_qc_layout(self):
        qc_layout = dmc.Stack([
            dmc.Title("Quality Control Overview", order=2, style={"marginBottom": "20px"}),

            # Add explanation for quality scores
            dmc.Alert(
                title="About Quality Scores (Phred Scores)",
                children=[
                    "Quality scores (Q-scores) indicate the probability of a base calling error:",
                    html.Br(),
                    html.Ul([
                        html.Li([
                            "Q20 score: Shows the percentage of bases with a Phred score ≥ 20 ",
                            "(1 error in 100 bases, 99% accuracy)"
                        ]),
                        html.Li([
                            "Q30 score: Shows the percentage of bases with a Phred score ≥ 30 ",
                            "(1 error in 1,000 bases, 99.9% accuracy)"
                        ]),
                        html.Li("Higher percentages indicate better sequencing quality")
                    ]),
                    html.Br(),
                    "Understanding the Violin Plot:",
                    html.Ul([
                        html.Li("The width of the violin shape shows how many samples have a particular quality score, GC content, or mean read lengths"),
                        html.Li("Wider sections indicate more samples with that score"),
                        html.Li("The box inside shows the median (middle line) and quartiles (box edges)"),
                        html.Li("Individual points show exact values for each sample")
                    ])
                ],
                color="blue",
                style={"marginBottom": "20px"}
            ),

            # Overview plots
            dmc.Grid([
                dcc.Graph(id='quality-overview', figure=self._create_quality_overview()),
                dcc.Graph(id='read-length-overview', figure=self._create_read_length_overview()),
                dcc.Graph(id='gc-overview', figure=self._create_gc_overview())
            ], gutter="xl"),

            dmc.Space(h="xl")
        ])

        return qc_layout

    def _create_empty_layout(self):
        return dmc.Text("QC app could not be loaded as no QC statistics were found.")

    def _create_quality_overview(self):
        """Create overview of quality scores across all samples."""
        if not self.stats.exists():
            return go.Figure()

        fig = go.Figure()

        # Add violin plot for quality scores
        fig.add_trace(go.Violin(
            y=self.stats_df['mean_quality_score'],
            name='Quality Score',
            box_visible=True,
            meanline_visible=True,
            points='all',
            text=self.stats_df['sample_name'],  # Add sample names for hover
            hovertemplate='Sample: %{text}<br>Quality Score: %{y:.1f}<extra></extra>'
        ))

        fig.update_layout(
            title='Quality Score Distribution Across All Samples',
            yaxis_title='Quality Score',
            showlegend=False,
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            
            height=500,   # Added fixed height for better proportions
            margin=dict(l=30, r=30, t=40, b=40)  # Reduced margins
        )

        return fig

    def _create_read_length_overview(self):
        """Create overview of read lengths across all samples."""
        if not self.stats.exists():
            return go.Figure()

        fig = go.Figure()

        # Add violin plot for read lengths
        fig.add_trace(go.Violin(
            y=self.stats_df['mean_read_length'],
            name='Read Length',
            box_visible=True,
            meanline_visible=True,
            points='all',
            text=self.stats_df['sample_name'],  # Add sample names for hover
            hovertemplate='Sample: %{text}<br>Mean Length: %{y:.0f}bp<extra></extra>'
        ))

        fig.update_layout(
            title='Read Length Distribution Across All Samples',
            yaxis_title='Read Length (bp)',
            showlegend=False,
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            height=500,   # Added fixed height for better proportions
            margin=dict(l=30, r=30, t=40, b=40)  # Reduced margins
        )

        return fig

    def _create_gc_overview(self):
        """Create overview of GC content across all samples."""
        if not self.stats.exists():
            return go.Figure()

        fig = go.Figure()

        # Add violin plot for GC content
        fig.add_trace(go.Violin(
            y=self.stats_df['mean_gc_content'],
            name='GC Content',
            box_visible=True,
            meanline_visible=True,
            points='all',
            text=self.stats_df['sample_name'],  # Add sample names for hover
            hovertemplate='Sample: %{text}<br>GC Content: %{y:.1f}%<extra></extra>'
        ))

        fig.update_layout(
            title='GC Content Distribution Across All Samples',
            yaxis_title='GC Content (%)',
            showlegend=False,
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
            height=500,   # Added fixed height for better proportions
            margin=dict(l=30, r=30, t=40, b=40)  # Reduced margins
        )

        return fig

    def update_qc_plots(self, selected_sample):
        """Update QC plots with caching."""
        logger.info(f"Starting update_qc_plots for sample: {selected_sample}")
        
        if not selected_sample:
            logger.warning("No sample selected, returning empty figures")
            return [go.Figure() for _ in range(4)]

        try:
            # Get sample data
            logger.info(f"Fetching stats for sample: {selected_sample}")
            stats = self.stats.filter(sample_name=selected_sample).first()
            if not stats:
                logger.error(f"No QC data found for sample {selected_sample}")
                raise ValueError(f"No QC data found for sample {selected_sample}")

            # Log the data we're working with
            logger.info(f"Sample data retrieved - Q20: {stats.q20_score}, Q30: {stats.q30_score}")
            logger.info("Checking data distributions...")
            logger.info(f"Read lengths data: {bool(stats.read_lengths)}")
            logger.info(f"Quality scores data: {bool(stats.avg_qualities)}")
            logger.info(f"GC content data: {bool(stats.gc_contents_per_sequence)}")

            # Create cache key from essential data
            data_hash = hashlib.md5(
                f"{stats.read_lengths}{stats.avg_qualities}{stats.gc_contents_per_sequence}".encode()
            ).hexdigest()
            
            # Check cache status
            cache_info = self._generate_plots_cached.cache_info()
            logger.info(f"Cache info - hits: {cache_info.hits}, misses: {cache_info.misses}, size: {cache_info.currsize}")
            
            # Generate plots using cached function
            logger.info("Generating plots...")
            plots = self._generate_plots_cached(
                selected_sample,
                data_hash,
                stats.read_lengths,
                stats.avg_qualities,
                stats.gc_contents_per_sequence,
                stats.q20_score,
                stats.q30_score
            )
            
            logger.info(f"Successfully generated {len(plots)} plots")
            return plots
            
        except Exception as e:
            logger.error(f"Error in update_qc_plots: {str(e)}")
            logger.error(f"Full traceback: {traceback.format_exc()}")
            return [go.Figure() for _ in range(4)]

    @staticmethod
    def _parse_json_data(data_str: str, default=None) -> dict:
        """Safely parse JSON data with error handling."""
        if not data_str:
            return default or {}
        try:
            return json.loads(data_str)
        except Exception as e:
            logger.error(f"Error parsing JSON data: {str(e)}")
            return default or {}

    @staticmethod
    def _create_base_layout(title: str, xaxis_title: str, yaxis_title: str) -> dict:
        """Create base layout for plots."""
        return {
            'title': title,
            'xaxis_title': xaxis_title,
            'yaxis_title': yaxis_title,
            'showlegend': False,
            'template': 'plotly_white',
            'margin': {'l': 40, 'r': 40, 't': 40, 'b': 40}
        }

    @staticmethod
    @lru_cache(maxsize=32)
    def _generate_plots_cached(sample_name: str, data_hash: str, read_lengths: str, avg_qualities: str, 
                             gc_content: str, q20_score: float, q30_score: float) -> List[go.Figure]:
        """Generate cached plots for a sample."""
        logger.info(f"Cache miss for {sample_name} - generating new plots")
        
        try:
            # Parse all data at once
            read_length_data = QC._parse_json_data(read_lengths)
            quality_data = QC._parse_json_data(avg_qualities)
            gc_data = QC._parse_json_data(gc_content)

            # Initialize all figures at once with base layouts
            figures = [
                go.Figure(layout=QC._create_base_layout(
                    f'Read Length Distribution for {sample_name}',
                    'Read Length (bp)',
                    'Count'
                )),
                go.Figure(layout=QC._create_base_layout(
                    f'Quality Score Distribution for {sample_name}',
                    'Quality Score',
                    'Count'
                )),
                go.Figure(layout=QC._create_base_layout(
                    f'GC Content Distribution for {sample_name}',
                    'GC Content (%)',
                    'Count'
                )),
                go.Figure(layout=QC._create_base_layout(
                    f'Base Composition for {sample_name}',
                    'Position',
                    'Percentage'
                ))
            ]

            # Add data to read length plot if available
            if read_length_data:
                x_vals = list(map(float, read_length_data.keys())) if isinstance(read_length_data, dict) else list(range(len(read_length_data)))
                y_vals = list(map(float, read_length_data.values())) if isinstance(read_length_data, dict) else read_length_data
                figures[0].add_trace(go.Bar(
                    x=x_vals,
                    y=y_vals,
                    name='Read Length Distribution',
                    marker_color='rgb(67, 147, 195)',
                    hovertemplate='Length: %{x} bp<br>Count: %{y}<extra></extra>'
                ))

            # Add data to quality score plot if available
            if quality_data:
                x_vals = list(map(float, quality_data.keys())) if isinstance(quality_data, dict) else list(range(len(quality_data)))
                y_vals = list(map(float, quality_data.values())) if isinstance(quality_data, dict) else quality_data
                figures[1].add_trace(go.Bar(
                    x=x_vals,
                    y=y_vals,
                    name='Quality Score Distribution',
                    marker_color='rgb(44, 160, 44)',
                    hovertemplate='Score: %{x}<br>Count: %{y}<extra></extra>'
                ))
                # Add Q20 and Q30 lines efficiently
                if q20_score or q30_score:
                    shapes = []
                    annotations = []
                    if q20_score:
                        shapes.append({
                            'type': 'line',
                            'x0': 20, 'x1': 20,
                            'y0': 0, 'y1': 1,
                            'yref': 'paper',
                            'line': {'dash': 'dash', 'color': 'red'}
                        })
                        annotations.append({
                            'x': 20,
                            'y': 1,
                            'text': f"Q20: {q20_score:.1f}%",
                            'showarrow': False,
                            'yref': 'paper'
                        })
                    if q30_score:
                        shapes.append({
                            'type': 'line',
                            'x0': 30, 'x1': 30,
                            'y0': 0, 'y1': 1,
                            'yref': 'paper',
                            'line': {'dash': 'dash', 'color': 'green'}
                        })
                        annotations.append({
                            'x': 30,
                            'y': 1,
                            'text': f"Q30: {q30_score:.1f}%",
                            'showarrow': False,
                            'yref': 'paper'
                        })
                    figures[1].update_layout(shapes=shapes, annotations=annotations)

            # Add data to GC content plot if available
            if gc_data:
                x_vals = list(map(float, gc_data.keys())) if isinstance(gc_data, dict) else list(range(len(gc_data)))
                y_vals = list(map(float, gc_data.values())) if isinstance(gc_data, dict) else gc_data
                figures[2].add_trace(go.Bar(
                    x=x_vals,
                    y=y_vals,
                    name='GC Content Distribution',
                    marker_color='rgb(214, 39, 40)',
                    hovertemplate='GC Content: %{x}%<br>Count: %{y}<extra></extra>'
                ))

            return figures

        except Exception as e:
            logger.error(f"Error in _generate_plots_cached: {str(e)}")
            logger.error(f"Full traceback: {traceback.format_exc()}")
            return [go.Figure() for _ in range(4)]

    def update_individual_sample_plots(self, selected_sample):
        """Update individual sample plots."""
        # Simply delegate to update_qc_plots since it handles all the caching
        return self.update_qc_plots(selected_sample)

    def _create_read_length_plot(self, data):
        """Create read length distribution plot."""
        fig = go.Figure()
        
        # Create histogram with custom styling
        fig.add_trace(go.Histogram(
            x=data['read_length'],
            y=data['count'],
            histfunc='sum',
            name='Read Length Distribution',
            marker=dict(
                color='rgb(67, 147, 195)',
                line=dict(color='white', width=1)
            ),
            hovertemplate='Read Length: %{x:,.0f}bp<br>Count: %{y:,.0f}<extra></extra>'
        ))

        fig.update_layout(
            title=dict(
                text='Read Length Distribution',
                x=0.5,
                font=dict(size=20)
            ),
            xaxis_title='Read Length (bp)',
            yaxis_title='Count',
            showlegend=False,
            plot_bgcolor='white',
            paper_bgcolor='white',
            xaxis=dict(
                gridcolor='lightgray',
                showline=True,
                linewidth=1,
                linecolor='black',
                mirror=True
            ),
            yaxis=dict(
                gridcolor='lightgray',
                showline=True,
                linewidth=1,
                linecolor='black',
                mirror=True
            ),
            margin=dict(l=50, r=50, t=80, b=50)
        )
        return fig

    def _create_quality_plot(self, data):
        """Create quality score distribution plot."""
        fig = go.Figure()
        
        # Create histogram with custom styling
        fig.add_trace(go.Histogram(
            x=data['mean_quality'],
            y=data['count'],
            histfunc='sum',
            name='Quality Distribution',
            marker=dict(
                color='rgb(44, 160, 44)',
                line=dict(color='white', width=1)
            ),
            hovertemplate='Quality Score: %{x:.1f}<br>Count: %{y:,.0f}<extra></extra>'
        ))

        # Add vertical lines for Q20 and Q30
        fig.add_vline(x=20, line_dash="dash", line_color="red", annotation_text="Q20")
        fig.add_vline(x=30, line_dash="dash", line_color="green", annotation_text="Q30")

        fig.update_layout(
            title=dict(
                text='Quality Score Distribution',
                x=0.5,
                font=dict(size=20)
            ),
            xaxis_title='Mean Quality Score (Phred)',
            yaxis_title='Count',
            showlegend=False,
            plot_bgcolor='white',
            paper_bgcolor='white',
            xaxis=dict(
                gridcolor='lightgray',
                showline=True,
                linewidth=1,
                linecolor='black',
                mirror=True,
                range=[0, 40]  # Typical range for quality scores
            ),
            yaxis=dict(
                gridcolor='lightgray',
                showline=True,
                linewidth=1,
                linecolor='black',
                mirror=True
            ),
            margin=dict(l=50, r=50, t=80, b=50)
        )
        
        # make fig wider
        fig.update_layout(
            width=800
        )
        return fig

    def _create_gc_plot(self, data):
        """Create GC content distribution plot."""
        fig = go.Figure()
        
        # Create histogram with custom styling
        fig.add_trace(go.Histogram(
            x=data['gc_content'],
            y=data['count'],
            histfunc='sum',
            name='GC Content Distribution',
            marker=dict(
                color='rgb(214, 39, 40)',
                line=dict(color='white', width=1)
            ),
            hovertemplate='GC Content: %{x:.1f}%<br>Count: %{y:,.0f}<extra></extra>'
        ))

        fig.update_layout(
            title=dict(
                text='GC Content Distribution',
                x=0.5,
                font=dict(size=20)
            ),
            xaxis_title='GC Content (%)',
            yaxis_title='Count',
            showlegend=False,
            plot_bgcolor='white',
            paper_bgcolor='white',
            xaxis=dict(
                gridcolor='lightgray',
                showline=True,
                linewidth=1,
                linecolor='black',
                mirror=True,
                range=[20, 80]  # Typical range for GC content
            ),
            yaxis=dict(
                gridcolor='lightgray',
                showline=True,
                linewidth=1,
                linecolor='black',
                mirror=True
            ),
            margin=dict(l=50, r=50, t=80, b=50)
        )
        return fig

    def _create_base_composition_plot(self, data):
        """Create base composition plot."""
        fig = go.Figure()

        colors = {
            'A': 'rgb(31, 119, 180)',  # Blue
            'C': 'rgb(44, 160, 44)',   # Green
            'G': 'rgb(214, 39, 40)',   # Red
            'T': 'rgb(148, 103, 189)'  # Purple
        }

        for base in ['A', 'C', 'G', 'T']:
            base_data = data[data['base'] == base]
            fig.add_trace(go.Scatter(
                x=base_data['position'],
                y=base_data['count'],
                name=base,
                mode='lines',
                line=dict(
                    color=colors[base],
                    width=2
                ),
                hovertemplate='Position: %{x}<br>' + f'{base}: ' + '%{y:,.0f}<extra></extra>'
            ))

        fig.update_layout(
            title=dict(
                text='Base Composition Along Reads',
                x=0.5,
                font=dict(size=20)
            ),
            xaxis_title='Position in Read',
            yaxis_title='Count',
            showlegend=True,
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.99,
                bgcolor='rgba(255,255,255,0.8)'
            ),
            plot_bgcolor='white',
            paper_bgcolor='white',
            xaxis=dict(
                gridcolor='lightgray',
                showline=True,
                linewidth=1,
                linecolor='black',
                mirror=True
            ),
            yaxis=dict(
                gridcolor='lightgray',
                showline=True,
                linewidth=1,
                linecolor='black',
                mirror=True
            ),
            margin=dict(l=50, r=50, t=80, b=50)
        )
        return fig