import base64
from typing import Tuple, List, Any
from functools import lru_cache

import dash_bootstrap_components as dbc
from dash import dcc
import dash_mantine_components as dmc
import pandas as pd
import plotly.express as px
from dash import dash_table
from dash import html
from dash_ag_grid import AgGrid
from dash_iconify import DashIconify
from django_plotly_dash import DjangoDash
from plotly import graph_objects as go
from users.models import NanoporeRecord
from users.models import SequencingStatistics
import numpy as np
from natsort import natsorted
import requests
import json
from typing import Dict, List
import os


class Taxonomy:
    """
    App to display the abundances of taxonomies in various formats.
    """

    def __init__(self, user_id):
        print("Initializing Taxonomy app")
        self.text_style = 'text-primary my-2'
        self.user_id = user_id
        app_name = 'taxonomy'
        try:
            # Initialize the Dash app
            self.app = DjangoDash(
                app_name,
                add_bootstrap_links=False,
                serve_locally=True  
            )
            
        except Exception as e:
            logger.error(f"Error initializing {app_name} app: {str(e)}", exc_info=True)
            # Create a minimal app to avoid breaking the whole dashboard
            self.app = DjangoDash(app_name, suppress_callback_exceptions=True)
            self.app.layout = html.Div(f"Error loading {app_name} app")
            
            # Initialize the app layout
            self.app.layout = self.create_layout()
        
        self._cache = {}
        
        # Colors initialization
        self.colors = [
            '#a9a9a9', '#2f4f4f', '#556b2f', '#8b4513', '#6b8e23', '#006400', '#708090',
            '#8b0000', '#3cb371', '#bc8f8f', '#663399', '#008080', '#bdb76b', '#4682b4',
            '#000080', '#9acd32', '#20b2aa', '#cd5c5c', '#32cd32', '#daa520', '#8fbc8f',
            '#8b008b', '#b03060', '#d2b48c', '#ff0000', '#ff8c00', '#ffd700', '#ffff00',
            '#c71585', '#00ff00', '#ba55d3', '#00ff7f', '#4169e1', '#e9967a', '#dc143c',
            '#00ffff', '#00bfff', '#f4a460', '#0000ff', '#a020f0', '#adff2f', '#ff7f50',
            '#ff00ff', '#db7093', '#eee8aa', '#6495ed', '#dda0dd', '#87ceeb', '#ff1493',
            '#f5f5dc', '#afeeee', '#ee82ee', '#98fb98', '#7fffd4', '#e6e6fa', '#ffc0cb'
        ]
        
        # Initialize data first
        self.records = self._get_cached_records()
        self.df = self._init_data()
        
        # Initialize derived data (always call to set defaults)
        self._init_derived_data()
        
        self._init_layout()

    def _get_cached_records(self):
        """Get records from database with caching"""
        if 'records' in self._cache:
            return self._cache['records']
            
        records = NanoporeRecord.objects.filter(user_id=self.user_id).select_related()
        self._cache['records'] = records
        return records

    def _init_data(self):
        """Initialize the main dataframe"""
        if 'df' in self._cache:
            return self._cache['df']
            
        if not self.records.exists():
            print("No records found")
            return pd.DataFrame()
            
        df = pd.DataFrame.from_records(self.records.values())
        df = df.drop_duplicates(subset=['sample_id', 'taxonomy'])
        
        # Clean data
        df.loc[df['taxonomy'].isin(['Not available', 'Not Available']), 'taxonomy'] = np.nan
        df.dropna(inplace=True)
        # Only use columns that actually exist in the mmonitor table
        drop_columns = ['sample_id', 'abundance', 'taxonomy', 'project_id', 'user_id']
        # Add optional columns if they exist
        if 'subproject' in df.columns:
            drop_columns.append('subproject')
        if 'date' in df.columns:
            drop_columns.append('date')
            
        df.drop_duplicates(
            subset=drop_columns,
            inplace=True)
        
        # Sort by available columns
        sort_columns = ['project_id', 'sample_id']
        if 'subproject' in df.columns:
            sort_columns.insert(1, 'subproject')  # Insert between project_id and sample_id
            
        self._cache['df'] = df.sort_values(by=sort_columns)
        return self._cache['df']

    def _init_derived_data(self):
        """Initialize all derived data after main dataframe is loaded"""
        # Initialize default empty values first
        self.unique_sample_ids = []
        self.unique_samples = []
        self.unique_projects_ids = []
        self.unique_subprojects = []
        self.unique_dates = []
        
        # If we have data, populate with actual values
        if not self.df.empty and self.records.exists():
            self.unique_sample_ids = self.df['sample_id'].unique().tolist()
            self.unique_samples = self.records.values('sample_id').distinct()
            self.unique_projects_ids = [str(item['project_id']) for item in self.records.values('project_id').distinct()]
            
            # These fields don't exist in the database - they're properties in the model
            # So we get them from the DataFrame instead
            if 'subproject' in self.df.columns:
                self.unique_subprojects = self.df['subproject'].unique().tolist()
            else:
                # Since subproject is a property that returns 'default', use that
                self.unique_subprojects = ['default']
            
            if 'date' in self.df.columns:
                self.unique_dates = [str(d) for d in self.df['date'].unique()]
            else:
                # Since date is a property that returns today's date, use that
                from datetime import date
                self.unique_dates = [str(date.today())]
            
            # Only access dataframe columns if data exists
            self.unique_species = self.df['taxonomy'].unique()
            
            # These taxonomy fields are properties in the model, not database fields
            # Since they all return 'Unknown', we'll use that
            self.unique_genera = ['Unknown']
            self.unique_families = ['Unknown']
            self.unique_classes = ['Unknown']
            self.unique_orders = ['Unknown']
            self.unique_phyla = ['Unknown']
        else:
            # Set empty arrays for taxonomy data when no data exists
            self.unique_species = []
            self.unique_genera = []
            self.unique_families = []
            self.unique_classes = []
            self.unique_orders = []
            self.unique_phyla = []
        
        # Only compute aggregations if we have data
        if not self.df.empty:
            self.df_full_for_diversity = self.df.pivot_table(
                index='sample_id',
                columns='taxonomy',
                values='abundance',
                fill_value=0
            )
        else:
            self.df_full_for_diversity = pd.DataFrame()
        
        # Project mapping - only if we have records
        if self.records.exists():
            sample_project_mapping = self.records.values('sample_id', 'project_id')
            self.sample_to_project_dict = {item['sample_id']: item['project_id'] for item in sample_project_mapping}
        else:
            self.sample_to_project_dict = {}
        
        # Sort and get counts - only if we have data
        if not self.df.empty:
            self.df_sorted = self.df.sort_values(by=["sample_id", "abundance"], ascending=[True, False])
            self.unique_counts = min(self.df.nunique()[1], 100)
        else:
            self.df_sorted = pd.DataFrame()
            self.unique_counts = 0
        
        # Color assignments - only if we have data
        if not self.df_sorted.empty:
            unique_species = self.df_sorted['taxonomy'].unique()
            self.species_colors = self.get_distinct_colors(len(unique_species))
            self.species_colors = {species: color for species, color in zip(unique_species, self.species_colors)}
        else:
            self.species_colors = {}
        
        if not self.df_sorted.empty:
            # Since tax_genus and tax_family don't exist as columns, just use empty dicts
            self.genus_colors = {}
            self.family_colors = {}
            self.combined_color_dict = {**self.species_colors}
        else:
            self.genus_colors = {}
            self.family_colors = {}
            self.combined_color_dict = {}

        self.legend_info = []
        self.abundance_lists = None
        self.shannon_diversity = None
        self.simpson_diversity = None

    def _init_layout(self):
        print("Unique dates:", self.unique_dates)
        
        
        app_explanation = dmc.Alert(
                    children=[
                        "Analyze and visualize the taxonomic composition of your samples:",
                        html.Ul([
                            html.Li("Compare abundance distributions across different samples"),
                            html.Li("Visualization options:"),
                            html.Ul([
                                html.Li("Change plot type to see different visualizations"),
                                html.Li("Change taxonomic rank to see different taxonomic levels"),
                                html.Li("Download high quality svg for publication or export data for custom plotting"),
                                html.Li("If samples share same dates filter by Project or Subproject before using date for plotting"),
                                
                            ]),
                            html.Li("Filter by project, date, or specific samples"),
                        ], style={"marginBottom": "0", "marginTop": "0", "paddingTop": "0"})
                    ],
                    color="blue",
                    variant="light",
                    style={"marginBottom": "5px", "padding": "5px 10px"}
                )
        sample_select_dropdown = dmc.MultiSelect(
            id='sample-select-value-tax',
            data=[{'value': sample, 'label': sample} for sample in self.unique_sample_ids],
            label="Samples to plot",
            placeholder="Select samples...",
            searchable=True,
            clearable=True,
            style={'width': '100%', 'maxHeight': '200px', 'overflowY': 'auto'},
            value=self.unique_sample_ids,
            maxDropdownHeight=200
        )

        project_dropdown = dmc.Select(
            id='project-dropdown-tax',
            data=[{'value': 'ALL', 'label': 'All Projects'}] + [{'value': project, 'label': project} for project in self.unique_projects_ids],
            label="Select by project",
            placeholder="Select a project...",
            clearable=True,
            value='ALL'
        )

        subproject_dropdown = dmc.Select(
            id='subproject-dropdown-tax',
            data=[{'value': 'ALL', 'label': 'All Subprojects'}] + [{'value': subproject, 'label': subproject} for subproject in self.unique_subprojects],
            label="Select by subproject",
            placeholder="Select a subproject...",
            clearable=True,
            value='ALL'
        )

        date_dropdown = dmc.Select(
            id='date-dropdown-tax',
            data=[{'value': 'ALL', 'label': 'All Dates'}] + [{'value': str(date), 'label': str(date)} for date in self.unique_dates],
            label="Select Samples by date",
            placeholder="Select date...",
            clearable=True,
            value='ALL'
        )

        plot_type_dropdown = dmc.Select(
            id='plot-type-dropdown-tax',
            data=[{'value': 'stackedbar', 'label': 'Stacked Barchart'},
                  {'value': 'groupedbar', 'label': 'Grouped Barchart'},
                  {'value': 'area', 'label': 'Area plot'},
                  {'value': 'line', 'label': 'Line plot'},
                  {'value': 'heatmap', 'label': 'Heatmap'},
                  {'value': 'scatter', 'label': 'Scatter plot'}],
            label="Plot type",
            value='stackedbar',
            placeholder="Select plot type...",
            clearable=True,
        )

        tax_rank_dropdown = dmc.Select(
            id='tax-rank-dropdown-tax',
            data=[
                {'value': 'taxonomy', 'label': 'Species'},
                # Note: Other taxonomy levels are not available in the current data structure
            ],
            label="Taxonomic rank for plot",
            value='taxonomy',
            placeholder="Select taxonomic rank...",
            clearable=True,
        )

        use_date_checkbox = dmc.Switch(
            id='use-date-checkbox-tax',
            label='Use Date for plotting',
            checked=False,
            size="sm"
        )

        legend_toggle_checkbox = dmc.Switch(
            id='legend-toggle-checkbox-tax',
            label='Toggle legend',
            checked=True,
            size="sm"
        )

        show_legend_modal_button = dmc.Button(
            "Show Legend in Window",
            id="show-legend-modal-button",
            variant="outline",
            size="sm"
        )

        download_legend_button = dmc.Button(
            "Download Legend",
            id="download-legend-button",
            variant="outline",
            size="sm"
        )

        # Get the total number of unique taxa for the default legend items value
        if not self.df.empty:
            total_taxa = len(self.df[self.df['abundance'] > 0]['taxonomy'].unique())
        else:
            total_taxa = 0

        download_component = dcc.Download(id="download-plot")

        export_controls = dmc.Paper(
            children=[
                dmc.Text("Export Data", fw=500, size="sm", style={"marginBottom": "10px"}),
                dmc.Group([
                    dmc.NumberInput(
                        id="width-input",
                        label="Plot Width",
                        value=1200,
                        min=100,
                        max=3000,
                        step=100,
                        size="sm",
                        style={"width": 100}
                    ),
                    dmc.NumberInput(
                        id="height-input",
                        label="Plot Height",
                        value=800,
                        min=100,
                        max=3000,
                        step=100,
                        size="sm",
                        style={"width": 100}
                    ),
                    dmc.NumberInput(
                        id="legend-items-input",
                        label="Legend Items",
                        value=total_taxa,
                        min=1,
                        max=total_taxa,
                        step=1,
                        size="sm",
                        style={"width": 100}
                    ),
                    dmc.Button(
                        "Export as SVG",
                        id="export-svg-button",
                        variant="outline",
                        size="sm"
                    ),
                    dmc.Group([
            dmc.Button('Export All', leftSection=DashIconify(icon="foundation:page-export-csv"), id="btn-download-csv-taxonomy", variant="outline", size="sm"),
            dmc.Button('Export Counts', leftSection=DashIconify(icon="foundation:page-export-csv"), id="btn-download-counts-taxonomy", variant="outline", size="sm"),
        ])
                ], gap="sm"),
                download_component
            ],
            p="md",
            shadow="sm",
            withBorder=True,
            style={'width': '100%', 'marginTop': '20px'}
        )
        
        

        # Move the graph container definition here
        graph_container = dmc.Paper(
            children=[
                dcc.Graph(
                    id="graph1",
                    figure=self.plot_stacked_bar(self.df_sorted, False, 'taxonomy'),
                    config={"displayModeBar": True, "displaylogo": False}
                )
            ],
            p="md",
            shadow="sm",
            withBorder=True,
            style={"width": "100%",  'height': '600px'},
        )

        # slider = dmc.Slider(
        #     id='slider-tax',
        #     min=1,
        #     max=self.unique_counts,
        #     value=self.unique_counts,
        #     marks=[{"value": i, "label": str(i)} for i in range(1, self.unique_counts + 1, max(1, self.unique_counts // 5))],
        #     step=1,
        #     labelAlwaysOn=False,
        #     style={"width": "100%"}
        # )


        download_components = html.Div([
            dcc.Download(id="download-csv"),
            dcc.Download(id="download-svg-taxonomy"),
            dcc.Download(id="download-counts"),
            dcc.Download(id="download-legend-svg"),
            download_component,
        ])

        legend_modal = dmc.Modal(
            title="Legend",
            id="legend-modal",
            size="70%",
            children=[html.Div(id="legend-content")]
        )


        layout = dmc.Stack([
            dmc.Title("Taxonomy Analysis", order=2),
            
            dmc.Accordion(
    children=[
        dmc.AccordionItem(
            [
                dmc.AccordionControl("About Taxonomy Analysis"),
                dmc.AccordionPanel(
                    app_explanation
                ),
            ],
                    value="customization",
                )
                
            ]),

            sample_select_dropdown,
            dmc.SimpleGrid(
                cols=3,
                spacing="md",
                children=[
                    project_dropdown,
                    subproject_dropdown,
                    date_dropdown,
                    plot_type_dropdown,
                    tax_rank_dropdown,
                    dmc.Group([use_date_checkbox, legend_toggle_checkbox]),
                ]
            ),
            dmc.Space(h=20),
            graph_container,
            export_controls,
            dmc.Space(h=20),
            # slider,
            dmc.Space(h=20),
            
            
            dmc.Space(h=20),
            dmc.Space(h=20),
            download_components,
            legend_modal,
            html.Div(id='taxonomy-error-message', style={'color': 'red'}),
        ])

        container = dmc.Container(
            layout,
            fluid=True,
            style={"width": "calc(100% - 20px)",  'maxWidth': '100%', 'padding': '0px', 'margin': '0px'}
        )

        self.app.layout = container

    def _generate_table_ag_grid(self, max_rows=40):
        """
        Generate data to populate a dash-ag-grid table with.
        dash-ag-grid requires the data in a list of dictionaries format and a collection of column definitions.

        Args:
        max_rows (int, optional): The maximum number of rows to include in the table. Defaults to 40.

        Returns:
        Dash AgGrid component with the specified data and columns.
        """
        # Limit the number of rows if needed
        if max_rows is not None:
            df_display = self.df.head(max_rows)
        else:
            df_display = self.df

        # Generate column definitions for dash-ag-grid
        column_defs = [{"headerName": col, "field": col, "sortable": True, "filter": True} for col in
                       df_display.columns]

        # Create the AgGrid component
        ag_grid_table = AgGrid(

            id="table-correlations",
            rowData=df_display.to_dict('records'),
            className="ag-theme-balham",

            columnDefs=column_defs)
        return ag_grid_table, [{'name': i, 'id': i} for i in self.df.columns]

    def _generate_table_data_cols(self, max_rows=40) -> Tuple[List[Any], List[Any]]:
        """
        Generate data to populate a dash data table with.
        A dash data table requires the data in a dict format as well as a collection of mapped column names.

        Args:
        max_rows (int, optional): The maximum number of rows to include in the table. Defaults to 40.

        Returns:
        Tuple[List[Any], List[Any]]: A tuple containing the data for the table in dict format and a collection of mapped column names.
        """
        # q = f"SELECT * FROM nanopore LIMIT {max_rows}"
        # self.df = self.df = pd.read_sql_query(q, self._engine)
        return self.df.to_dict('records'), [{'name': i, 'id': i} for i in self.df.columns]


    
    def plot_stacked_bar(self, df, use_date, taxonomic_rank):
        x_axis = "date" if use_date else "sample_id"
        
        # Sort the sample_ids alphanumerically
        sorted_samples = natsorted(df[x_axis].unique())
        
        # Create a categorical type with the sorted order
        df[x_axis] = pd.Categorical(df[x_axis], categories=sorted_samples, ordered=True)
        
        # Sort the dataframe
        df = df.sort_values(x_axis)

        total_abundance = df.groupby(taxonomic_rank)['abundance'].sum().sort_values(ascending=False)
        
        fig = px.bar(df, x=x_axis, y="abundance", color=taxonomic_rank, barmode="stack",
                     category_orders={x_axis: sorted_samples},
                     color_discrete_map=self.combined_color_dict)
        
        fig.update_layout(legend={'traceorder': 'normal'}, xaxis_title="Sample ID" if not use_date else "Date", yaxis_title="Relative Abundance")
        fig.data = sorted(fig.data, key=lambda trace: total_abundance[trace.name], reverse=True)
        
        return self.apply_common_layout(fig)

    def plot_grouped_bar(self, df, use_date, taxonomic_rank):
        x_axis = "date" if use_date else "sample_id"
        
        sorted_samples = natsorted(df[x_axis].unique())
        df[x_axis] = pd.Categorical(df[x_axis], categories=sorted_samples, ordered=True)
        df = df.sort_values(x_axis)

        fig = px.bar(df, x=x_axis, y="abundance", color=taxonomic_rank, barmode="group",
                     category_orders={x_axis: sorted_samples},
                     color_discrete_map=self.combined_color_dict)
        
        return self.apply_common_layout(fig)

    def plot_area(self, df, use_date, taxonomic_rank):
        x_axis = "date" if use_date else "sample_id"
        
        sorted_samples = natsorted(df[x_axis].unique())
        df[x_axis] = pd.Categorical(df[x_axis], categories=sorted_samples, ordered=True)
        df = df.sort_values(x_axis)

        fig = px.area(df, x=x_axis, y="abundance", color=taxonomic_rank,
                      category_orders={x_axis: sorted_samples},
                      color_discrete_map=self.combined_color_dict)
        
        return self.apply_common_layout(fig)

    def plot_line(self, df, use_date, taxonomic_rank):
        x_axis = "date" if use_date else "sample_id"
        
        sorted_samples = natsorted(df[x_axis].unique())
        df[x_axis] = pd.Categorical(df[x_axis], categories=sorted_samples, ordered=True)
        df = df.sort_values(x_axis)

        fig = px.line(df, x=x_axis, y="abundance", color=taxonomic_rank,
                      category_orders={x_axis: sorted_samples},
                      color_discrete_map=self.combined_color_dict)
        
        return self.apply_common_layout(fig)

    def plot_heatmap(self, df, use_date, taxonomic_rank):
        x_axis = "date" if use_date else "sample_id"
        
        sorted_samples = natsorted(df[x_axis].unique())
        df[x_axis] = pd.Categorical(df[x_axis], categories=sorted_samples, ordered=True)
        df = df.sort_values(x_axis)

        pivot_df = df.pivot(index=taxonomic_rank, columns=x_axis, values='abundance')
        pivot_df = pivot_df.reindex(columns=sorted_samples)
        
        fig = px.imshow(pivot_df, aspect="auto",
                        color_continuous_scale='Viridis')
        
        return self.apply_common_layout(fig)

    def plot_scatter(self, df, use_date, taxonomic_rank):
        x_axis = "date" if use_date else "sample_id"
        
        sorted_samples = natsorted(df[x_axis].unique())
        df[x_axis] = pd.Categorical(df[x_axis], categories=sorted_samples, ordered=True)
        df = df.sort_values(x_axis)

        fig = px.scatter(df, x=x_axis, y="abundance", color=taxonomic_rank,
                         category_orders={x_axis: sorted_samples},
                         color_discrete_map=self.combined_color_dict)
        
        return self.apply_common_layout(fig)

    def plot_scatter_3d(self, df, taxonomic_rank):
        fig = px.scatter_3d(df, x='taxonomy', y='abundance', z='sample_id', color=taxonomic_rank)
        return self.apply_common_layout(fig)

    def plot_pie(self, df, sample_value_piechart, taxonomic_rank):
        filtered_df = df[df["sample_id"] == sample_value_piechart]
        aggregated_df = filtered_df.groupby(taxonomic_rank)['abundance'].sum().reset_index()
        fig = px.pie(aggregated_df, values='abundance', names=taxonomic_rank)
        # Customize hover template if necessary
        # fig.update_traces(
        #     hovertemplate="<br>".join([
        #                                   f"{x_axis}: %{{x}}",
        #                                   f"Abundance: %{{y}}",
        #                               ] + [f"{col}: %{{customdata[{i}]}}" for i, col in enumerate(df.columns) if
        #                                    col not in [x_axis, 'abundance']])
        # )

        # Add lines for project_id and subproject_id changes as before

        return fig

    def plot_scatter_3d(self, df, taxonomic_rank):
        # Plotting code for scatter 3D goes here...
        fig1 = px.scatter_3d(df, x='taxonomy', y='abundance', z='sample_id', color=taxonomic_rank)
        # fig2 = px.scatter_3d(df, x='abundance', y='taxonomy', z='sample_id', color='abundance',hover_data=self.hover_data)
        return fig1

    def plot_pie(self, df, sample_value_piechart, taxonomic_rank):
        filtered_df = df[df["sample_id"] == sample_value_piechart]
        aggregated_df = filtered_df.groupby(taxonomic_rank)['abundance'].sum().reset_index()
        fig1 = px.pie(aggregated_df, values='abundance', names=taxonomic_rank,
                      title=f'Pie chart of bioreactor taxonomy of sample {sample_value_piechart}')
        piechart_style = {'display': 'block'}
        return fig1, piechart_style

    # def split_alphanumeric(self, text):
    #     matches = re.findall(r'(\d+|\D+)', text)
    #     numbers = [int(m) for m in matches if m.isdigit()]
    #     non_numbers = [m for m in matches if not m.isdigit()]
    #     number = numbers[0] if numbers else float('inf')
    #     non_number = non_numbers[0] if non_numbers else ''
    #     return (number, non_number)


    def calculate_normalized_counts(self):
        counts_df = pd.DataFrame.from_records(NanoporeRecord.objects.filter(user_id=self.user_id).values())
        stats_df = pd.DataFrame.from_records(SequencingStatistics.objects.filter(user_id=self.user_id).values())

        # Merge the counts and stats dataframes on sample_id
        merged_df = pd.merge(counts_df, stats_df, left_on='sample_id', right_on="sample_name")

        # Calculate normalized counts by dividing each count by the numebr of bases sequenced and then multiplying by
        # 1 million to make normalized counts more similar to normal counts in value range
        merged_df['normalized_count'] = (merged_df['count'] / merged_df['total_bases_sequenced']) * 10000000

        # Selecting only necessary columns for the final DataFrame
        normalized_counts_df = merged_df[['sample_id', 'taxonomy', 'abundance', 'count', 'normalized_count']]
        return normalized_counts_df
    
    def apply_common_layout(self, fig):
        fig.update_layout(
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
        font=dict(family="Arial, sans-serif", size=16, color="black"),
        margin=dict(l=0, r=0, t=0, b=0),
        autosize=True,
        height=550,
        xaxis=dict(
            showline=True,
            showgrid=False,
            showticklabels=True,
            linecolor='rgb(204, 204, 204)',
            linewidth=2,
            ticks='outside',
            tickfont=dict(
                family='Arial',
                size=16,
                color='rgb(82, 82, 82)',
            ),
        ),
        yaxis=dict(
            showline=True,
            showgrid=False,
            showticklabels=True,
            linecolor='rgb(204, 204, 204)',
            linewidth=2,
            ticks='outside',
            tickfont=dict(
                family='Arial',
                size=18,
                color='rgb(82, 82, 82)',
            ),
    ))
        return fig

    def assign_color_to_taxonomy(self, df, taxonomy, color_dict):
        taxonomy_colors = {}
        for item in df[taxonomy].unique():
            species_in_taxonomy = df[df[taxonomy] == item]['taxonomy'].iloc[0]
            taxonomy_colors[item] = color_dict.get(species_in_taxonomy, self.colors[
                len(taxonomy_colors) % len(self.colors)])
        return taxonomy_colors

    def get_distinct_colors(self, n_colors):
        if n_colors <= len(self.colors):
            return self.colors[:n_colors]

        colors = list(self.colors)
        while len(colors) < n_colors:
            h = len(colors) * 137.508
            s = 0.75 + np.random.random() * 0.25
            v = 0.75 + np.random.random() * 0.25

            h = h % 360
            h_i = int(h / 60)
            f = h / 60 - h_i
            p = v * (1 - s)
            q = v * (1 - f * s)
            t = v * (1 - (1 - f) * s)

            if h_i == 0:
                r, g, b = v, t, p
            elif h_i == 1:
                r, g, b = q, v, p
            elif h_i == 2:
                r, g, b = p, v, t
            elif h_i == 3:
                r, g, b = p, q, v
            elif h_i == 4:
                r, g, b = t, p, v
            else:
                r, g, b = v, p, q

            new_color = "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))
            colors.append(new_color)

        return colors[:n_colors]
