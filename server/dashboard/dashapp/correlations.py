from typing import Any, List, Union

import dash_bootstrap_components as dbc
from dash import dcc
from dash import html
import dash_mantine_components as dmc
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from dash_iconify import DashIconify
from dash import dash_table
from django_plotly_dash import DjangoDash
import traceback
from users.models import NanoporeRecord, Metadata
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import io
import base64
import scipy.stats as stats


class Correlations:
    """
    App to display correlations between abundance and metadata and their probabilities.
    The graph and data table are independent.
    The data table content may be downloaded in csv format.

    WARNING: - correlations can only be calculated for complete and sanitized data
             - bootstrapping is slow
    """

    def __init__(self, user_id):
        print("Initializing Correlations app")
        self.text_style = 'text-primary my-2'
        self.user_id = user_id
        self.app = DjangoDash('correlations', add_bootstrap_links=True)
        self.correlation_matrix = None
        self.pvalue_matrix = None
        self.corr_method = 'pearson'  # Default correlation method

        print("Getting data")
        self.nanopore_df = self.get_data()
        print("Getting metadata")
        self.meta_df = self.get_all_meta()
        print("Computing correlations")
        self.correlation_scores = self.compute_correlations_for_taxonomies(self.nanopore_df, self.meta_df,
                                                                           self.corr_method, 'taxonomy')

        self._methods = [
            'Pearson',
            'Spearman',
            'Kendall',
        ]
        self._tests = [
            'T-Test',
            'Bootstrap',
            'Fisher Z Transform',
        ]
        self._table_df = None

        print("Initializing layout")
        self._init_layout()
        # Callbacks are now in index.py
        print("Correlations app initialization complete")

    def get_data(self):
        print("Fetching NanoporeRecord data")
        records = NanoporeRecord.objects.filter(user_id=self.user_id)
        self.records = records

        if not records.exists():
            print("No records found")
            return pd.DataFrame()  # Return an empty DataFrame if no records are found
        df = pd.DataFrame.from_records(records.values())
        self._taxonomies = df['taxonomy'].unique()
        print(f"Found {len(self._taxonomies)} unique taxonomies")
        return df

    def _init_layout(self) -> None:
        try:
            # App explanation
            correlations_layout = dmc.Stack([
                dmc.Title("Correlations Analysis", order=2),
                
                # Update explanation for correlations analysis
                dmc.Alert(
                    title="About Correlations Analysis",
                    children=[
                        "Analyze relationships between different taxa in your samples:",
                        html.Br(),
                        html.Ul([
                            html.Li("The heatmap shows correlation coefficients between pairs of taxa"),
                            html.Li([
                                "Correlation values range from -1 to 1:",
                                html.Ul([
                                    html.Li("Positive values (red): Taxa tend to increase together"),
                                    html.Li("Negative values (blue): As one taxon increases, the other decreases"),
                                    html.Li("Values near 0 (white): No clear relationship")
                                ])
                            ]),
                            html.Li("Use the controls to filter taxa and adjust the visualization"),
                            html.Li("Statistically significant correlations can be highlighted based on p-value")
                        ])
                    ],
                    color="blue",
                    variant="light",
                    style={"marginBottom": "20px"}
                ),
            ], gap="sm", p=0)
            
            # Control section
            control_section = dmc.Paper(
                children=[
                    dmc.Title("Analysis Controls", order=3),
                    dmc.Grid([
                        dmc.GridCol([
                            dmc.Select(
                                id="correlation-method",
                                label="Correlation Method",
                                description="Statistical method used to calculate correlations",
                                data=[
                                    {"value": "pearson", "label": "Pearson"},
                                    {"value": "spearman", "label": "Spearman"},
                                    {"value": "kendall", "label": "Kendall"},
                                ],
                                value=self.corr_method,
                                style={"width": "100%"},
                            ),
                        ], span=4),
                        dmc.GridCol([
                            dmc.NumberInput(
                                id="pvalue-cutoff",
                                label="P-value significance cutoff",
                                description="Display only correlations with p-value below this threshold",
                                value=0.05,
                                min=0.001,
                                max=1,
                                decimalSeparator=".",
                                decimalScale=3,
                                step=0.001,
                                stepHoldDelay=500,
                                stepHoldInterval=100,
                                style={"width": "100%"},
                            ),
                        ], span=4),
                        dmc.GridCol([
                            dmc.Switch(
                                id="filter-significant",
                                label="Show only significant correlations",
                                size="md",
                                radius="sm",
                                checked=False,
                                style={"marginTop": "25px"},
                            ),
                        ], span=4),
                    ]),
                    dmc.Space(h=20),
                    dmc.Button(
                        "Apply Settings",
                        id="apply-correlation-settings",
                        leftSection=[DashIconify(icon="mdi:refresh")],
                        color="blue",
                        variant="filled",
                        size="md",
                    ),
                ],
                p="md",
                withBorder=True,
                shadow="sm",
                radius="md",
                mb="md"
            )

            # File operations section
            file_section = dmc.Paper(
                children=[
                    dmc.Title("File Operations", order=3),
                    dmc.Group([
                        dmc.Button(
                            "Download Metadata CSV Template",
                            leftSection=[DashIconify(icon="mdi:file-download-outline")],
                            id="download-csv-button",
                        ),
                        dmc.Button(
                            "Download Correlations",
                            leftSection=[DashIconify(icon="mdi:chart-scatter-plot")],
                            id="btn-download-corr",
                        ),
                        dmc.Button(
                            "Download Plots as SVG",
                            leftSection=[DashIconify(icon="mdi:file-image-outline")],
                            id="btn-download-svg-correlations",
                        ),
                    ]),
                    dmc.Space(h=20),
                    dcc.Upload(
                        id='upload-data',
                        children=dmc.Button(
                            "Upload Metadata CSV",
                            leftSection=[DashIconify(icon="mdi:file-upload-outline")],
                            variant="outline",
                            fullWidth=True,
                        ),
                        multiple=False
                    ),
                ],
                p="md",
                withBorder=True,
                shadow="sm",
                radius="md",
                mb="md"
            )

            # Correlation heatmap
            heatmap_section = dmc.Paper(
                children=[
                    dmc.Title("Correlation Heatmap", order=3),
                    dcc.Graph(
                        id='heatmap-dendrogram',
                        figure=self.create_heatmap(),
                        style={"width": "100%", "height": "800px"}
                    )
                ],
                p="md",
                withBorder=True,
                shadow="sm",
                radius="md",
                mb="md",
                style={"width": "100%"}
            )

            # Combine all sections in a container with no padding
            container = dmc.Container(
                [
                    correlations_layout,
                    control_section,
                    file_section,
                    heatmap_section,
                    html.Div(id='notification-output'),
                    dcc.Download(id='download-metadata-csv'),
                    dcc.Download(id="download-corr-csv"),
                    dcc.Download(id="download-svg-diversity"),
                ],
                fluid=True,
                p=0,
                style={"maxWidth": "100%", "width": "100%", "margin": 0, "padding": 0}
            )

            self.app.layout = dmc.MantineProvider(
                theme={
                    "colorScheme": "light",
                    "primaryColor": "indigo",
                    "components": {
                        "Button": {"styles": {"root": {"fontWeight": 400}}},
                        "Title": {"styles": {"root": {"fontWeight": 500}}},
                        "Container": {"styles": {"root": {"maxWidth": "100%", "padding": 0, "margin": 0}}},
                        "Paper": {"styles": {"root": {"width": "100%"}}}
                    }
                },
                children=[container]
            )
            print("Correlations layout initialized")
        except Exception as e:
            error_msg = f"Error initializing Correlations layout: {str(e)}\n{traceback.format_exc()}"
            print(error_msg)
            self.app.layout = html.Div([
                html.H1("Error loading Correlations app"),
                html.Pre(error_msg)
            ])

    def compute_correlations_for_taxonomies(self, taxonomy_df, metadata_df, method, taxonomic_level='taxonomy'):
        taxonomy_df = taxonomy_df.drop_duplicates(subset=['sample_id', 'taxonomy'])

        if metadata_df.empty or taxonomy_df.empty:
            print("No metadata or taxonomy data found")
            return pd.DataFrame()

        # Merging the DataFrames
        merged_df = pd.merge(taxonomy_df, metadata_df, on=['sample_id', 'user_id'], how='right')

        # Extract metadata keys from 'data' column
        first_valid_dict = next((item for item in merged_df['data'] if isinstance(item, dict)), None)
        if first_valid_dict is None:
            print("No valid metadata found in 'data' column")
            return pd.DataFrame()
        metadata_keys = list(first_valid_dict.keys())

        # Extract metadata keys into separate columns
        for key in metadata_keys:
            merged_df[key] = merged_df['data'].apply(lambda x: x.get(key) if isinstance(x, dict) else None)

        # Initialize a list to collect correlation data
        correlation_data = []
        pvalue_data = []

        # Calculate correlation for each taxonomy
        unique_taxonomies = merged_df[taxonomic_level].unique()
        for taxonomy in unique_taxonomies:
            # Select abundance for the current taxonomy, fill NaN with 0
            taxonomy_abundance = merged_df[merged_df[taxonomic_level] == taxonomy]['abundance'].fillna(0)
            corr_row = {'taxonomy': taxonomy}
            pval_row = {'taxonomy': taxonomy}
            
            for meta_key in metadata_keys:
                if pd.api.types.is_numeric_dtype(merged_df[meta_key]):
                    # Check if there is variance in the data
                    if taxonomy_abundance.std() != 0 and merged_df[meta_key].std() != 0:
                        try:
                            # Get valid pairs (drop NaN values)
                            valid_data = pd.DataFrame({
                                'abundance': taxonomy_abundance,
                                'metadata': merged_df[meta_key]
                            }).dropna()
                            
                            # Calculate correlation coefficient
                            if method == 'pearson':
                                corr, p_value = stats.pearsonr(valid_data['abundance'], valid_data['metadata'])
                            elif method == 'spearman':
                                corr, p_value = stats.spearmanr(valid_data['abundance'], valid_data['metadata'])
                            elif method == 'kendall':
                                corr, p_value = stats.kendalltau(valid_data['abundance'], valid_data['metadata'])
                            else:
                                # Default to pandas corr method as fallback
                                corr = taxonomy_abundance.corr(merged_df[meta_key], method=method)
                                p_value = np.nan
                                
                            corr_row[meta_key] = corr
                            pval_row[meta_key] = p_value
                            
                        except Exception as e:
                            print(f"Error calculating correlation for {taxonomy} and {meta_key}: {str(e)}")
                            corr_row[meta_key] = np.nan
                            pval_row[meta_key] = np.nan
                    else:
                        corr_row[meta_key] = 0  # Assign 0 for zero variance
                        pval_row[meta_key] = 1.0  # p-value of 1.0 for zero variance (no correlation)
                else:
                    corr_row[meta_key] = np.nan
                    pval_row[meta_key] = np.nan
                    
            correlation_data.append(corr_row)
            pvalue_data.append(pval_row)

        # Create DataFrames from the collected data
        correlation_matrix = pd.DataFrame(correlation_data)
        pvalue_matrix = pd.DataFrame(pvalue_data)
        
        self.correlation_matrix = correlation_matrix
        self.pvalue_matrix = pvalue_matrix

        return correlation_matrix

    def get_all_meta(self):
        meta = Metadata.objects.filter(user_id=self.user_id)
        if not meta.exists():
            print(f"Found no metadata for user id {self.user_id}")
            return pd.DataFrame()  # Return an empty DataFrame if no records are found
        meta_df = pd.DataFrame.from_records(meta.values())
        return meta_df

    def create_heatmap(self, p_value_cutoff=0.05, filter_significant=False):
        if self.correlation_matrix is None or self.correlation_matrix.empty:
            return go.Figure()  # Return an empty figure if no data is available
            
        # Print debug info
        print(f"Creating heatmap with p_value_cutoff={p_value_cutoff}, filter_significant={filter_significant}")
            
        # Apply p-value filtering if requested
        filtered_corr_matrix = self.correlation_matrix.copy()
        
        if hasattr(self, 'pvalue_matrix') and filter_significant:
            print(f"Applying p-value filtering with cutoff {p_value_cutoff}")
            # Create mask where p-values are > cutoff (not significant)
            mask = self.pvalue_matrix.copy()
            
            # Convert dataframes to have the same structure
            mask = mask.set_index('taxonomy')
            filtered_corr_matrix = filtered_corr_matrix.set_index('taxonomy')
            
            # Create mask of non-significant correlations
            for col in mask.columns:
                # Get count before filtering
                valid_before = filtered_corr_matrix[col].notna().sum()
                
                # Set correlations to NaN where p-value > cutoff
                mask_invalid = mask[col] > p_value_cutoff
                filtered_corr_matrix.loc[mask_invalid, col] = np.nan
                
                # Get count after filtering
                valid_after = filtered_corr_matrix[col].notna().sum()
                print(f"Column {col}: {valid_before} values before filtering, {valid_after} after filtering")
                
            # Reset index to keep taxonomy as column
            filtered_corr_matrix = filtered_corr_matrix.reset_index()
            
            print(f"After filtering: {filtered_corr_matrix.notna().sum().sum()} valid correlation values")
        
        # Get top 30 taxa by absolute correlation values, ignoring NaN values
        filtered_corr_matrix = filtered_corr_matrix.dropna(axis=0, how='all')
        if len(filtered_corr_matrix) > 0:
            abs_corr = filtered_corr_matrix.iloc[:, 1:].abs().mean(axis=1, skipna=True)
            # Get top taxa (or all if less than 30)
            n_taxa = min(30, len(abs_corr))
            if n_taxa > 0:
                top_indices = abs_corr.nlargest(n_taxa).index
                filtered_matrix = filtered_corr_matrix.iloc[top_indices]
            else:
                filtered_matrix = filtered_corr_matrix  # Use all if filtering resulted in empty data
        else:
            print("No data after filtering or empty correlation matrix")
            # Return an empty figure with a message
            empty_fig = go.Figure()
            empty_fig.add_annotation(
                text="No significant correlations found with the current settings",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=18)
            )
            return empty_fig
        
        matrix_for_heatmap = filtered_matrix.set_index('taxonomy')
        
        print(f"Final heatmap shape: {matrix_for_heatmap.shape}")
        
        # Reverse the y-axis labels for better visualization 
        y_labels = [f'<i>{elem}</i> ' for elem in matrix_for_heatmap.index][::-1]
        # Reverse the z data matrix vertically
        z_data = matrix_for_heatmap.values[::-1]
        
        # Add a trace with custom colorbar title
        heatmap_trace = go.Heatmap(
            z=z_data, 
            x=matrix_for_heatmap.columns, 
            y=y_labels, 
            colorscale='Portland',
            colorbar=dict(
                title=dict(
                    text=f'{self.corr_method.capitalize()} Correlation',
                    side='right'
                )
            ),
            hovertemplate='<b>%{y}</b> vs <b>%{x}</b><br>Correlation: %{z:.3f}<extra></extra>'
        )

        fig = go.Figure(data=heatmap_trace)
        
        # Add a more informative title based on filtering
        title = 'Heatmap of Taxonomy-Metadata Correlations'
        if filter_significant:
            title += f' (p < {p_value_cutoff})'
        
        fig.update_layout(
            xaxis=dict(type='category', automargin=True, tickfont=dict(size=14)),
            yaxis=dict(type='category', automargin=True, tickfont=dict(size=14), autorange="reversed"),
            height=800,
            title=title,
            xaxis_title='Metadata',
            yaxis_title='Taxonomy',
            plot_bgcolor='rgba(0,0,0,0)',
            font=dict(size=14),
            margin=dict(l=10, r=10, t=80, b=10),
            autosize=True
        )
        return fig

    def handle_template_download(self):
        """Helper function for downloading template, to be called from index.py"""
        df = pd.DataFrame(columns=['sample_id', 'meta_1', 'meta_2', 'meta_3'])
        return dcc.send_data_frame(df.to_csv, "metadata_template.csv", index=False)

    def handle_correlations_download(self):
        """Helper function for downloading correlations, to be called from index.py"""
        if self.correlation_matrix is None:
            return None
            
        if hasattr(self, 'pvalue_matrix'):
            # Create a combined dataframe that includes both correlations and p-values
            corr_df = self.correlation_matrix.copy().set_index('taxonomy')
            pval_df = self.pvalue_matrix.copy().set_index('taxonomy')
            
            # Rename columns in p-value dataframe to distinguish them
            pval_df.columns = [f"{col}_pvalue" for col in pval_df.columns]
            
            # Combine correlations and p-values
            combined_df = pd.concat([corr_df, pval_df], axis=1)
            combined_df = combined_df.reset_index()
            
            return dcc.send_data_frame(combined_df.to_csv, f"{self.corr_method}_correlations_with_pvalues.csv", index=False)
        else:
            return dcc.send_data_frame(self.correlation_matrix.to_csv, f"{self.corr_method}_correlations.csv", index=False)

    def handle_file_upload(self, contents, filename):
        """Helper function for handling file upload, to be called from index.py"""
        if contents is None:
            return dmc.Alert("No file uploaded", color="red")

        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        try:
            if 'csv' in filename:
                df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
            elif 'xls' in filename:
                df = pd.read_excel(io.BytesIO(decoded))
            else:
                return dmc.Alert("Unsupported file type", color="red")

            # Process the uploaded file
            self.meta_df = df
            
            # Recalculate correlations with current method
            self.correlation_scores = self.compute_correlations_for_taxonomies(
                self.nanopore_df, self.meta_df, self.corr_method, 'taxonomy')

            return dmc.Alert("File uploaded and processed successfully. Click 'Apply Settings' to update the visualization.", color="green")
        except Exception as e:
            return dmc.Alert(f"Error processing file: {str(e)}", color="red")

    def handle_svg_download(self):
        """Helper function for downloading SVG, to be called from index.py"""
        # Placeholder for SVG download
        return None


def _force_list(x: Union[Any, List]) -> List[Any]:
    """
    Force object into a list if object is not a list.
    """

    return x if isinstance(x, list) else [x]
