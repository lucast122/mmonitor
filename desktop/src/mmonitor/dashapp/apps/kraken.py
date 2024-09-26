import base64
from io import StringIO
from typing import List, Tuple, Set

import dash
import pandas as pd
import plotly.express as px
from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from plotly.graph_objects import Figure

from mmonitor.dashapp.app import app
from mmonitor.dashapp.base_app import BaseApp
from mmonitor.database.mmonitor_db import MMonitorDBInterface


class Kraken(BaseApp):
    """
    App to upload and display taxonomies computed using Kraken2 in fastq.gz_report format.
    """

    def __init__(self, sql: MMonitorDBInterface):
        super().__init__(sql)
        self._init_layout()
        self._init_callbacks()
        self._raw_dfs: List[pd.DataFrame] = []
        self._ranks: Set[str] = set()

    def _init_layout(self) -> None:
        """
        Kraken layout consists of an upload field and the resulting graph.
        """

        # data upload
        header = html.H1(children='Bioreactor taxonomy computed using Kraken2')
        upload_header = html.H4(children='Please input kraken2 reports')
        upload = dcc.Upload(
            id='upload-data',
            children='Drag and Drop or Select Files',
            style={
                'width': '100%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px'
            },
            multiple=True
        )

        # pipeline_out content
        graph = dcc.Graph(id='kraken-graph')

        # slider to adjust entities per bar
        # capped max value for performance reasons
        slider_label = html.Div("Entities per Bar:", style={'margin': '5px'})
        slider = dcc.Slider(id='slider', min=1, max=50, step=1, value=10, tooltip={"placement": "bottom"})
        slider_container = html.Div([slider], style={'width': '25%', 'margin': '5px'})
        slider_content = html.Div(
            [slider_label, slider_container],
            style={'display': 'flex'}
        )

        # dropdown menu to dynamically select taxonomic rank
        tax_rank_label = html.Div("Taxonomic rank:", style={'margin': '5px'})
        tax_rank = dcc.Dropdown(
            id='tax-rank-dd',
            clearable=False,
        )
        tax_rank_content = html.Div(
            [tax_rank_label, tax_rank],
            style={'display': 'flex'}
        )

        # hide content before first upload
        content = html.Div([graph, slider_content, tax_rank_content], id='content', style={'display': 'none'})
        container = html.Div([header, upload_header, upload, content])
        self.layout = container

    def _init_callbacks(self) -> None:
        @app.callback(
            Output('kraken-graph', 'figure'),
            Output('tax-rank-dd', 'options'),
            Output('tax-rank-dd', 'value'),
            Output('content', 'style'),
            Input('upload-data', 'contents'),
            Input('slider', 'value'),
            Input('tax-rank-dd', 'value'),
            State('upload-data', 'filename'),
            State('tax-rank-dd', 'options'),
            State('content', 'style'),
        )
        def _update_output(contents, slider_value, tax_rank, filenames, tax_rank_options, content_style):
            """
            Read the uploaded files and display the pipeline_out plot.
            The plot is a stacked bar chart with abundances
            projected to 100%. The number of entities per bar
            is adjustable with different elements.
            """

            # check if the callback was triggered by the slider/dropdown
            # that means no file was uploaded -> don't recompute raw data
            ctx = dash.callback_context
            if ctx.triggered:
                element_id = ctx.triggered[0]['prop_id'].split('.')[0]
                if element_id == 'slider' or element_id == 'tax-rank-dd':
                    kraken_figure = self._generate_kraken_figure(slider_value, tax_rank)
                    return kraken_figure, tax_rank_options, tax_rank, content_style

            # otherwise files were uploaded -> (re)compute raw data
            # catch cancelled uploads and wrong formats
            if contents is None:
                raise PreventUpdate
            if any(not f.endswith('fastq.gz_report') for f in filenames):
                raise PreventUpdate

            # decode and parse uploaded content
            contents = [_decode_base64content(c) for c in contents]
            content_samples = [(c, f.split('.')[0]) for c, f in zip(contents, filenames)]
            self._parse_data(content_samples)

            # extract present tax ranks
            tax_rank_options = [{'label': t, 'value': t} for t in self._ranks]
            tax_rank_options.sort(key=lambda t: t['label'])

            # select a tax rank, preferably species
            tax_rank = 'S' if 'S' in self._ranks else tax_rank_options[0]['value']

            # generate graph
            kraken_figure = self._generate_kraken_figure(slider_value, tax_rank)

            content_style['display'] = 'block'
            return kraken_figure, tax_rank_options, tax_rank, content_style

    def _parse_data(self, content_sample_names: List[Tuple[str, str]]):
        """
        Parse uploaded data into a raw data format that can later be
        used to compute data to generate the graph.
        """

        self._raw_dfs = []
        self._ranks = set()

        for content, sample_name in content_sample_names:
            # read the contents and add a header to the table
            df = pd.read_csv(
                StringIO(content),
                sep='\t',
                header=None,
                usecols=[1, 3, 5],
                names=['Count', 'Rank', 'Name']
            )

            # sort by value
            df = df.sort_values('Count', ascending=False)

            # format name
            df['Name'] = df['Name'].apply(lambda s: s.strip())

            # add sample name
            df['Sample'] = sample_name

            # append to raw data and available tax ranks
            self._raw_dfs.append(df)
            self._ranks |= set(df['Rank'].unique())

    def _generate_kraken_figure(self, entities_per_bar: int = 10, tax_rank: str = 'S') -> Figure:
        """
        Generate a plot from content/sample name tuples. The plot is a stacked bar chart
        with the y-axis projected to 100%. That means a taxonomy is represented by how
        many percent of the total abundances the abundance of the taxonomy accounts for.
        """

        dfs = []
        for raw_df in self._raw_dfs:
            df = raw_df.copy()

            # filter for the specified taxonomic rank,
            # and pop specified amount of entries
            df = df[df['Rank'] == tax_rank]
            df = df.drop(columns='Rank')
            df = df.head(entities_per_bar)

            # replace counts by percentage distribution in relation to total counts
            count_sum = sum(df['Count'])
            df['Count'] = df['Count'].apply(lambda c: c / count_sum * 100)
            df['Percent'] = df['Count'].rename('Percent')

            # append and finally concatenate all dataframes
            dfs.append(df)
        df = pd.concat(dfs)

        fig = px.bar(df, x='Sample', y='Percent', color='Name')
        return fig


def _decode_base64content(content: str) -> str:
    """
    Decode base64 content according to wise computer sorcerers
    https://stackabuse.com/encoding-and-decoding-base64-strings-in-python/
    """

    base64_message = content.split(',')[1]
    base64_bytes = base64_message.encode('utf-8')
    message_bytes = base64.b64decode(base64_bytes)
    message = message_bytes.decode('utf-8')
    return message
