from typing import Tuple, Any, List, Iterable, Dict, Union

import dash
import plotly.express as px
from dash import dash_table
from dash import dcc as dcc
from dash import html as html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from pandas import DataFrame
from plotly.graph_objects import Figure

from mmonitor.calculations.stats import scipy_correlation
from mmonitor.dashapp.app import app
from mmonitor.dashapp.base_app import BaseApp
from mmonitor.database.mmonitor_db import MMonitorDBInterface


class Correlations(BaseApp):
    """
    App to display correlations between abundance and metadata and their probabilities.
    The graph and data table are independent.
    The data table content may be downloaded in csv format.

    WARNING: - correlations can only be calculated for complete and sanitized data
             - bootstrapping is slow
    """

    def __init__(self, sql: MMonitorDBInterface):
        super().__init__(sql)
        self._taxonomies = self._sql.get_unique_taxonomies()
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
        self._init_layout()
        self._init_callbacks()

    def _init_layout(self) -> None:
        """
        Correlations layout consists of
            - Graphs containing correlation and test scores
            - Data table containing correlations and test scores
        """

        # graph dropdown menus in a flex box
        header = html.H1('Correlations of abundances and metadata')

        taxonomy_dd = dcc.Dropdown(
            id='taxonomy-dd',
            options=[{'label': t, 'value': t} for t in self._taxonomies],
            value=self._taxonomies[0],
            style={'flex-grow': '1'},
            clearable=False,
        )
        methods_dd = dcc.Dropdown(
            id='methods-dd',
            options=[{'label': m, 'value': m} for m in self._methods],
            value=self._methods[0],
            style={'flex-grow': '1'},
            clearable=False,
        )
        tests_dd = dcc.Dropdown(
            id='tests-dd',
            options=[{'label': t, 'value': t} for t in self._tests],
            value=self._tests[0],
            style={'flex-grow': '1'},
            clearable=False,
        )
        dropdowns = html.Div(
            [taxonomy_dd, methods_dd, tests_dd],
            style={'display': 'flex', 'justify-content': 'space-between'}
        )

        # graphs
        graph_score = dcc.Graph(id='graph-score', figure={'data': []})
        graph_test = dcc.Graph(id='graph-test', figure={'data': []})

        # data table dropdown menus in a flexbox
        header_tb = html.H2('Select correlation data and export')
        taxonomy_dd_tb = dcc.Dropdown(
            id='taxonomy-dd-tb',
            options=[{'label': t, 'value': t} for t in self._taxonomies],
            value=self._taxonomies,
            style={'flex-grow': '1'},
            clearable=False,
            multi=True,
        )
        methods_dd_tb = dcc.Dropdown(
            id='methods-dd-tb',
            options=[{'label': m, 'value': m} for m in self._methods],
            value=self._methods[0],
            style={'flex-grow': '1'},
            clearable=False,
            multi=True,
        )
        tests_dd_tb = dcc.Dropdown(
            id='tests-dd-tb',
            options=[{'label': t, 'value': t} for t in self._tests],
            value=[],
            style={'flex-grow': '1'},
            clearable=False,
            multi=True,
        )
        dropdowns_tb = html.Div(
            [taxonomy_dd_tb, methods_dd_tb, tests_dd_tb],
            style={'display': 'flex', 'justify-content': 'space-between'}
        )

        # table interaction buttons in flexbox
        apply_btn_tb = html.Button(
            'Apply',
            id='apply-btn-tb',
            style={'margin': '5px'}
        )
        select_all_btn_tb = html.Button(
            'Select All',
            id='select-all-btn-tb',
            style={'margin': '5px'}
        )
        clear_selection_btn_tb = html.Button(
            'Clear Selection',
            id='clear-selection-btn-tb',
            style={'margin': '5px'}
        )
        export_btn_tb = html.Button(
            'Export',
            id='export-btn-tb',
            style={'margin': '5px'}
        )
        selections_tb = html.Div(
            [apply_btn_tb, select_all_btn_tb, clear_selection_btn_tb, export_btn_tb],
            style={'display': 'flex', 'padding': '10px'}
        )

        # data table and its download element
        correlations_tb = dash_table.DataTable(id='table-correlations', data=[], columns=[])
        download_tb = dcc.Download(id='download-tb')

        CONTENT_STYLE = {

            "margin-right": "2rem",
            "padding": "2rem 1rem",
            'margin-bottom': '200px',
            'font-size': '21px'
        }

        container = html.Div(
            [header, dropdowns, graph_score, graph_test, header_tb, dropdowns_tb,
             selections_tb, correlations_tb, download_tb],
            style=CONTENT_STYLE)

        self.layout = container

    def _init_callbacks(self) -> None:
        @app.callback(
            Output('graph-score', 'figure'),
            Output('graph-test', 'figure'),
            Input('taxonomy-dd', 'value'),
            Input('methods-dd', 'value'),
            Input('tests-dd', 'value')
        )
        def _update_graphs(taxonomy: str, method: str, test: str) -> Tuple[Figure, Figure]:
            """
            Populate the correlations and probability graphs
            according to dropdown selections.
            """

            # request data from database and filter columns containing metadata
            # abundance and metadata tables are joined to align entries by their sample_ids
            df = self._sql.get_abundance_meta_by_taxonomy(taxonomy)
            meta_cols = [col for col in df.columns if col not in ['sample_id', 'abundance']]

            # for each metadata calculate the correlation and probability
            # transpose the list of tuples and extract each a correlation and a probability list:
            # [(c1, p1), (c2, p2), (c3, p3), ...] -> [c1, c2, c3, ...], [p1, p2, p3, ...]
            y_axis_tuple = [scipy_correlation(df['abundance'], df[meta], test, method) for meta in meta_cols]
            y_axis_score, y_axis_test = zip(*y_axis_tuple)
            y_axis_test = [xs[test] for xs in y_axis_test]

            # generate plots
            fig_score = px.scatter(x=meta_cols, y=y_axis_score, labels={'x': 'Metric', 'y': f'{method} Score'})
            fig_test = px.scatter(x=meta_cols, y=y_axis_test, labels={'x': 'Metric', 'y': f'{test} Score'})

            return fig_score, fig_test

        @app.callback(
            Output('table-correlations', 'columns'),
            Output('table-correlations', 'data'),
            Output('taxonomy-dd-tb', 'value'),
            Input('apply-btn-tb', 'n_clicks'),
            Input('select-all-btn-tb', 'n_clicks'),
            Input('clear-selection-btn-tb', 'n_clicks'),
            State('taxonomy-dd-tb', 'value'),
            State('methods-dd-tb', 'value'),
            State('tests-dd-tb', 'value'),
            State('table-correlations', 'columns'),
            State('table-correlations', 'data'),
        )
        def _update_table(x: int, y: int, z: int,
                          taxonomies: Union[str, List[str]],
                          methods: Union[str, List[str]],
                          tests: Union[str, List[str]],
                          tb_columns: Any,
                          tb_data: Any,
                          ) -> Tuple[Iterable, Dict, List]:
            """
            Populate the data table with correlation scores and their probabilities
            according to the dropdown selections. This table contains data more suited
            for higher dimensions, but in order to keep it displayable and exportable
            it is complexly condensed into two dimensions.

            The final dataframe/table will have this form:

            Taxonomy  | Method or Test    | Meta1  | Meta2  | Meta3  | ...
            ---------------------------------------------------------------
            taxonomy1 | Pearson           | score1 | score2 | score3 | ...
            taxonomy1 | Pearson T-Test    | score1 | score2 | score3 | ...
            taxonomy1 | Pearson Bootstrap | score1 | score2 | score3 | ...
            taxonomy1 | Kendall           | score1 | score2 | score3 | ...
            taxonomy1 | Kendall T-Test    | score1 | score2 | score3 | ...
            taxonomy1 | Kendall Bootstrap | score1 | score2 | score3 | ...
            taxonomy2 | Pearson           | score1 | score2 | score3 | ...
            taxonomy2 | Pearson T-Test    | score1 | score2 | score3 | ...
            taxonomy2 | Pearson Bootstrap | score1 | score2 | score3 | ...
            taxonomy2 | Kendall           | score1 | score2 | score3 | ...
            taxonomy2 | Kendall T-Test    | score1 | score2 | score3 | ...
            taxonomy2 | Kendall Bootstrap | score1 | score2 | score3 | ...
            ...
            """

            # figure out which button was clicked
            # default behaviour is 'apply'
            ctx = dash.callback_context
            button_id = 'apply-btn-tb'
            if ctx.triggered:
                button_id = ctx.triggered[0]['prop_id'].split('.')[0]

            # keep table, adjust selection in taxonomies dropdown menu
            if button_id == 'select-all-btn-tb':
                return tb_columns, tb_data, self._taxonomies
            elif button_id == 'clear-selection-btn-tb':
                return tb_columns, tb_data, []

            # format dropdown selections and eject if selections are invalid
            taxonomies: List[str] = _force_list(taxonomies)
            methods: List[str] = _force_list(methods)
            tests: List[str] = _force_list(tests)
            if len(taxonomies) == 0 or len(methods) == 0:
                raise PreventUpdate

            # request data from database and filter for metadata columns
            meta_df = self._sql.get_all_meta()
            meta_cols = [meta for meta in meta_df.columns if meta != 'sample_id']

            # prepare the resulting dataframe/table
            columns = ['Taxonomy', 'Method or Test', *meta_cols]
            df = DataFrame(columns=columns, index=None)

            # for each selected taxonomy
            for taxonomy in taxonomies:

                # request corresponding abundances from database
                tax_df = self._sql.get_abundance_by_taxonomy(taxonomy)
                joined_df = tax_df.join(meta_df, on='sample_id', how='outer', lsuffix='_t', rsuffix='_m').dropna()

                # for each selected method
                for method in methods:

                    # prepare the row containing correlation scores per metadata
                    # and the rows containing the corresponding probabilities
                    # of the correlation score per metadata according to the selected tests
                    score_row = {'Taxonomy': taxonomy, 'Method or Test': method}
                    test_rows = {t: {'Taxonomy': taxonomy, 'Method or Test': ' '.join([method, t])} for t in tests}

                    # for each metadata available in the database
                    for meta in meta_cols:

                        # calculate the correlation score and its probability
                        # according to the selected tests
                        score, test_scores = scipy_correlation(joined_df['abundance'], joined_df[meta], tests, method)

                        # populate the row containing correlation scores per metadata
                        score_row[meta] = score

                        # then populate the rows containing the corresponding probabilities
                        # of the correlation score per metadata according to the selected tests
                        for test, test_score in test_scores.items():
                            test_rows[test][meta] = test_score

                    # finally append the correlation row and the probability rows the the dataframe
                    df = df.append(score_row, ignore_index=True)
                    for test in tests:
                        df = df.append(test_rows[test], ignore_index=True)

            # return data for the dash data table according to
            # https://dash.plotly.com/datatable
            self._table_df = df
            return [{"name": i, "id": i} for i in df.columns], df.to_dict('records'), taxonomies

        @app.callback(
            Output('download-tb', 'data'),
            Input('export-btn-tb', 'n_clicks'),
        )
        def _export_table(n_clicks) -> Dict[Any, str]:
            """
            Download data table content.
            """

            if n_clicks is None:
                raise PreventUpdate
            return dict(content=self._table_df.to_csv(index=False), filename='correlations.csv')


def _force_list(x: Union[Any, List]) -> List[Any]:
    """
    Force object into a list if object is not a list.
    """

    return x if isinstance(x, list) else [x]
