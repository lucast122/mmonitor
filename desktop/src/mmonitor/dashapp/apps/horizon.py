import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
from plotly.subplots import make_subplots

from mmonitor.dashapp.app import app
from mmonitor.dashapp.base_app import BaseApp
from mmonitor.database.mmonitor_db import MMonitorDBInterface


class Horizon(BaseApp):
    """
    App to display the taxonomy abundances in a horizon plot
    """
    def __init__(self, sql: MMonitorDBInterface):
        super().__init__(sql)
        # self.app = dash.Dash(__name__)
        q = 'SELECT sample_id, taxonomy, abundance FROM mmonitor'
        self.df = self._sql.query_to_dataframe(q)
        self._init_layout()
        self._init_callbacks()

    def _init_layout(self):
        div = html.Div([
            dcc.Graph(id='horizon-plot'),
            dcc.Interval(
                id='interval-component',
                interval=1 * 1000,  # in milliseconds
                n_intervals=0
            )
        ])

        container = dbc.Container(
            [
                div
            ],
            fluid=True
        )

        self.layout = container

    def _init_callbacks(self):
        @app.callback(
            Output('horizon-plot', 'figure'),
            [Input('interval-component', 'n_intervals')]
        )
        def update_graph_live(n):
            q = 'SELECT sample_id, taxonomy, abundance FROM mmonitor'
            self.df = self._sql.query_to_dataframe(q)

            fig = make_subplots(rows=1, cols=1)

            # Assuming the 'abundance' column contains the data to be plotted

            fig.add_trace(go.Heatmap(x=self.df['sample_id'], y=self.df['taxonomy'], z=self.df['abundance']), row=1,
                          col=1)
            print(self.df['abundance'])

            fig.update_layout(height=1200, width=600, title_text="Horizon Plot")

            return fig

        # self.app.run_server(debug=True)
