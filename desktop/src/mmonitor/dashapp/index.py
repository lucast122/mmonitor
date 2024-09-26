from dash import dcc
from dash import html
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
from flask import request

from mmonitor.dashapp.app import app
from mmonitor.dashapp.apps import correlations, taxonomy, horizon, kegg
from mmonitor.dashapp.base_app import BaseApp
from mmonitor.database.mmonitor_db import MMonitorDBInterface


class Index(BaseApp):
    """
    Landing page of the dash application.
    Contains navigation to the various application pages.
    """

    def __init__(self, sql: MMonitorDBInterface):
        super().__init__(sql)
        self._app = app
        self._init_apps()
        self._init_layout()
        app.layout = self.layout
        self._init_callbacks()

    def _init_apps(self) -> None:
        """
        Register apps by their urls, names and instances.
        This is the only place you need to add an app.
        """

        self._apps = {
            '/apps/taxonomy': {
                'name': 'Taxonomy',
                'app': taxonomy.Taxonomy(self._sql)

            },

            '/apps/horizon': {
                'name': 'Horizon',
                'app': horizon.Horizon(self._sql)
            },

            '/apps/correlations': {
                'name': 'Correlations',
                'app': correlations.Correlations(self._sql)
                # },
                # '/apps/kraken2': {
                #     'name': 'Kraken2',
                #     'app': kraken.Kraken(self._sql)
            },
            '/apps/kegg': {
                'name': 'KEGG',
                'app': kegg.Kegg()
            },
            # '/apps/genome_browser': {
            #     'name': 'Gene Browser',
            #     'app': (genome_browser.GenomeBrowser(self._sql))
            #     }

        }

    def _init_layout(self) -> None:
        """
        The index page's layout consists of the app navigation and
        the currently selected app's page content.
        """

        location = dcc.Location(id='url', refresh=False)
        navigation = html.Div([
            dcc.Link(values['name'], href=url, style={'padding': '10px', 'font-size': "30px", "font-weight" : "bold",  "hover" : "#B22222:"})
            for url, values in self._apps.items()
        ], className="row")
        page_content = html.Div(id='page-content', children=[])
        container = html.Div([location, navigation, page_content])
        self.layout = container

    def _init_callbacks(self) -> None:
        @app.callback(
            Output('page-content', 'children'),
            Input('url', 'pathname')
        )
        def display_page(pathname) -> html:
            """
            Change page content to the selected app
            if the url is valid.
            """

            # empty urls have no effect
            if pathname is None or pathname == '/':
                raise PreventUpdate
            # attempt to change to selected app
            elif pathname in self._apps:
                return self._apps[pathname]['app'].layout
            # otherwise it's page not found
            else:
                return "404 Page Error! Please choose a link"

        @app.server.route('/shutdown', methods=['POST'])
        def shutdown():
            """
            Terminate the dash app.
            """
            func = request.environ.get('werkzeug.server.shutdown')
            if func is None:
                raise RuntimeError('Not running with the Werkzeug Server')
            func()
            return 'Server is shutting down...'

    def run_server(self, debug: bool) -> None:
        """
        Runs the dash app.
        """
        self._app.run_server(debug=debug)
