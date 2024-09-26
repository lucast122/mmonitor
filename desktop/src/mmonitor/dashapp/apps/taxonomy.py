import base64
from io import StringIO
from typing import Tuple, List, Any, Dict

import dash_bootstrap_components as dbc
import pandas as pd
import plotly.express as px
from dash import dash_table
from dash import dcc
from dash import html, dash
from dash.dependencies import Input, Output
from plotly.graph_objects import Figure

from mmonitor.dashapp.app import app
from mmonitor.dashapp.base_app import BaseApp
from mmonitor.database.mmonitor_db import MMonitorDBInterface


class Taxonomy(BaseApp):
    """
    App to display the abundances of taxonomies in various formats.
    """

    def __init__(self, sql: MMonitorDBInterface):
        """
          Initialize the Taxonomy app.

          Args:
              sql (MMonitorDBInterface): An interface to the SQL database.
        """
        super().__init__(sql)
        # initialize bootstrap
        app = dash.Dash(__name__, external_stylesheets=[dbc.themes.SANDSTONE])


        q = f"SELECT * FROM mmonitor"
        self.df = self._sql.query_to_dataframe(q)
        self.df_sorted = self.df.sort_values(by=["sample_id", "abundance"], ascending=[True, False])
        self.unique_sample_ids = self._sql.query_to_dataframe("SELECT DISTINCT sample_id FROM mmonitor")[
            "sample_id"].tolist()

        # Get the number of unique values in each column
        # get number of unique taxonomies for creation of slider. limit max taxa to plot to 100
        self.unique_counts = min(self.df.nunique()[1], 100)

        self._init_layout()
        self._init_callbacks()


    def _init_layout(self) -> None:
        """
        Taxonomy app layout consists of two graphs visualizing the abundances
        of taxonomies and an optional pie chart as well as a data table for
        debugging purposes.
        """

        header = dbc.Row(
            dbc.Col(html.H1(children='Taxonomic composition'), width={'size': 12, 'offset': 0}),
            justify="center",
        )
        sample_select_dropdown = dbc.Row([
            dbc.Col(
                dbc.Label("Samples to plot:"),
                width={"size": 1, "offset": 0}
            ),
            dbc.Col(
                dcc.Dropdown(
                    id='sample_select_value',
                    options=[{'label': i, 'value': i} for i in self.unique_sample_ids],
                    multi=True  # allow multiple selections
                ),
                width={"size": 8, "offset": 0}
            ),
        ])

        # dropdown menu to select chart type
        # initialize as stacked bar chart
        demo_dd = dbc.Row([
            dbc.Col(
                dbc.Label("Plot type:"),
                width={"size": 1, "offset": 0}
            ),
            dbc.Col(
                dcc.Dropdown(
                    id='dropdown',
                    options=[
                        {'label': 'Stacked Barchart', 'value': 'stackedbar'},
                        {'label': 'Grouped Barchart', 'value': 'groupedbar'},
                        {'label': 'Area plot', 'value': 'area'},
                        {'label': 'Pie chart', 'value': 'pie'},
                        {'label': 'Scatter plot', 'value': 'scatter'},
                        {'label': 'Scatter 3D', 'value': 'scatter3d'}],
                    value='stackedbar'
                ),
                width={"size": 3, "offset": 0}
            ),
        ], justify="start")

        # graph elements
        graph1 = dcc.Graph(id='graph1', figure={'data': []})
        graph2 = dcc.Graph(id='graph2', figure={'data': []})

        header_pie_chart_sample_select_dbc = dbc.Row(
            dbc.Col(html.H1(children='Select a sample to display.'), width={'size': 10, 'offset': 0}),
            justify="center", style={'display': 'none'}, id='header_pie_chart_sample_select_dbc')

        slider_header = dbc.Row(
            dbc.Col(html.H1(children='Select the number of taxa to display:'), width={'size': 12, 'offset': 0}),
            justify="start")

        slider = dbc.Row(dbc.Col(dcc.Slider(
            id='slider',
            min=1,
            max=self.unique_counts,
            value=10,
            marks={i: str(i) if i in [1, self.unique_counts // 2, self.unique_counts] else "" for i in
                   range(1, self.unique_counts + 1)},

            step=1
        ), width={'size': 10, 'offset': 0}), justify='start')

        pie_chart_input = dcc.Dropdown(
            id='number_input_piechart',
            options=[{'label': t, 'value': t} for t in self._sql.get_unique_samples()],
            value=self._sql.get_unique_samples()[0],
            style={'display': 'none'},
            clearable=False,
        )

        # data table for debugging

        db_header = dbc.Row(dbc.Col(html.H4(children="Database"), width={'size': 12, 'offset': 0}),
                            justify="center")
        data, columns = self._generate_table_data_cols()
        data_tb = dbc.Row(dbc.Col(dash_table.DataTable(id='table-correlations', data=data, columns=columns),
                                  width={'size': 12, 'offset': 0}), justify="center")

        download_button = dbc.Button("Download CSV", id="btn-download")
        download_component = dcc.Download(id="download-csv")

        container = dbc.Container(
            [
                header,
                sample_select_dropdown,
                demo_dd,
                graph1,
                graph2,
                slider_header,
                slider,
                header_pie_chart_sample_select_dbc,
                pie_chart_input,
                db_header,
                download_button,
                download_component,
                data_tb
            ],
            fluid=True
        )

        self.layout = container

    def _generate_table_data_cols(self, max_rows=40) -> Tuple[List[Any], List[Any]]:
        """
        Generate data to populate a dash data table with.
        A dash data table requires the data in a dict format as well as a collection of mapped column names.

        Args:
        max_rows (int, optional): The maximum number of rows to include in the table. Defaults to 40.

        Returns:
        Tuple[List[Any], List[Any]]: A tuple containing the data for the table in dict format and a collection of mapped column names.
        """

        q = f"SELECT * FROM mmonitor LIMIT {max_rows}"
        df = self._sql.query_to_dataframe(q)
        return df.to_dict('records'), [{'name': i, 'id': i} for i in df.columns]

    def plot_stacked_bar(self, df):
        fig1 = px.bar(df, x="sample_id", y="abundance", color="taxonomy", barmode="stack")
        fig1.update_layout(
            legend=dict(
                orientation="v",
                y=1,
                x=1.1
            ),
            margin=dict(
                l=100,  # Add left margin to accommodate the legend
                r=100,  # Add right margin to accommodate the legend
                b=100,  # Add bottom margin
                t=100  # Add top margin
            ),
            width=2000,
            height=1000
        )

        fig2 = px.bar(df, x="taxonomy", y="abundance", color="sample_id", barmode="stack")
        return fig1, fig2

    def plot_grouped_bar(self, df):
        # Plotting code for grouped bar goes here...
        fig1 = px.bar(df, x="sample_id", y="abundance", color="taxonomy", barmode="group")
        fig2 = px.bar(df, x="taxonomy", y="abundance", color="sample_id", barmode="group")

        return fig1, fig2

    def plot_scatter(self, df):
        fig1 = px.scatter(df, x="sample_id", y="abundance", size="abundance", color="sample_id")
        # fig2 = px.scatter(df, x="taxonomy", y="abundance", size="sample_id", color="sample_id")
        return fig1

    def plot_area(self, df):
        # Plotting code for area plot goes here...
        fig1 = px.area(df, x="sample_id", y="abundance", color="taxonomy", line_group="taxonomy")
        return fig1

    def plot_scatter_3d(self, df):
        # Plotting code for scatter 3D goes here...
        fig1 = px.scatter_3d(df, x='taxonomy', y='abundance', z='sample_id', color='taxonomy')
        fig2 = px.scatter_3d(df, x='abundance', y='taxonomy', z='sample_id', color='sample_id')
        return fig1, fig2

    def plot_pie(self, df,sample_value_piechart):
        # Plotting code for pie chart goes here...
        pie_values = df.loc[df["sample_id"] == sample_value_piechart, 'abundance']
        pie_names = df.loc[df["sample_id"] == sample_value_piechart, 'taxonomy']
        fig1 = px.pie(df, values=pie_values, names=pie_names,
                      title=f'Pie chart of bioreactor taxonomy of sample {sample_value_piechart}')
        piechart_style = {'display': 'block'}
        fig2_style = {'display': 'none'}
        return fig1, piechart_style, fig2_style

    def _init_callbacks(self) -> None:
        """
        Initialize the callbacks for the app.

        The callbacks include updating the figures when the dropdown value, piechart value, or slider value changes,
        and downloading the CSV when the download button is clicked.
        """

        # this updates the graph based on the samples selected

        @app.callback(
            Output('graph1', 'figure'),
            Output('graph2', 'figure'),
            # this output hides the pie chart number input when no pie chart is plotted
            Output('number_input_piechart', 'style'),
            # hides 2nd plot for pie chart
            Output('graph2', 'style'),
            Input('dropdown', 'value'),
            Input('number_input_piechart', 'value'),
            Input('slider', 'value'),
            Input('sample_select_value', 'value')

        )
        def plot_selected_figure(value, sample_value_piechart, slider_value, sample_select_value) -> Tuple[
            Figure, Figure, Dict[str, str]]:
            """
            Update the figures based on the selected value from the dropdown menu, the selected sample value for the piechart,
            and the selected number of taxa from the slider.

            Args:
            value (str): The selected value from the dropdown menu.
            sample_value_piechart (str): The selected sample value for the piechart.
            slider_value (int): The selected number of taxa from the slider.
            sample_select_value: (str) output of multi select dropdown (selected samples)

            Returns:
                Tuple[Figure, Figure, Dict[str, str]]: A tuple containing two figures for the graphs and a dict for the piechart style.
            """
            # fallback values
            fig1 = {'data': []}
            fig2 = {'data': []}
            piechart_style = {'display': 'none'}

            fig2_style = {'display': 'block'}
            # request necessary data from database
            # q = "SELECT sample_id, taxonomy, abundance FROM mmonitor"
            # df = self._sql.query_to_dataframe(q)
            df_sorted = self.df.sort_values("abundance", ascending=False)
            df_grouped = df_sorted.groupby("sample_id")
            result_df = pd.DataFrame()
            for name, group in df_grouped:
                # get most abundant for slider_value number of selected taxa for each sample to display in taxonomy plots
                sample_rows = group.head(slider_value)
                result_df = pd.concat([result_df, sample_rows])

            result_df = result_df.reset_index(drop=True)
            self.result_df = result_df
            self.df_selected = self.result_df[self.result_df['sample_id'].astype(str).isin(sample_select_value)]

            if value == 'stackedbar':
                fig1, fig2 = self.plot_stacked_bar(self.df_selected)

            elif value == 'groupedbar':
                fig1, fig2 = self.plot_grouped_bar(self.df_selected)

            elif value == 'scatter':
                fig1, fig2_style = self.plot_scatter(self.df_selected)

            elif value == 'area':
                fig1, fig2_style = self.plot_area(self.df_selected)

            elif value == 'scatter3d':
                fig1, fig2 = self.plot_scatter_3d(self.df_selected)

            elif value == "pie":
                fig1, piechart_style, fig2_style = self.plot_pie(self.result_df, sample_value_piechart)

            return fig1, fig2, piechart_style, fig2_style

        # @self.app.callback(
        # Output('graph1', 'figure'),
        # Input('sample_select_dropdown', 'value')
        # )
        # def update_graph(selected_samples):
        #     # Filter dataframe based on selected samples
        #     df_selected = self.df[self.df['sample_id'].isin(selected_samples)]
        #     print()

        #     # Update graph based on selected samples
        #     fig = self.plot_stacked_bar(df_selected)

        #     return fig

        # Add a new callback that updates the header's style based on the dropdown's value
        @app.callback(
            Output('header_pie_chart_sample_select_dbc', 'style'),
            Input('dropdown', 'value')
        )
        def show_hide_element(value):
            if value == "pie":
                return {'display': 'block'}  # Show the header if 'pie' is selected
            else:
                return {'display': 'none'}  # Hide the header for other options

        @app.callback(
            Output("download-csv", "data"),
            [Input("btn-download", "n_clicks")]
        )
        def download_csv(n_clicks):
            if n_clicks is not None:
                # Create a sample DataFrame for demonstration

                # Convert DataFrame to CSV string
                csv_string = self.df_sorted.to_csv(index=False)

                # Create a BytesIO object to hold the CSV data
                csv_bytes = StringIO()
                csv_bytes.write(csv_string)

                # Seek to the beginning of the BytesIO stream
                csv_bytes.seek(0)

                # Base64 encode the CSV data
                csv_base64 = base64.b64encode(csv_bytes.read().encode()).decode()

                # Construct the download link
                csv_href = f"data:text/csv;base64,{csv_base64}"

                # Specify the filename for the download
                filename = "data.csv"

                # Return the download link and filename
                return dcc.send_data_frame(self.df_sorted.to_csv, filename=filename)
