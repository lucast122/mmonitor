import base64
import glob
from pathlib import Path

import dash
from dash import dcc
from dash import html
from mmonitor.dashapp.app import app
from mmonitor.dashapp.base_app import BaseApp


class Kegg(BaseApp):

    def __init__(self):
        self._init_layout()
        self._init_callbacks()

    def _init_layout(self) -> None:
        self._static_image_route = "/static/"
        base_path = Path(__file__)
        base_path = Path(base_path.parent.absolute().parent.absolute().parent.absolute().parent.absolute())
        self._samples = glob.glob(f"{base_path}/resources/pipeline_out/*/")
        header = html.H1("KEGG pathways in metagenome")
        sample_dropdown = dcc.Dropdown(id="sample-dropdown", options=[{'label': t, 'value': t} for t in self._samples],
                                       value=self._samples[0],style={"max-width": "70%", "height": "auto"})

        pathway_dropdown = dcc.Dropdown(id="pathway-dropdown",style={"max-width": "70%", "height": "auto"})
        dd = html.Div(id='dd-output-container')

        CONTENT_STYLE = {

            "margin-right": "2rem",
            "padding": "2rem 1rem",
            'margin-bottom': '200px',
            'font-size': '15px'
        }


        image = html.Img(
            id='kegg-plot',
            alt='Please select a kegg plot',
            style={'display': 'block', 'margin-right': 'auto', 'padding': '20px', "max-width": "99%", "height": "auto"})



        container = html.Div([header, sample_dropdown, pathway_dropdown, dd, image],style=CONTENT_STYLE)

        self._layout = container

    def _init_callbacks(self) -> None:
        @app.callback(
            dash.dependencies.Output('pathway-dropdown', 'options'),
            [dash.dependencies.Input('sample-dropdown', 'value')]
        )
        def update_pathway_dropdown(sample):
            return [{'label': i.split("/")[-1].removesuffix(".png"), 'value': i} for i in glob.glob(f"{sample}/*.png")]

        @app.callback(
            dash.dependencies.Output('kegg-plot', 'src'),
            dash.dependencies.Input('pathway-dropdown', 'value')
        )
        def update_image_src(image_path):
            print('current image_path = {}'.format(image_path))
            encoded_image = base64.b64encode(open(image_path, 'rb').read())
            return 'data:image/png;base64,{}'.format(encoded_image.decode())

            return value
