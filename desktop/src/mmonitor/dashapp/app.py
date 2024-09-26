from dash import Dash

"""
Configure app/server instances of dash.
"""

# meta_tags are required for the app layout to be mobile responsive
if 'app' not in locals():
    app = Dash(
        __name__,
        suppress_callback_exceptions=True,
        meta_tags=[{'name': 'viewport', 'content': 'width=device-width, initial-scale=1.0'}]
    )
    server = app.server
