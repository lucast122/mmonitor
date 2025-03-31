# Register apps with Django Dash
for app_name, app_instance in app_registry.items():
    # Set suppress_callback_exceptions to True for all apps to avoid console errors
    app_instance.app._suppress_callback_exceptions = True
    
    # Log app registration
    print(f"Registering app: {app_name}")
    
    # Add error handling middleware if not already present
    if not hasattr(app_instance.app, 'error_handler'):
        @app_instance.app.server.errorhandler(Exception)
        def handle_error(e):
            print(f"Error in app {app_name}: {str(e)}")
            return f"An error occurred in the {app_name} app: {str(e)}", 500
        app_instance.app.error_handler = True

# Special handling for main app
main_app._suppress_callback_exceptions = True

# Add a callback to handle global errors
@main_app.callback(
    Output("global-error-container", "children"),
    [Input("url", "pathname")],
    prevent_initial_call=True
)
def handle_global_errors(pathname):
    # This callback is just to ensure the error container exists
    # and to demonstrate global error handling
    return [] 