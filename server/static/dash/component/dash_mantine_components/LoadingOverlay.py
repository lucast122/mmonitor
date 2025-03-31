# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class LoadingOverlay(Component):
    """A LoadingOverlay component.
Similar to dcc.Loading, overlay over given container with centered Loader from Mantine Theme. For more information, see: https://mantine.dev/core/loading-overlay/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    Content.

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- loader (a list of or a singular dash component, string or number; optional):
    Provide custom loader.

- loaderProps (dict; optional):
    Loader component props.

    `loaderProps` is a dict with keys:

    - color (a value equal to: "dark", "gray", "red", "pink", "grape", "violet", "indigo", "blue", "cyan", "teal", "green", "lime", "yellow", "orange"; optional):
        Loader color from theme.

    - size (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
        Defines width of loader.

    - variant (a value equal to: "bars", "oval", "dots"; optional):
        Loader appearance.

- loading_state (dict; optional):
    Object that holds the loading state object coming from
    dash-renderer.

    `loading_state` is a dict with keys:

    - component_name (string; optional):
        Holds the name of the component that is loading.

    - is_loading (boolean; optional):
        Determines if the component is loading or not.

    - prop_name (string; optional):
        Holds which property is loading.

- overlayColor (string; optional):
    Sets overlay color, defaults to theme.white in light theme and to
    theme.colors.dark[5] in dark theme.

- overlayOpacity (number; optional):
    Sets overlay opacity.

- radius (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Value from theme.radius or number to set border-radius in px.

- style (dict; optional):
    Inline style override.

- transitionDuration (number; optional):
    Appear and disappear animation duration.

- zIndex (number; optional):
    Loading overlay z-index."""
    @_explicitize_args
    def __init__(self, children=None, class_name=Component.UNDEFINED, id=Component.UNDEFINED, loader=Component.UNDEFINED, loaderProps=Component.UNDEFINED, loading_state=Component.UNDEFINED, overlayColor=Component.UNDEFINED, overlayOpacity=Component.UNDEFINED, radius=Component.UNDEFINED, style=Component.UNDEFINED, transitionDuration=Component.UNDEFINED, zIndex=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'id', 'class_name', 'loader', 'loaderProps', 'loading_state', 'overlayColor', 'overlayOpacity', 'radius', 'style', 'transitionDuration', 'zIndex']
        self._type = 'LoadingOverlay'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'id', 'class_name', 'loader', 'loaderProps', 'loading_state', 'overlayColor', 'overlayOpacity', 'radius', 'style', 'transitionDuration', 'zIndex']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(LoadingOverlay, self).__init__(children=children, **args)
