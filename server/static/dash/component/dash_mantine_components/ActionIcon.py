# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class ActionIcon(Component):
    """An ActionIcon component.
Icon ActionIcon to indicate secondary action. For more information, see: https://mantine.dev/core/action-icon/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    Primary content.

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- color (a value equal to: "dark", "gray", "red", "pink", "grape", "violet", "indigo", "blue", "cyan", "teal", "green", "lime", "yellow", "orange"; optional):
    ActionIcon color from theme.

- loaderProps (dict; optional):
    Props spread to Loader component.

    `loaderProps` is a dict with keys:

    - color (a value equal to: "dark", "gray", "red", "pink", "grape", "violet", "indigo", "blue", "cyan", "teal", "green", "lime", "yellow", "orange"; optional):
        Loader color from theme.

    - size (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
        Defines width of loader.

    - variant (a value equal to: "bars", "oval", "dots"; optional):
        Loader appearance.

- loading (boolean; default False):
    Indicate loading state.

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

- n_clicks (number; default 0):
    An integer that represents the number of times that this element
    has been clicked on.

- radius (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    ActionIcon border-radius from theme or number to set border-radius
    in px.

- size (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Predefined ActionIcon size.

- style (dict; optional):
    Inline style override.

- variant (a value equal to: "transparent", "default", "hover", "filled", "light", "gradient", "outline"; optional):
    Controls ActionIcon appearance."""
    @_explicitize_args
    def __init__(self, children=None, class_name=Component.UNDEFINED, color=Component.UNDEFINED, id=Component.UNDEFINED, loaderProps=Component.UNDEFINED, loading=Component.UNDEFINED, loading_state=Component.UNDEFINED, n_clicks=Component.UNDEFINED, radius=Component.UNDEFINED, size=Component.UNDEFINED, style=Component.UNDEFINED, variant=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'id', 'class_name', 'color', 'loaderProps', 'loading', 'loading_state', 'n_clicks', 'radius', 'size', 'style', 'variant']
        self._type = 'ActionIcon'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'id', 'class_name', 'color', 'loaderProps', 'loading', 'loading_state', 'n_clicks', 'radius', 'size', 'style', 'variant']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(ActionIcon, self).__init__(children=children, **args)
