# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class Alert(Component):
    """An Alert component.
Attract user attention with important static message. For more information, see: https://mantine.dev/core/alert/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    Alert message.

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- color (a value equal to: "dark", "gray", "red", "pink", "grape", "violet", "indigo", "blue", "cyan", "teal", "green", "lime", "yellow", "orange"; optional):
    Alert title and line colors from theme.

- duration (number; optional):
    Duration in milliseconds after which the Alert dismisses itself.

- hide (boolean; default False):
    Whether to hide the alert.

- icon (boolean | number | string | dict | list; optional):
    Icon displayed next to title.

- radius (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Radius from theme.radius, or number to set border-radius in px.

- style (dict; optional):
    Inline style override.

- title (boolean | number | string | dict | list; optional):
    Optional alert title.

- variant (a value equal to: "filled", "outline", "light"; default "light"):
    Controls Alert background, color and border styles.

- withCloseButton (boolean; optional):
    Display close button."""
    @_explicitize_args
    def __init__(self, children=None, class_name=Component.UNDEFINED, color=Component.UNDEFINED, duration=Component.UNDEFINED, hide=Component.UNDEFINED, icon=Component.UNDEFINED, id=Component.UNDEFINED, radius=Component.UNDEFINED, title=Component.UNDEFINED, style=Component.UNDEFINED, variant=Component.UNDEFINED, withCloseButton=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'id', 'class_name', 'color', 'duration', 'hide', 'icon', 'radius', 'style', 'title', 'variant', 'withCloseButton']
        self._type = 'Alert'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'id', 'class_name', 'color', 'duration', 'hide', 'icon', 'radius', 'style', 'title', 'variant', 'withCloseButton']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(Alert, self).__init__(children=children, **args)
