# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class ColorPicker(Component):
    """A ColorPicker component.
Alternative to Select and RadioGroup. For more information, see: https://mantine.dev/core/ColorPicker/

Keyword arguments:

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- focusable (boolean; optional):
    Should interactive elements be focusable.

- format (a value equal to: "hex", "rgba", "rgb", "hsl", "hsla"; optional):
    Color format.

- fullWidth (boolean; optional):
    Force picker to take 100% width of its container.

- persisted_props (list of a value equal to: "value"s; default ["value"]):
    Properties whose user interactions will persist after refreshing
    the component or the page.

- persistence (boolean | string | number; optional):
    Used to allow user interactions in this component to be persisted
    when the component - or the page - is refreshed. If `persisted` is
    truthy and hasn't changed from its previous value, a `value` that
    the user has changed while using the app will keep that change, as
    long as the new `value` also matches what was given originally.
    Used in conjunction with `persistence_type`.

- persistence_type (a value equal to: "local", "session", "memory"; default "local"):
    Where persisted user changes will be stored: memory: only kept in
    memory, reset on page refresh. local: window.localStorage, data is
    kept after the browser quit. session: window.sessionStorage, data
    is cleared once the browser quit.

- size (a value equal to: "xs", "sm", "md", "lg", "xl"; optional):
    Predefined component size.

- style (dict; optional):
    Inline style override.

- swatches (list of strings; optional):
    Predefined colors.

- swatchesPerRow (number; optional):
    Number of swatches displayed in one row.

- value (string; optional):
    Controlled component value.

- withPicker (boolean; optional):
    Set to False to display swatches only."""
    @_explicitize_args
    def __init__(self, class_name=Component.UNDEFINED, focusable=Component.UNDEFINED, format=Component.UNDEFINED, fullWidth=Component.UNDEFINED, id=Component.UNDEFINED, persistence=Component.UNDEFINED, persisted_props=Component.UNDEFINED, persistence_type=Component.UNDEFINED, size=Component.UNDEFINED, style=Component.UNDEFINED, swatches=Component.UNDEFINED, swatchesPerRow=Component.UNDEFINED, value=Component.UNDEFINED, withPicker=Component.UNDEFINED, **kwargs):
        self._prop_names = ['id', 'class_name', 'focusable', 'format', 'fullWidth', 'persisted_props', 'persistence', 'persistence_type', 'size', 'style', 'swatches', 'swatchesPerRow', 'value', 'withPicker']
        self._type = 'ColorPicker'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['id', 'class_name', 'focusable', 'format', 'fullWidth', 'persisted_props', 'persistence', 'persistence_type', 'size', 'style', 'swatches', 'swatchesPerRow', 'value', 'withPicker']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(ColorPicker, self).__init__(**args)
