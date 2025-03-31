# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class Code(Component):
    """A Code component.
Inline or block code without syntax highlight. For more information, see: https://mantine.dev/core/code/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    Code content.

- block (boolean; optional):
    True for code block, False for inline code.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- color (a value equal to: "dark", "gray", "red", "pink", "grape", "violet", "indigo", "blue", "cyan", "teal", "green", "lime", "yellow", "orange"; optional):
    Code color and background from theme, defaults to gray in light
    theme and to dark in dark theme.

- style (dict; optional):
    Inline style override."""
    @_explicitize_args
    def __init__(self, children=None, block=Component.UNDEFINED, class_name=Component.UNDEFINED, color=Component.UNDEFINED, style=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'block', 'class_name', 'color', 'style']
        self._type = 'Code'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'block', 'class_name', 'color', 'style']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(Code, self).__init__(children=children, **args)
