# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class Kbd(Component):
    """A Kbd component.
Display keyboard button or keys combination. For more information, see: https://mantine.dev/core/kbd/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    Keyboard key.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- style (dict; optional):
    Inline style override."""
    @_explicitize_args
    def __init__(self, children=None, class_name=Component.UNDEFINED, style=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'class_name', 'style']
        self._type = 'Kbd'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'class_name', 'style']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(Kbd, self).__init__(children=children, **args)
