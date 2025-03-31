# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class MantineProvider(Component):
    """A MantineProvider component.
MantineProvider component allows you to change theme globally. It is not required if you decide to use default theme. For more information, see: https://mantine.dev/theming/mantine-provider/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    Children.

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- inherit (boolean; optional):
    inherit from one level up MantineProvider.

- styles (dict; optional):
    Styles.

- theme (dict; optional):
    Theme.

- withGlobalStyles (boolean; optional):
    With global styles.

- withNormalizeCSS (boolean; optional):
    Normalize CSS?."""
    @_explicitize_args
    def __init__(self, children=None, id=Component.UNDEFINED, inherit=Component.UNDEFINED, styles=Component.UNDEFINED, theme=Component.UNDEFINED, withNormalizeCSS=Component.UNDEFINED, withGlobalStyles=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'id', 'inherit', 'styles', 'theme', 'withGlobalStyles', 'withNormalizeCSS']
        self._type = 'MantineProvider'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'id', 'inherit', 'styles', 'theme', 'withGlobalStyles', 'withNormalizeCSS']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(MantineProvider, self).__init__(children=children, **args)
