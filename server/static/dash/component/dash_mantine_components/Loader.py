# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class Loader(Component):
    """A Loader component.
Indicate loading state. For more information, see: https://mantine.dev/core/loader/

Keyword arguments:

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- color (a value equal to: "dark", "gray", "red", "pink", "grape", "violet", "indigo", "blue", "cyan", "teal", "green", "lime", "yellow", "orange"; optional):
    Loader color from theme.

- size (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Defines width of loader.

- style (dict; optional):
    Inline style override.

- variant (a value equal to: "bars", "oval", "dots"; optional):
    Loader appearance."""
    @_explicitize_args
    def __init__(self, class_name=Component.UNDEFINED, id=Component.UNDEFINED, color=Component.UNDEFINED, size=Component.UNDEFINED, style=Component.UNDEFINED, variant=Component.UNDEFINED, **kwargs):
        self._prop_names = ['id', 'class_name', 'color', 'size', 'style', 'variant']
        self._type = 'Loader'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['id', 'class_name', 'color', 'size', 'style', 'variant']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(Loader, self).__init__(**args)
