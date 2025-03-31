# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class Container(Component):
    """A Container component.
Center content horizontally with predefined max-width. For more information, see: https://mantine.dev/core/container/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    Content that should be centered vertically and horizontally.

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- fluid (boolean; optional):
    If fluid is set to True, size prop is ignored and Container always
    take 100% of width.

- m (string | number; optional):
    Margin props.

- mb (string | number; optional):
    Margin props.

- ml (string | number; optional):
    Margin props.

- mr (string | number; optional):
    Margin props.

- mt (string | number; optional):
    Margin props.

- mx (string | number; optional):
    Margin props.

- my (string | number; optional):
    Margin props.

- p (string | number; optional):
    Padding props.

- pb (string | number; optional):
    Padding props.

- pl (string | number; optional):
    Padding props.

- pr (string | number; optional):
    Padding props.

- pt (string | number; optional):
    Padding props.

- px (string | number; optional):
    Padding props.

- py (string | number; optional):
    Padding props.

- size (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Predefined container max-width or number for max-width in px.

- style (dict; optional):
    Inline style override."""
    @_explicitize_args
    def __init__(self, children=None, class_name=Component.UNDEFINED, fluid=Component.UNDEFINED, id=Component.UNDEFINED, size=Component.UNDEFINED, style=Component.UNDEFINED, m=Component.UNDEFINED, my=Component.UNDEFINED, mx=Component.UNDEFINED, mt=Component.UNDEFINED, mb=Component.UNDEFINED, ml=Component.UNDEFINED, mr=Component.UNDEFINED, p=Component.UNDEFINED, py=Component.UNDEFINED, px=Component.UNDEFINED, pt=Component.UNDEFINED, pb=Component.UNDEFINED, pl=Component.UNDEFINED, pr=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'id', 'class_name', 'fluid', 'm', 'mb', 'ml', 'mr', 'mt', 'mx', 'my', 'p', 'pb', 'pl', 'pr', 'pt', 'px', 'py', 'size', 'style']
        self._type = 'Container'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'id', 'class_name', 'fluid', 'm', 'mb', 'ml', 'mr', 'mt', 'mx', 'my', 'p', 'pb', 'pl', 'pr', 'pt', 'px', 'py', 'size', 'style']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(Container, self).__init__(children=children, **args)
