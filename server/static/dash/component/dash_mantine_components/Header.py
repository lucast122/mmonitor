# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class Header(Component):
    """A Header component.
Header. For more information, see: https://mantine.dev/core/app-shell/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    Header content.

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- fixed (boolean; optional):
    Changes position to fixed, controlled by AppShell component if
    rendered inside.

- height (string | number; optional):
    Header height.

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

- position (dict; optional):
    Control top, left, right or bottom position values, controlled by
    AppShell component if rendered inside.

    `position` is a dict with keys:

    - bottom (number; optional)

    - left (number; optional)

    - right (number; optional)

    - top (number; optional)

- pr (string | number; optional):
    Padding props.

- pt (string | number; optional):
    Padding props.

- px (string | number; optional):
    Padding props.

- py (string | number; optional):
    Padding props.

- style (dict; optional):
    Inline style override.

- zIndex (number; optional):
    z-index."""
    @_explicitize_args
    def __init__(self, children=None, class_name=Component.UNDEFINED, fixed=Component.UNDEFINED, height=Component.UNDEFINED, id=Component.UNDEFINED, position=Component.UNDEFINED, style=Component.UNDEFINED, zIndex=Component.UNDEFINED, m=Component.UNDEFINED, my=Component.UNDEFINED, mx=Component.UNDEFINED, mt=Component.UNDEFINED, mb=Component.UNDEFINED, ml=Component.UNDEFINED, mr=Component.UNDEFINED, p=Component.UNDEFINED, py=Component.UNDEFINED, px=Component.UNDEFINED, pt=Component.UNDEFINED, pb=Component.UNDEFINED, pl=Component.UNDEFINED, pr=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'id', 'class_name', 'fixed', 'height', 'm', 'mb', 'ml', 'mr', 'mt', 'mx', 'my', 'p', 'pb', 'pl', 'position', 'pr', 'pt', 'px', 'py', 'style', 'zIndex']
        self._type = 'Header'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'id', 'class_name', 'fixed', 'height', 'm', 'mb', 'ml', 'mr', 'mt', 'mx', 'my', 'p', 'pb', 'pl', 'position', 'pr', 'pt', 'px', 'py', 'style', 'zIndex']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(Header, self).__init__(children=children, **args)
