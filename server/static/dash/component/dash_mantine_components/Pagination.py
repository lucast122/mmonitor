# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class Pagination(Component):
    """A Pagination component.
Display active page and navigate between multiple pages. For more information, see: https://mantine.dev/core/pagination/

Keyword arguments:

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- align (a value equal to: "stretch", "center", "flex-end", "flex-start"; optional):
    Defines align-items css property.

- boundaries (number; optional):
    Amount of elements visible on left/right edges.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- color (a value equal to: "dark", "gray", "red", "pink", "grape", "violet", "indigo", "blue", "cyan", "teal", "green", "lime", "yellow", "orange"; optional):
    Active item color from theme, defaults to theme.primaryColor.

- direction (a value equal to: "row", "column"; optional):
    Defines flex-direction property, row for horizontal, column for
    vertical.

- grow (boolean; optional):
    Defines flex-grow property for each element, True -> 1, False ->
    0.

- noWrap (boolean; optional):
    Defined flex-wrap property.

- page (number; optional):
    Controlled active page number.

- position (a value equal to: "right", "center", "left", "apart"; optional):
    Defines justify-content property.

- radius (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Predefined item radius or number to set border-radius in px.

- siblings (number; optional):
    Siblings amount on left/right side of selected page.

- size (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Predefined item size or number to set width and height in px.

- spacing (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Spacing between items from theme or number to set value in px,
    defaults to theme.spacing.xs / 2.

- style (dict; optional):
    Inline style.

- total (number; required):
    Total amount of pages.

- withControls (boolean; optional):
    Show/hide prev/next controls.

- withEdges (boolean; optional):
    Show/hide jump to start/end controls."""
    @_explicitize_args
    def __init__(self, align=Component.UNDEFINED, boundaries=Component.UNDEFINED, class_name=Component.UNDEFINED, color=Component.UNDEFINED, direction=Component.UNDEFINED, grow=Component.UNDEFINED, id=Component.UNDEFINED, noWrap=Component.UNDEFINED, page=Component.UNDEFINED, position=Component.UNDEFINED, radius=Component.UNDEFINED, siblings=Component.UNDEFINED, size=Component.UNDEFINED, spacing=Component.UNDEFINED, style=Component.UNDEFINED, total=Component.REQUIRED, withControls=Component.UNDEFINED, withEdges=Component.UNDEFINED, **kwargs):
        self._prop_names = ['id', 'align', 'boundaries', 'class_name', 'color', 'direction', 'grow', 'noWrap', 'page', 'position', 'radius', 'siblings', 'size', 'spacing', 'style', 'total', 'withControls', 'withEdges']
        self._type = 'Pagination'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['id', 'align', 'boundaries', 'class_name', 'color', 'direction', 'grow', 'noWrap', 'page', 'position', 'radius', 'siblings', 'size', 'spacing', 'style', 'total', 'withControls', 'withEdges']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in ['total']:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(Pagination, self).__init__(**args)
