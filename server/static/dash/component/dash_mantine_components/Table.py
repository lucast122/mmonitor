# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class Table(Component):
    """A Table component.
A simple table component. For more information, see: https://mantine.dev/core/table/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    Component children, specifically an HTML representation of the
    table.

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- captionSide (a value equal to: "bottom", "top"; optional):
    Table caption position.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- fontSize (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Sets font size of all text inside table.

- highlightOnHover (boolean; optional):
    If True row will have hover color.

- horizontalSpacing (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Horizontal cells spacing from theme.spacing or number to set value
    in px.

- striped (boolean; optional):
    If True every odd row of table will have gray background color.

- style (dict; optional):
    Inline style override.

- verticalSpacing (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Vertical cells spacing from theme.spacing or number to set value
    in px."""
    @_explicitize_args
    def __init__(self, children=None, captionSide=Component.UNDEFINED, class_name=Component.UNDEFINED, fontSize=Component.UNDEFINED, id=Component.UNDEFINED, highlightOnHover=Component.UNDEFINED, horizontalSpacing=Component.UNDEFINED, striped=Component.UNDEFINED, style=Component.UNDEFINED, verticalSpacing=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'id', 'captionSide', 'class_name', 'fontSize', 'highlightOnHover', 'horizontalSpacing', 'striped', 'style', 'verticalSpacing']
        self._type = 'Table'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'id', 'captionSide', 'class_name', 'fontSize', 'highlightOnHover', 'horizontalSpacing', 'striped', 'style', 'verticalSpacing']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(Table, self).__init__(children=children, **args)
