# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class List(Component):
    """A List component.
Display ordered or unordered list, see: https://mantine.dev/core/list/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    <List.Item /> components only.

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- center (boolean; optional):
    Center items with icon.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- icon (boolean | number | string | dict | list; optional):
    Icon that should replace list item dot.

- listStyleType (string; optional):
    List style type such as circle, square, upper-roman, lower-alpha,
    etc.

- size (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Font size from theme or number to set value in px.

- spacing (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Spacing between items from theme or number to set value in px.

- style (dict; optional):
    Inline style override.

- type (a value equal to: "ordered", "unordered"; optional):
    List type: ol or ul.

- withPadding (boolean; optional):
    Include padding-left to offset list from main content."""
    @_explicitize_args
    def __init__(self, children=None, center=Component.UNDEFINED, class_name=Component.UNDEFINED, icon=Component.UNDEFINED, id=Component.UNDEFINED, listStyleType=Component.UNDEFINED, size=Component.UNDEFINED, spacing=Component.UNDEFINED, style=Component.UNDEFINED, type=Component.UNDEFINED, withPadding=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'id', 'center', 'class_name', 'icon', 'listStyleType', 'size', 'spacing', 'style', 'type', 'withPadding']
        self._type = 'List'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'id', 'center', 'class_name', 'icon', 'listStyleType', 'size', 'spacing', 'style', 'type', 'withPadding']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(List, self).__init__(children=children, **args)
