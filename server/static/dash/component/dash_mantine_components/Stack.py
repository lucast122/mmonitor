# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class Stack(Component):
    """A Stack component.
Compose elements and components in vertical flex container. For more information, see: https://mantine.dev/core/Stack/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    Primary content inside the stack.

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- align (a value equal to: "stretch", "center", "flex-end", "flex-start"; optional):
    Sets text-align css property.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- justify (a value equal to: "space-between", "space-around", "center", "flex-end", "flex-start"; optional):
    Set grid justify-content property.

- spacing (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Spacing between children.

- style (dict; optional):
    Inline style."""
    @_explicitize_args
    def __init__(self, children=None, align=Component.UNDEFINED, class_name=Component.UNDEFINED, id=Component.UNDEFINED, justify=Component.UNDEFINED, spacing=Component.UNDEFINED, style=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'id', 'align', 'class_name', 'justify', 'spacing', 'style']
        self._type = 'Stack'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'id', 'align', 'class_name', 'justify', 'spacing', 'style']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(Stack, self).__init__(children=children, **args)
