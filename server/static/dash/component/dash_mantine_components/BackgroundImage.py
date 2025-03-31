# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class BackgroundImage(Component):
    """A BackgroundImage component.
Use when you need to display image below any content. For more information, see: https://mantine.dev/core/image/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    Contents.

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- radius (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Predefined border-radius value from theme.radius or number for
    border-radius in px.

- src (string; optional):
    Image src."""
    @_explicitize_args
    def __init__(self, children=None, id=Component.UNDEFINED, radius=Component.UNDEFINED, src=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'id', 'radius', 'src']
        self._type = 'BackgroundImage'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'id', 'radius', 'src']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(BackgroundImage, self).__init__(children=children, **args)
