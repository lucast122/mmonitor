# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class Timeline(Component):
    """A Timeline component.
Display list of events in chronological order. For more information, see: https://mantine.dev/core/timeline/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    dmc.TimelineItem components only.

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- active (number; optional):
    Index of active element.

- align (a value equal to: "left", "right"; optional):
    Timeline alignment.

- bulletSize (number; optional):
    Bullet size in px.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- color (a value equal to: "dark", "gray", "red", "pink", "grape", "violet", "indigo", "blue", "cyan", "teal", "green", "lime", "yellow", "orange"; optional):
    Active color from theme.

- lineWidth (number; optional):
    Line width in px.

- radius (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Radius from theme.radius, or number to set border-radius in px.

- reverseActive (boolean; optional):
    Reverse active direction without reversing items."""
    @_explicitize_args
    def __init__(self, children=None, active=Component.UNDEFINED, align=Component.UNDEFINED, bulletSize=Component.UNDEFINED, class_name=Component.UNDEFINED, color=Component.UNDEFINED, id=Component.UNDEFINED, lineWidth=Component.UNDEFINED, radius=Component.UNDEFINED, reverseActive=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'id', 'active', 'align', 'bulletSize', 'class_name', 'color', 'lineWidth', 'radius', 'reverseActive']
        self._type = 'Timeline'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'id', 'active', 'align', 'bulletSize', 'class_name', 'color', 'lineWidth', 'radius', 'reverseActive']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(Timeline, self).__init__(children=children, **args)
