# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class TimelineItem(Component):
    """A TimelineItem component.
Display list of events in chronological order. For more information, see: https://mantine.dev/core/timeline/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    React node that will be rendered after title.

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- align (a value equal to: "left", "right"; optional):
    Line and bullet position relative to item content, controlled by
    Timeline component.

- bullet (boolean | number | string | dict | list; optional):
    Bullet size in px.

- bulletSize (number; optional):
    Bullet size in px.

- color (a value equal to: "dark", "gray", "red", "pink", "grape", "violet", "indigo", "blue", "cyan", "teal", "green", "lime", "yellow", "orange"; optional):
    Active color from theme.

- lineVariant (a value equal to: "dashed", "dotted", "solid"; optional):
    Line border style.

- lineWidth (number; optional):
    Line border width in px, controlled by Timeline component.

- radius (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Radius from theme.radius, or number to set border-radius in px.

- title (boolean | number | string | dict | list; optional):
    Item title, rendered next to bullet."""
    @_explicitize_args
    def __init__(self, children=None, align=Component.UNDEFINED, bullet=Component.UNDEFINED, bulletSize=Component.UNDEFINED, color=Component.UNDEFINED, id=Component.UNDEFINED, lineVariant=Component.UNDEFINED, lineWidth=Component.UNDEFINED, radius=Component.UNDEFINED, title=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'id', 'align', 'bullet', 'bulletSize', 'color', 'lineVariant', 'lineWidth', 'radius', 'title']
        self._type = 'TimelineItem'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'id', 'align', 'bullet', 'bulletSize', 'color', 'lineVariant', 'lineWidth', 'radius', 'title']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(TimelineItem, self).__init__(children=children, **args)
