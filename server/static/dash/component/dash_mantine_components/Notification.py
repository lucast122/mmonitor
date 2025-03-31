# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class Notification(Component):
    """A Notification component.
Show dynamic notifications and alerts to user, part of notifications system. For more information, see: https://mantine.dev/others/notifications/

Keyword arguments:

- id (string; required):
    The ID of this component, used to identify dash components in
    callbacks.

- action (a value equal to: "show", "update", "hide"; required):
    Action.

- autoClose (boolean | number; optional):
    Whether autoclose and if True then duration.

- color (a value equal to: "dark", "gray", "red", "pink", "grape", "violet", "indigo", "blue", "cyan", "teal", "green", "lime", "yellow", "orange"; optional):
    Notification line or icon color.

- disallowClose (boolean; optional):
    Removes close button.

- icon (boolean | number | string | dict | list; optional):
    Notification icon, replaces color line.

- loading (boolean; optional):
    Replaces colored line or icon with Loader component.

- message (boolean | number | string | dict | list; optional):
    Notification Body.

- radius (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Radius from theme.radius, or number to set border-radius in px.

- title (boolean | number | string | dict | list; optional):
    Notification title, displayed before body."""
    @_explicitize_args
    def __init__(self, action=Component.REQUIRED, autoClose=Component.UNDEFINED, color=Component.UNDEFINED, disallowClose=Component.UNDEFINED, icon=Component.UNDEFINED, id=Component.REQUIRED, loading=Component.UNDEFINED, message=Component.UNDEFINED, radius=Component.UNDEFINED, title=Component.UNDEFINED, **kwargs):
        self._prop_names = ['id', 'action', 'autoClose', 'color', 'disallowClose', 'icon', 'loading', 'message', 'radius', 'title']
        self._type = 'Notification'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['id', 'action', 'autoClose', 'color', 'disallowClose', 'icon', 'loading', 'message', 'radius', 'title']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in ['action', 'id']:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(Notification, self).__init__(**args)
