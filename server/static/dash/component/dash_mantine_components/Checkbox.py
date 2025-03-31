# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class Checkbox(Component):
    """A Checkbox component.
Capture boolean input from user. For more information, see: https://mantine.dev/core/checkbox/

Keyword arguments:

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- checked (boolean; default False):
    State of check box.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- color (a value equal to: "dark", "gray", "red", "pink", "grape", "violet", "indigo", "blue", "cyan", "teal", "green", "lime", "yellow", "orange"; optional):
    Checkbox color.

- disabled (boolean; optional):
    A checkbox can show it is currently unable to be interacted with.

- label (string; optional):
    Checkbox label.

- persisted_props (list of a value equal to: "checked"s; default ["checked"]):
    Properties whose user interactions will persist after refreshing
    the component or the page.

- persistence (boolean | string | number; optional):
    Used to allow user interactions in this component to be persisted
    when the component - or the page - is refreshed. If `persisted` is
    truthy and hasn't changed from its previous value, a `value` that
    the user has changed while using the app will keep that change, as
    long as the new `value` also matches what was given originally.
    Used in conjunction with `persistence_type`.

- persistence_type (a value equal to: "local", "session", "memory"; default "local"):
    Where persisted user changes will be stored: memory: only kept in
    memory, reset on page refresh. local: window.localStorage, data is
    kept after the browser quit. session: window.sessionStorage, data
    is cleared once the browser quit.

- radius (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Radius from theme.radius, or number to set border-radius in px.

- size (a value equal to: "xs", "sm", "md", "lg", "xl"; optional):
    Predefined label font-size and checkbox width and height in px.

- style (dict; optional):
    Inline style override.

- transitionDuration (number; optional):
    Check/uncheck transition duration, set to 0 to disable all
    transitions."""
    @_explicitize_args
    def __init__(self, checked=Component.UNDEFINED, class_name=Component.UNDEFINED, color=Component.UNDEFINED, disabled=Component.UNDEFINED, id=Component.UNDEFINED, label=Component.UNDEFINED, persistence=Component.UNDEFINED, persisted_props=Component.UNDEFINED, persistence_type=Component.UNDEFINED, radius=Component.UNDEFINED, size=Component.UNDEFINED, style=Component.UNDEFINED, transitionDuration=Component.UNDEFINED, **kwargs):
        self._prop_names = ['id', 'checked', 'class_name', 'color', 'disabled', 'label', 'persisted_props', 'persistence', 'persistence_type', 'radius', 'size', 'style', 'transitionDuration']
        self._type = 'Checkbox'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['id', 'checked', 'class_name', 'color', 'disabled', 'label', 'persisted_props', 'persistence', 'persistence_type', 'radius', 'size', 'style', 'transitionDuration']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(Checkbox, self).__init__(**args)
