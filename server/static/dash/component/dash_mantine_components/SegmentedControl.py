# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class SegmentedControl(Component):
    """A SegmentedControl component.
Horizontal control made of multiple segments, alternative to RadioGroup. For more information, see: https://mantine.dev/core/segmented-control/

Keyword arguments:

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- color (a value equal to: "dark", "gray", "red", "pink", "grape", "violet", "indigo", "blue", "cyan", "teal", "green", "lime", "yellow", "orange"; optional):
    Active control color from theme.colors, defaults to white in light
    color scheme and theme.colors.dark[9] in dark.

- data (list of dicts; required):
    Data based on which controls are rendered.

    `data` is a list of dicts with keys:

    - label (string; required):
        The option's label.

    - value (string; required):
        Option's value. | list of strings

- disabled (boolean; optional):
    Disabled input state.

- fullWidth (boolean; optional):
    True if component should have 100% width.

- orientation (a value equal to: "horizontal", "vertical"; optional):
    Display Vertically or horizontally.

- persisted_props (list of a value equal to: "value"s; default ["value"]):
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
    Border-radius from theme or number to set border-radius in px.

- size (a value equal to: "xs", "sm", "md", "lg", "xl"; optional):
    Controls font-size, paddings and height.

- style (dict; optional):
    Inline style override.

- transitionDuration (number; optional):
    Transition duration in ms, set to 0 to turn off transitions.

- transitionTimingFunction (string; optional):
    Transition timing function for all transitions, defaults to
    theme.transitionTimingFunction.

- value (string; optional):
    Current selected value."""
    @_explicitize_args
    def __init__(self, class_name=Component.UNDEFINED, color=Component.UNDEFINED, data=Component.REQUIRED, disabled=Component.UNDEFINED, fullWidth=Component.UNDEFINED, id=Component.UNDEFINED, orientation=Component.UNDEFINED, persistence=Component.UNDEFINED, persisted_props=Component.UNDEFINED, persistence_type=Component.UNDEFINED, radius=Component.UNDEFINED, size=Component.UNDEFINED, style=Component.UNDEFINED, transitionDuration=Component.UNDEFINED, transitionTimingFunction=Component.UNDEFINED, value=Component.UNDEFINED, **kwargs):
        self._prop_names = ['id', 'class_name', 'color', 'data', 'disabled', 'fullWidth', 'orientation', 'persisted_props', 'persistence', 'persistence_type', 'radius', 'size', 'style', 'transitionDuration', 'transitionTimingFunction', 'value']
        self._type = 'SegmentedControl'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['id', 'class_name', 'color', 'data', 'disabled', 'fullWidth', 'orientation', 'persisted_props', 'persistence', 'persistence_type', 'radius', 'size', 'style', 'transitionDuration', 'transitionTimingFunction', 'value']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in ['data']:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(SegmentedControl, self).__init__(**args)
