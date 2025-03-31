# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class RadioGroup(Component):
    """A RadioGroup component.
Capture user feedback limited to small set of options. For more information, see: https://mantine.dev/core/radio-group/

Keyword arguments:

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- color (a value equal to: "dark", "gray", "red", "pink", "grape", "violet", "indigo", "blue", "cyan", "teal", "green", "lime", "yellow", "orange"; optional):
    Active radio color from theme.colors.

- data (list of dicts; optional):
    RadioGroup options.

    `data` is a list of dicts with keys:

    - label (string; required):
        The option's label.

    - value (string; required):
        Option's value.

- description (string; optional):
    Input description, displayed after label.

- error (string; optional):
    Displays error message after input.

- label (string; optional):
    Input label, displayed before input.

- orientation (a value equal to: "horizontal", "vertical"; optional):
    Horizontal or vertical orientation.

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

- required (boolean; optional):
    Adds red asterisk on the right side of label.

- size (a value equal to: "xs", "sm", "md", "lg", "xl"; optional):
    Predefined label fontSize, radio width, height and border-radius.

- spacing (a value equal to: "xs", "sm", "md", "lg", "xl"; optional):
    Spacing between radios in horizontal variant.

- style (dict; optional):
    Inline style override.

- value (string; optional):
    Value of currently selected radio (controlled)."""
    @_explicitize_args
    def __init__(self, class_name=Component.UNDEFINED, color=Component.UNDEFINED, data=Component.UNDEFINED, description=Component.UNDEFINED, error=Component.UNDEFINED, id=Component.UNDEFINED, label=Component.UNDEFINED, orientation=Component.UNDEFINED, persistence=Component.UNDEFINED, persisted_props=Component.UNDEFINED, persistence_type=Component.UNDEFINED, required=Component.UNDEFINED, size=Component.UNDEFINED, spacing=Component.UNDEFINED, style=Component.UNDEFINED, value=Component.UNDEFINED, **kwargs):
        self._prop_names = ['id', 'class_name', 'color', 'data', 'description', 'error', 'label', 'orientation', 'persisted_props', 'persistence', 'persistence_type', 'required', 'size', 'spacing', 'style', 'value']
        self._type = 'RadioGroup'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['id', 'class_name', 'color', 'data', 'description', 'error', 'label', 'orientation', 'persisted_props', 'persistence', 'persistence_type', 'required', 'size', 'spacing', 'style', 'value']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(RadioGroup, self).__init__(**args)
