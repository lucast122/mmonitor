# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class TimeInput(Component):
    """A TimeInput component.
Capture time input from user. For more information, see: https://mantine.dev/dates/time-input/

Keyword arguments:

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- amPmPlaceholder (string; optional):
    Placeholder for am/pm input.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- clearable (boolean; optional):
    Allow to clear item.

- description (boolean | number | string | dict | list; optional):
    Input description, displayed after label [PropTypes.node].

- disabled (boolean; optional):
    Disabled input state.

- error (boolean | number | string | dict | list; optional):
    Displays error message after input [PropTypes.node].

- format (a value equal to: "12", "24"; optional):
    Time format.

- icon (boolean | number | string | dict | list; optional):
    Adds icon on the left side of input [PropTypes.node].

- iconWidth (number; optional):
    Width of icon section in px.

- label (boolean | number | string | dict | list; optional):
    Input label, displayed before input [PropTypes.node].

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
    Input border-radius from theme or number to set border-radius in
    px.

- required (boolean; optional):
    Adds red asterisk on the right side of label.

- rightSection (boolean | number | string | dict | list; optional):
    Right section of input, similar to icon but on the right
    [PropTypes.node].

- rightSectionWidth (number; optional):
    Width of right section, is used to calculate input padding-right.

- size (a value equal to: "xs", "sm", "md", "lg", "xl"; optional):
    Input size.

- style (dict; optional):
    Inline style override.

- timePlaceholder (string; optional):
    Placeholder for hours/minutes/seconds inputs.

- value (string; optional):
    Value for controlled input.

- variant (a value equal to: "default", "filled", "unstyled", "headless"; optional):
    Defines input appearance, defaults to default in light color
    scheme and filled in dark.

- withSeconds (boolean; optional):
    Display seconds input."""
    @_explicitize_args
    def __init__(self, amPmPlaceholder=Component.UNDEFINED, class_name=Component.UNDEFINED, clearable=Component.UNDEFINED, description=Component.UNDEFINED, disabled=Component.UNDEFINED, error=Component.UNDEFINED, format=Component.UNDEFINED, icon=Component.UNDEFINED, iconWidth=Component.UNDEFINED, id=Component.UNDEFINED, label=Component.UNDEFINED, persistence=Component.UNDEFINED, persisted_props=Component.UNDEFINED, persistence_type=Component.UNDEFINED, radius=Component.UNDEFINED, required=Component.UNDEFINED, rightSection=Component.UNDEFINED, rightSectionWidth=Component.UNDEFINED, size=Component.UNDEFINED, style=Component.UNDEFINED, timePlaceholder=Component.UNDEFINED, value=Component.UNDEFINED, variant=Component.UNDEFINED, withSeconds=Component.UNDEFINED, **kwargs):
        self._prop_names = ['id', 'amPmPlaceholder', 'class_name', 'clearable', 'description', 'disabled', 'error', 'format', 'icon', 'iconWidth', 'label', 'persisted_props', 'persistence', 'persistence_type', 'radius', 'required', 'rightSection', 'rightSectionWidth', 'size', 'style', 'timePlaceholder', 'value', 'variant', 'withSeconds']
        self._type = 'TimeInput'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['id', 'amPmPlaceholder', 'class_name', 'clearable', 'description', 'disabled', 'error', 'format', 'icon', 'iconWidth', 'label', 'persisted_props', 'persistence', 'persistence_type', 'radius', 'required', 'rightSection', 'rightSectionWidth', 'size', 'style', 'timePlaceholder', 'value', 'variant', 'withSeconds']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(TimeInput, self).__init__(**args)
