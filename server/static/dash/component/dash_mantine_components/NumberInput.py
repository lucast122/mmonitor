# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class NumberInput(Component):
    """A NumberInput component.
Capture number input from user. For more information, see: https://mantine.dev/core/number-input/

Keyword arguments:

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- decimalSeparator (string; optional):
    The decimal separator.

- description (boolean | number | string | dict | list; optional):
    Input description, displayed after label.

- disabled (boolean; optional):
    The component can show it is currently unable to be interacted
    with.

- error (boolean | number | string | dict | list; optional):
    Displays error message after input.

- hideControls (boolean; optional):
    Removes increment/decrement controls.

- icon (boolean | number | string | dict | list; optional):
    Adds icon on the left side of input.

- iconWidth (number; optional):
    Width of icon section in px.

- label (boolean | number | string | dict | list; optional):
    Input label, displayed before input.

- max (number; optional):
    Maximum possible value.

- min (number; optional):
    Minimal possible value.

- noClampOnBlur (boolean; optional):
    Prevent value clamp on blur.

- precision (number; optional):
    Amount of digits after the decimal point.

- radius (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Input border-radius from theme or number to set border-radius in
    px.

- required (boolean; optional):
    Adds red asterisk on the right side of label.

- rightSection (boolean | number | string | dict | list; optional):
    Right section of input, similar to icon but on the right.

- rightSectionWidth (number; optional):
    Width of right section, is used to calculate input padding-right.

- size (a value equal to: "xs", "sm", "md", "lg", "xl"; optional):
    Input size.

- step (number; optional):
    Number by which value will be incremented/decremented with
    controls and up/down arrows.

- stepHoldDelay (number; optional):
    Initial delay in milliseconds before stepping the value.

- stepHoldInterval (number; optional):
    Delay before stepping the value. Number of milliseconds.

- style (dict; optional):
    Inline style override.

- value (a list of or a singular dash component, string or number; optional):
    Input value.

- variant (a value equal to: "default", "filled", "unstyled", "headless"; optional):
    Defines input appearance, defaults to default in light color
    scheme and filled in dark."""
    @_explicitize_args
    def __init__(self, class_name=Component.UNDEFINED, decimalSeparator=Component.UNDEFINED, description=Component.UNDEFINED, disabled=Component.UNDEFINED, error=Component.UNDEFINED, hideControls=Component.UNDEFINED, icon=Component.UNDEFINED, iconWidth=Component.UNDEFINED, id=Component.UNDEFINED, label=Component.UNDEFINED, max=Component.UNDEFINED, min=Component.UNDEFINED, noClampOnBlur=Component.UNDEFINED, precision=Component.UNDEFINED, radius=Component.UNDEFINED, required=Component.UNDEFINED, rightSection=Component.UNDEFINED, rightSectionWidth=Component.UNDEFINED, size=Component.UNDEFINED, step=Component.UNDEFINED, style=Component.UNDEFINED, stepHoldDelay=Component.UNDEFINED, stepHoldInterval=Component.UNDEFINED, value=Component.UNDEFINED, variant=Component.UNDEFINED, **kwargs):
        self._prop_names = ['id', 'class_name', 'decimalSeparator', 'description', 'disabled', 'error', 'hideControls', 'icon', 'iconWidth', 'label', 'max', 'min', 'noClampOnBlur', 'precision', 'radius', 'required', 'rightSection', 'rightSectionWidth', 'size', 'step', 'stepHoldDelay', 'stepHoldInterval', 'style', 'value', 'variant']
        self._type = 'NumberInput'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['id', 'class_name', 'decimalSeparator', 'description', 'disabled', 'error', 'hideControls', 'icon', 'iconWidth', 'label', 'max', 'min', 'noClampOnBlur', 'precision', 'radius', 'required', 'rightSection', 'rightSectionWidth', 'size', 'step', 'stepHoldDelay', 'stepHoldInterval', 'style', 'value', 'variant']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(NumberInput, self).__init__(**args)
