# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class JsonInput(Component):
    """A JsonInput component.
Capture json data from user. For more information, see: https://mantine.dev/core/json-input/

Keyword arguments:

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- autosize (boolean; optional):
    If True textarea will grow with content until maxRows are reached.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- description (boolean | number | string | dict | list; optional):
    Input description, displayed after label [PropTypes.node].

- disabled (boolean; optional):
    Disabled input state.

- error (boolean | number | string | dict | list; optional):
    Displays error message after input [PropTypes.node].

- formatOnBlur (boolean; optional):
    Format json on blur.

- icon (boolean | number | string | dict | list; optional):
    Adds icon on the left side of input [PropTypes.node].

- iconWidth (number; optional):
    Width of icon section in px.

- label (boolean | number | string | dict | list; optional):
    Input label, displayed before input [PropTypes.node].

- maxRows (number; optional):
    Defines maxRows in autosize variant, not applicable to regular
    variant.

- minRows (number; optional):
    Defines minRows in autosize variant, not applicable to regular
    variant.

- placeholder (string; optional):
    Placeholder.

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

- validationError (boolean | number | string | dict | list; optional):
    Error message shown when json is not valid [PropTypes.node].

- value (string; optional):
    Value for controlled input.

- variant (a value equal to: "default", "filled", "unstyled", "headless"; optional):
    Defines input appearance, defaults to default in light color
    scheme and filled in dark."""
    @_explicitize_args
    def __init__(self, autosize=Component.UNDEFINED, class_name=Component.UNDEFINED, description=Component.UNDEFINED, disabled=Component.UNDEFINED, error=Component.UNDEFINED, formatOnBlur=Component.UNDEFINED, icon=Component.UNDEFINED, iconWidth=Component.UNDEFINED, id=Component.UNDEFINED, label=Component.UNDEFINED, maxRows=Component.UNDEFINED, minRows=Component.UNDEFINED, placeholder=Component.UNDEFINED, radius=Component.UNDEFINED, required=Component.UNDEFINED, rightSection=Component.UNDEFINED, rightSectionWidth=Component.UNDEFINED, size=Component.UNDEFINED, style=Component.UNDEFINED, validationError=Component.UNDEFINED, value=Component.UNDEFINED, variant=Component.UNDEFINED, **kwargs):
        self._prop_names = ['id', 'autosize', 'class_name', 'description', 'disabled', 'error', 'formatOnBlur', 'icon', 'iconWidth', 'label', 'maxRows', 'minRows', 'placeholder', 'radius', 'required', 'rightSection', 'rightSectionWidth', 'size', 'style', 'validationError', 'value', 'variant']
        self._type = 'JsonInput'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['id', 'autosize', 'class_name', 'description', 'disabled', 'error', 'formatOnBlur', 'icon', 'iconWidth', 'label', 'maxRows', 'minRows', 'placeholder', 'radius', 'required', 'rightSection', 'rightSectionWidth', 'size', 'style', 'validationError', 'value', 'variant']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(JsonInput, self).__init__(**args)
