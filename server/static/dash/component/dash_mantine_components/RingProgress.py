# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class RingProgress(Component):
    """A RingProgress component.
Give user feedback for status of the task with circle diagram. For more information, see: https://mantine.dev/core/ring-progress/

Keyword arguments:

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- label (boolean | number | string | dict | list; optional):
    Label displayed in the center of the ring [PropTypes.node].

- roundCaps (boolean; optional):
    Sets whether the edges of the progress circle are rounded.

- sections (list of dicts; required):
    Ring sections.

    `sections` is a list of dicts with keys:

    - color (a value equal to: "dark", "gray", "red", "pink", "grape", "violet", "indigo", "blue", "cyan", "teal", "green", "lime", "yellow", "orange" | string; required)

    - value (number; required)

- size (number; optional):
    Width and height of the progress ring in px.

- style (dict; optional):
    Inline style override.

- thickness (number; optional):
    Ring thickness."""
    @_explicitize_args
    def __init__(self, class_name=Component.UNDEFINED, id=Component.UNDEFINED, label=Component.UNDEFINED, roundCaps=Component.UNDEFINED, sections=Component.REQUIRED, size=Component.UNDEFINED, style=Component.UNDEFINED, thickness=Component.UNDEFINED, **kwargs):
        self._prop_names = ['id', 'class_name', 'label', 'roundCaps', 'sections', 'size', 'style', 'thickness']
        self._type = 'RingProgress'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['id', 'class_name', 'label', 'roundCaps', 'sections', 'size', 'style', 'thickness']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in ['sections']:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(RingProgress, self).__init__(**args)
