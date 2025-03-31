# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class Tooltip(Component):
    """A Tooltip component.
Renders tooltip at given element on mouse over or any other event. For more information, see: https://mantine.dev/core/tooltip/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    Any react node that should trigger tooltip.

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- arrowDistance (number; optional):
    Arrow distance to the left/right * arrowSize.

- arrowSize (number; optional):
    Arrow size in px.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- closeDelay (number; optional):
    Close delay in ms, 0 to disable delay.

- color (a value equal to: "dark", "gray", "red", "pink", "grape", "violet", "indigo", "blue", "cyan", "teal", "green", "lime", "yellow", "orange"; optional):
    Any color from theme.colors, defaults to gray in light color
    scheme and dark in dark colors scheme.

- disabled (boolean; optional):
    True to disable tooltip.

- exitTransitionDuration (number; optional):
    Unmount transition duration in ms.

- gutter (number; optional):
    Spacing between element and popper in px.

- label (boolean | number | string | dict | list; required):
    Tooltip content.

- openDelay (number; optional):
    Open delay in ms, 0 to disable delay.

- opened (boolean; optional):
    Tooltip opened state for controlled variant.

- placement (a value equal to: "center", "end", "start"; optional):
    Placement relative to reference element.

- position (a value equal to: "bottom", "left", "right", "top"; optional):
    Position relative to reference element.

- radius (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Radius from theme.radius, or number to set border-radius in px.

- style (dict; optional):
    Inline style override.

- transition (a value equal to: "fade", "skew-up", "skew-down", "rotate-right", "rotate-left", "slide-down", "slide-up", "slide-right", "slide-left", "scale-y", "scale-x", "scale", "pop", "pop-top-left", "pop-top-right", "pop-bottom-left", "pop-bottom-right"; optional):
    Customize mount/unmount transition.

- transitionDuration (number; optional):
    Mount transition duration in ms.

- transitionTimingFunction (string; optional):
    Mount/unmount transition timing function, defaults to
    theme.transitionTimingFunction.

- width (number | a value equal to: "auto"; optional):
    Tooltip width in px or auto.

- withArrow (boolean; optional):
    Renders arrow if True.

- wrapLines (boolean; optional):
    Allow multiline tooltip content.

- zIndex (number; optional):
    Popper z-index."""
    @_explicitize_args
    def __init__(self, children=None, arrowDistance=Component.UNDEFINED, arrowSize=Component.UNDEFINED, class_name=Component.UNDEFINED, closeDelay=Component.UNDEFINED, color=Component.UNDEFINED, disabled=Component.UNDEFINED, exitTransitionDuration=Component.UNDEFINED, gutter=Component.UNDEFINED, id=Component.UNDEFINED, label=Component.REQUIRED, openDelay=Component.UNDEFINED, opened=Component.UNDEFINED, placement=Component.UNDEFINED, position=Component.UNDEFINED, radius=Component.UNDEFINED, style=Component.UNDEFINED, transition=Component.UNDEFINED, transitionDuration=Component.UNDEFINED, transitionTimingFunction=Component.UNDEFINED, width=Component.UNDEFINED, withArrow=Component.UNDEFINED, wrapLines=Component.UNDEFINED, zIndex=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'id', 'arrowDistance', 'arrowSize', 'class_name', 'closeDelay', 'color', 'disabled', 'exitTransitionDuration', 'gutter', 'label', 'openDelay', 'opened', 'placement', 'position', 'radius', 'style', 'transition', 'transitionDuration', 'transitionTimingFunction', 'width', 'withArrow', 'wrapLines', 'zIndex']
        self._type = 'Tooltip'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'id', 'arrowDistance', 'arrowSize', 'class_name', 'closeDelay', 'color', 'disabled', 'exitTransitionDuration', 'gutter', 'label', 'openDelay', 'opened', 'placement', 'position', 'radius', 'style', 'transition', 'transitionDuration', 'transitionTimingFunction', 'width', 'withArrow', 'wrapLines', 'zIndex']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in ['label']:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(Tooltip, self).__init__(children=children, **args)
