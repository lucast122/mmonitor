# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class Modal(Component):
    """A Modal component.
Modal with optional header. For more information, see: https://mantine.dev/core/modal/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    Content that should be centered vertically and horizontally.

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- centered (boolean; optional):
    Controls if modal should be centered.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- closeOnClickOutside (boolean; optional):
    Should modal be closed when outside click was registered?.

- closeOnEscape (boolean; optional):
    Should modal be closed when escape is pressed?.

- opened (boolean; default False):
    Mounts modal if True.

- overflow (a value equal to: "inside", "outside"; optional):
    Control vertical overflow behavior.

- overlayColor (string; optional):
    Overlay below modal color, defaults to theme.black in light theme
    and to theme.colors.dark[9] in dark theme.

- overlayOpacity (number; optional):
    Overlay below modal opacity, defaults to 0.75 in light theme and
    to 0.85 in dark theme.

- padding (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Modal padding from theme or number value for padding in px.

- radius (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Modal padding from theme or number value for padding in px.

- shadow (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Modal shadow from theme or css value.

- size (a value equal to: "xs", "sm", "md", "lg", "xl" | number | string; optional):
    Modal body width.

- style (dict; optional):
    Inline style override.

- title (boolean | number | string | dict | list; optional):
    Modal title, displayed in header before close button.

- transition (a value equal to: "fade", "skew-up", "skew-down", "rotate-right", "rotate-left", "slide-down", "slide-up", "slide-right", "slide-left", "scale-y", "scale-x", "scale", "pop", "pop-top-left", "pop-top-right", "pop-bottom-left", "pop-bottom-right"; optional):
    Modal body transition.

- transitionDuration (number; optional):
    Duration in ms of modal transitions, set to 0 to disable all
    animations.

- transitionTimingFunction (string; optional):
    Modal body transitionTimingFunction, defaults to
    theme.transitionTimingFunction.

- trapFocus (boolean; optional):
    Disables focus trap.

- withCloseButton (boolean; optional):
    Hides close button if set to False, modal still can be closed with
    escape key and by clicking outside.

- zIndex (number; optional):
    Popper zIndex."""
    @_explicitize_args
    def __init__(self, children=None, centered=Component.UNDEFINED, class_name=Component.UNDEFINED, closeOnClickOutside=Component.UNDEFINED, closeOnEscape=Component.UNDEFINED, id=Component.UNDEFINED, opened=Component.UNDEFINED, overflow=Component.UNDEFINED, overlayColor=Component.UNDEFINED, overlayOpacity=Component.UNDEFINED, padding=Component.UNDEFINED, radius=Component.UNDEFINED, shadow=Component.UNDEFINED, size=Component.UNDEFINED, style=Component.UNDEFINED, title=Component.UNDEFINED, transition=Component.UNDEFINED, transitionDuration=Component.UNDEFINED, transitionTimingFunction=Component.UNDEFINED, trapFocus=Component.UNDEFINED, withCloseButton=Component.UNDEFINED, zIndex=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'id', 'centered', 'class_name', 'closeOnClickOutside', 'closeOnEscape', 'opened', 'overflow', 'overlayColor', 'overlayOpacity', 'padding', 'radius', 'shadow', 'size', 'style', 'title', 'transition', 'transitionDuration', 'transitionTimingFunction', 'trapFocus', 'withCloseButton', 'zIndex']
        self._type = 'Modal'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'id', 'centered', 'class_name', 'closeOnClickOutside', 'closeOnEscape', 'opened', 'overflow', 'overlayColor', 'overlayOpacity', 'padding', 'radius', 'shadow', 'size', 'style', 'title', 'transition', 'transitionDuration', 'transitionTimingFunction', 'trapFocus', 'withCloseButton', 'zIndex']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(Modal, self).__init__(children=children, **args)
