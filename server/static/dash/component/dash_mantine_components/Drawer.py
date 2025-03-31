# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class Drawer(Component):
    """A Drawer component.
Display overlay area at any side of the screen. For more information, see: https://mantine.dev/core/drawer/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    Drawer children components.

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- closeOnClickOutside (boolean; optional):
    Disable onClick trigger for outside events.

- closeOnEscape (boolean; optional):
    Disable onClick trigger for escape key press.

- hideCloseButton (boolean; optional):
    Hides close button, modal still can be closed with escape key and
    by clicking outside.

- lockScroll (boolean; optional):
    Disables scroll lock.

- opened (boolean; default False):
    If True drawer is mounted to the dom.

- overlayColor (string; optional):
    Sets overlay color, defaults to theme.black in light theme and to
    theme.colors.dark[9] in dark theme.

- overlayOpacity (number; optional):
    Sets overlay opacity, defaults to 0.75 in light theme and to 0.85
    in dark theme.

- padding (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Drawer body padding from theme or number for padding in px.

- position (a value equal to: "left", "right", "top", "bottom"; optional):
    Drawer body position.

- shadow (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Drawer body shadow from theme or any css shadow value.

- size (a value equal to: "xs", "sm", "md", "lg", "xl" | string | number; optional):
    Drawer body width (right | left position) or height (top | bottom
    position), cannot exceed 100vh for height and 100% for width.

- title (boolean | number | string | dict | list; optional):
    Drawer title, displayed in header before close button.

- transition (a value equal to: "fade", "skew-up", "skew-down", "rotate-right", "rotate-left", "slide-down", "slide-up", "slide-right", "slide-left", "scale-y", "scale-x", "scale", "pop", "pop-top-left", "pop-top-right", "pop-bottom-left", "pop-bottom-right"; optional):
    Drawer appear and disappear transition, see Transition component
    for full documentation.

- transitionDuration (number; optional):
    Transition duration in ms.

- transitionTimingFunction (string; optional):
    Drawer transitionTimingFunction css property.

- trapFocus (boolean; optional):
    Disables focus trap.

- withCloseButton (boolean; optional):
    Hides close button if set to False, drawer still can be closed
    with escape key and by clicking outside.

- withOverlay (boolean; optional):
    Removes overlay entirely.

- zIndex (number; optional):
    Popper zIndex."""
    @_explicitize_args
    def __init__(self, children=None, class_name=Component.UNDEFINED, closeOnClickOutside=Component.UNDEFINED, closeOnEscape=Component.UNDEFINED, hideCloseButton=Component.UNDEFINED, id=Component.UNDEFINED, lockScroll=Component.UNDEFINED, opened=Component.UNDEFINED, overlayColor=Component.UNDEFINED, overlayOpacity=Component.UNDEFINED, padding=Component.UNDEFINED, position=Component.UNDEFINED, shadow=Component.UNDEFINED, size=Component.UNDEFINED, title=Component.UNDEFINED, transition=Component.UNDEFINED, transitionDuration=Component.UNDEFINED, transitionTimingFunction=Component.UNDEFINED, trapFocus=Component.UNDEFINED, withCloseButton=Component.UNDEFINED, withOverlay=Component.UNDEFINED, zIndex=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'id', 'class_name', 'closeOnClickOutside', 'closeOnEscape', 'hideCloseButton', 'lockScroll', 'opened', 'overlayColor', 'overlayOpacity', 'padding', 'position', 'shadow', 'size', 'title', 'transition', 'transitionDuration', 'transitionTimingFunction', 'trapFocus', 'withCloseButton', 'withOverlay', 'zIndex']
        self._type = 'Drawer'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'id', 'class_name', 'closeOnClickOutside', 'closeOnEscape', 'hideCloseButton', 'lockScroll', 'opened', 'overlayColor', 'overlayOpacity', 'padding', 'position', 'shadow', 'size', 'title', 'transition', 'transitionDuration', 'transitionTimingFunction', 'trapFocus', 'withCloseButton', 'withOverlay', 'zIndex']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(Drawer, self).__init__(children=children, **args)
