# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class Menu(Component):
    """A Menu component.
Combine a list of secondary actions into single interactive area. For more information, see: https://mantine.dev/core/menu/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    dmc.MenuItem, dmc.MenuLabel, and dmc.Divider components only,
    children are passed to MenuBody component.

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- arrowDistance (number; optional):
    Arrow distance to the left/right * arrowSize.

- arrowSize (number; optional):
    Arrow size in px.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- closeOnItemClick (boolean; optional):
    Should menu close on item click.

- closeOnScroll (boolean; optional):
    Close menu on scroll.

- delay (number; optional):
    Close delay for hover trigger.

- exitTransitionDuration (number; optional):
    Unmount transition duration in ms.

- gutter (number; optional):
    Spacing between element and popper in px.

- placement (a value equal to: "center", "end", "start"; optional):
    Placement relative to reference element.

- position (a value equal to: "bottom", "left", "right", "top"; optional):
    Position relative to reference element.

- radius (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Menu body and items border-radius.

- shadow (a value equal to: "xs", "sm", "md", "lg", "xl" | string; optional):
    Predefined shadow from theme or box-shadow value.

- size (a value equal to: "auto", "xs", "sm", "md", "lg", "xl" | number; optional):
    Predefined menu width or number for width in px.

- style (dict; optional):
    Inline style override.

- transition (a value equal to: "fade", "skew-up", "skew-down", "rotate-right", "rotate-left", "slide-down", "slide-up", "slide-right", "slide-left", "scale-y", "scale-x", "scale", "pop", "pop-top-left", "pop-top-right", "pop-bottom-left", "pop-bottom-right"; optional):
    Customize mount/unmount transition.

- transitionDuration (number; optional):
    Mount transition duration in ms.

- transitionTimingFunction (string; optional):
    Mount/unmount transition timing function, defaults to
    theme.transitionTimingFunction.

- trapFocus (boolean; optional):
    Should focus be trapped when menu is opened.

- trigger (a value equal to: "hover", "click"; optional):
    Event which should open menu.

- withArrow (boolean; optional):
    Renders arrow if True.

- zIndex (number; optional):
    Menu body z-index."""
    @_explicitize_args
    def __init__(self, children=None, arrowDistance=Component.UNDEFINED, arrowSize=Component.UNDEFINED, class_name=Component.UNDEFINED, closeOnItemClick=Component.UNDEFINED, closeOnScroll=Component.UNDEFINED, delay=Component.UNDEFINED, exitTransitionDuration=Component.UNDEFINED, gutter=Component.UNDEFINED, id=Component.UNDEFINED, placement=Component.UNDEFINED, position=Component.UNDEFINED, radius=Component.UNDEFINED, shadow=Component.UNDEFINED, size=Component.UNDEFINED, style=Component.UNDEFINED, transition=Component.UNDEFINED, transitionDuration=Component.UNDEFINED, transitionTimingFunction=Component.UNDEFINED, trapFocus=Component.UNDEFINED, trigger=Component.UNDEFINED, withArrow=Component.UNDEFINED, zIndex=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'id', 'arrowDistance', 'arrowSize', 'class_name', 'closeOnItemClick', 'closeOnScroll', 'delay', 'exitTransitionDuration', 'gutter', 'placement', 'position', 'radius', 'shadow', 'size', 'style', 'transition', 'transitionDuration', 'transitionTimingFunction', 'trapFocus', 'trigger', 'withArrow', 'zIndex']
        self._type = 'Menu'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'id', 'arrowDistance', 'arrowSize', 'class_name', 'closeOnItemClick', 'closeOnScroll', 'delay', 'exitTransitionDuration', 'gutter', 'placement', 'position', 'radius', 'shadow', 'size', 'style', 'transition', 'transitionDuration', 'transitionTimingFunction', 'trapFocus', 'trigger', 'withArrow', 'zIndex']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(Menu, self).__init__(children=children, **args)
