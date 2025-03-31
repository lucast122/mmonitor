# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class Accordion(Component):
    """An Accordion component.
Divide content into collapsible sections. For more information, see: https://mantine.dev/core/accordion/

Keyword arguments:

- children (a list of or a singular dash component, string or number; optional):
    <AccordionItem /> components only.

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- disableIconRotation (boolean; optional):
    Should icon rotation be disabled.

- icon (boolean | number | string | dict | list; optional):
    Replace icon on all items.

- iconPosition (a value equal to: "right", "left"; optional):
    Change icon position: left or right.

- iconSize (number; optional):
    Icon width in px.

- multiple (boolean; optional):
    Allow multiple items to be opened at the same time.

- offsetIcon (boolean; optional):
    Should icon be offset with padding, applicable only when
    iconPosition is left.

- order (a value equal to: 1, 2, 3, 4, 5, 6; optional):
    Heading level used for items.

- state (dict with strings as keys and values of type boolean; optional):
    Controlled state (controls opened state of accordion items).

- style (dict; optional):
    Inline style.

- transitionDuration (number; optional):
    Open/close item transition duration in ms."""
    @_explicitize_args
    def __init__(self, children=None, class_name=Component.UNDEFINED, disableIconRotation=Component.UNDEFINED, icon=Component.UNDEFINED, iconPosition=Component.UNDEFINED, iconSize=Component.UNDEFINED, id=Component.UNDEFINED, multiple=Component.UNDEFINED, offsetIcon=Component.UNDEFINED, order=Component.UNDEFINED, state=Component.UNDEFINED, style=Component.UNDEFINED, transitionDuration=Component.UNDEFINED, **kwargs):
        self._prop_names = ['children', 'id', 'class_name', 'disableIconRotation', 'icon', 'iconPosition', 'iconSize', 'multiple', 'offsetIcon', 'order', 'state', 'style', 'transitionDuration']
        self._type = 'Accordion'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['children', 'id', 'class_name', 'disableIconRotation', 'icon', 'iconPosition', 'iconSize', 'multiple', 'offsetIcon', 'order', 'state', 'style', 'transitionDuration']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(Accordion, self).__init__(children=children, **args)
