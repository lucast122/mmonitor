# AUTO GENERATED FILE - DO NOT EDIT

from dash.development.base_component import Component, _explicitize_args


class DatePicker(Component):
    """A DatePicker component.
Capture date input from user. For more information, see: https://mantine.dev/dates/date-picker/

Keyword arguments:

- id (string; optional):
    The ID of this component, used to identify dash components in
    callbacks.

- allowFreeInput (boolean; optional):
    Allow free input.

- allowLevelChange (boolean; optional):
    Allow to change level (date – month – year).

- amountOfMonths (number; optional):
    Amount of displayed months.

- class_name (string; optional):
    Often used with CSS to style elements with common properties.

- clearable (boolean; optional):
    Allow to clear value.

- closeCalendarOnChange (boolean; optional):
    Set to False to force dropdown to stay open after date was
    selected.

- closeDropdownOnScroll (boolean; optional):
    Set to True to disable dropdown closing on scroll.

- description (boolean | number | string | dict | list; optional):
    Input description, displayed after label.

- disableOutsideEvents (boolean; optional):
    When True dates that are outside of given month cannot be clicked
    or focused.

- disabled (boolean; optional):
    A Datepicker can show it is currently unable to be interacted
    with.

- dropdownType (a value equal to: "modal", "popover"; optional):
    Where to show calendar in modal or popover.

- error (boolean | number | string | dict | list; optional):
    Displays error message after input.

- firstDayOfWeek (a value equal to: "sunday", "monday"; optional):
    Set first day of the week.

- fixOnBlur (boolean; optional):
    call onChange with last valid value onBlur.

- focusable (boolean; optional):
    Should focusable days have tabIndex={0}?.

- fullWidth (boolean; optional):
    Set to True to make calendar take 100% of container width.

- hideOutsideDates (boolean; optional):
    Remove outside dates.

- hideWeekdays (boolean; optional):
    Set to False to remove weekdays row.

- icon (boolean | number | string | dict | list; optional):
    Adds icon on the left side of input.

- iconWidth (number; optional):
    Width of icon section in px.

- initialLevel (a value equal to: "date", "month", "year"; optional):
    Initial date selection level.

- initialMonth (string; optional):
    Initial selected month.

- initiallyOpened (boolean; optional):
    Control initial dropdown opened state.

- inputFormat (string; optional):
    dayjs input format.

- invalid (boolean; optional):
    Sets border color to red and aria-invalid=True on input element.

- label (boolean | number | string | dict | list; optional):
    Input label, displayed before input.

- locale (string; optional):
    Locale used for all labels formatting.

- maxDate (string; optional):
    Maximum possible date.

- minDate (string; optional):
    Minimum possible date.

- modalZIndex (number; optional):
    Modal z-index.

- multiline (boolean; optional):
    Will input have multiple lines?.

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

- placeholder (string; optional):
    Placeholder, displayed when date is not selected.

- preventFocus (boolean; optional):
    Prevent focusing upon clicking.

- radius (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Input border-radius from theme or number to set border-radius in
    px.

- required (boolean; optional):
    Adds red asterisk on the right side of label.

- rightSection (boolean | number | string | dict | list; optional):
    Right section of input, similar to icon but on the right.

- rightSectionWidth (number; optional):
    Width of right section, is used to calculate input padding-right.

- shadow (a value equal to: "xs", "sm", "md", "lg", "xl" | number; optional):
    Dropdown shadow from theme or css value for custom box-shadow.

- size (a value equal to: "xs", "sm", "md", "lg", "xl"; optional):
    Input size.

- style (dict; optional):
    Inline style override.

- transition (a value equal to: "fade", "skew-up", "skew-down", "rotate-right", "rotate-left", "slide-down", "slide-up", "slide-right", "slide-left", "scale-y", "scale-x", "scale", "pop", "pop-top-left", "pop-top-right", "pop-bottom-left", "pop-bottom-right"; optional):
    Dropdown appear/disappear transition.

- transitionDuration (number; optional):
    Dropdown appear/disappear transition duration.

- value (string; optional):
    Selected date.

- variant (a value equal to: "default", "filled", "unstyled", "headless"; optional):
    Defines input appearance, defaults to default in light color
    scheme and filled in dark.

- withinPortal (boolean; optional):
    Whether to render the dropdown in a Portal.

- zIndex (number; optional):
    Popper zIndex."""
    @_explicitize_args
    def __init__(self, allowFreeInput=Component.UNDEFINED, allowLevelChange=Component.UNDEFINED, amountOfMonths=Component.UNDEFINED, class_name=Component.UNDEFINED, clearable=Component.UNDEFINED, closeCalendarOnChange=Component.UNDEFINED, closeDropdownOnScroll=Component.UNDEFINED, description=Component.UNDEFINED, disableOutsideEvents=Component.UNDEFINED, disabled=Component.UNDEFINED, dropdownType=Component.UNDEFINED, error=Component.UNDEFINED, firstDayOfWeek=Component.UNDEFINED, fixOnBlur=Component.UNDEFINED, focusable=Component.UNDEFINED, fullWidth=Component.UNDEFINED, hideOutsideDates=Component.UNDEFINED, hideWeekdays=Component.UNDEFINED, icon=Component.UNDEFINED, iconWidth=Component.UNDEFINED, id=Component.UNDEFINED, initialLevel=Component.UNDEFINED, initialMonth=Component.UNDEFINED, initiallyOpened=Component.UNDEFINED, inputFormat=Component.UNDEFINED, invalid=Component.UNDEFINED, label=Component.UNDEFINED, locale=Component.UNDEFINED, maxDate=Component.UNDEFINED, minDate=Component.UNDEFINED, modalZIndex=Component.UNDEFINED, multiline=Component.UNDEFINED, persistence=Component.UNDEFINED, persisted_props=Component.UNDEFINED, persistence_type=Component.UNDEFINED, placeholder=Component.UNDEFINED, preventFocus=Component.UNDEFINED, radius=Component.UNDEFINED, required=Component.UNDEFINED, rightSection=Component.UNDEFINED, rightSectionWidth=Component.UNDEFINED, shadow=Component.UNDEFINED, size=Component.UNDEFINED, style=Component.UNDEFINED, transition=Component.UNDEFINED, transitionDuration=Component.UNDEFINED, value=Component.UNDEFINED, variant=Component.UNDEFINED, withinPortal=Component.UNDEFINED, zIndex=Component.UNDEFINED, **kwargs):
        self._prop_names = ['id', 'allowFreeInput', 'allowLevelChange', 'amountOfMonths', 'class_name', 'clearable', 'closeCalendarOnChange', 'closeDropdownOnScroll', 'description', 'disableOutsideEvents', 'disabled', 'dropdownType', 'error', 'firstDayOfWeek', 'fixOnBlur', 'focusable', 'fullWidth', 'hideOutsideDates', 'hideWeekdays', 'icon', 'iconWidth', 'initialLevel', 'initialMonth', 'initiallyOpened', 'inputFormat', 'invalid', 'label', 'locale', 'maxDate', 'minDate', 'modalZIndex', 'multiline', 'persisted_props', 'persistence', 'persistence_type', 'placeholder', 'preventFocus', 'radius', 'required', 'rightSection', 'rightSectionWidth', 'shadow', 'size', 'style', 'transition', 'transitionDuration', 'value', 'variant', 'withinPortal', 'zIndex']
        self._type = 'DatePicker'
        self._namespace = 'dash_mantine_components'
        self._valid_wildcard_attributes =            []
        self.available_properties = ['id', 'allowFreeInput', 'allowLevelChange', 'amountOfMonths', 'class_name', 'clearable', 'closeCalendarOnChange', 'closeDropdownOnScroll', 'description', 'disableOutsideEvents', 'disabled', 'dropdownType', 'error', 'firstDayOfWeek', 'fixOnBlur', 'focusable', 'fullWidth', 'hideOutsideDates', 'hideWeekdays', 'icon', 'iconWidth', 'initialLevel', 'initialMonth', 'initiallyOpened', 'inputFormat', 'invalid', 'label', 'locale', 'maxDate', 'minDate', 'modalZIndex', 'multiline', 'persisted_props', 'persistence', 'persistence_type', 'placeholder', 'preventFocus', 'radius', 'required', 'rightSection', 'rightSectionWidth', 'shadow', 'size', 'style', 'transition', 'transitionDuration', 'value', 'variant', 'withinPortal', 'zIndex']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(DatePicker, self).__init__(**args)
