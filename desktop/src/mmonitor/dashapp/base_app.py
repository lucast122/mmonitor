from dash import html

from mmonitor.database.mmonitor_db import MMonitorDBInterface


class BaseApp:
    """
    Base class for applications that can be registered in
    the main dash app. Requires the initialization of the
    database instance, the page content (layout), and the
    page's callbacks.
    """

    def __init__(self, sql: MMonitorDBInterface):
        self._sql = sql
        self._layout = None

    def _init_layout(self) -> None:
        pass

    def _init_callbacks(self) -> None:
        pass

    @property
    def layout(self) -> html:
        if self._layout is None:
            raise RuntimeError('Layout has not been initialized yet!')
        else:
            return self._layout

    @layout.setter
    def layout(self, value):
        self._layout = value
