from abc import ABC, abstractmethod
import pandas as pd
from json import loads, dumps

class BaseDashApp(ABC):

    def __init__(self, common_property):
        self.common_property = common_property


    # Common method shared across all child classes
    def common_method(self):
        print("This is a common method used by all Dash apps.")

    def _explode_metadata(df):
        return pd.concat([df, df['data'].apply(_parse_dict)], axis=1).drop(columns='data')

    def _parse_dict(x):
        return pd.Series(loads(x))

    def _init_mysql(self):
        conn_settings = settings.DATABASES['mmonitor']
        dialect = 'mysql'
        user = conn_settings['USER']
        password = conn_settings['PASSWORD']
        host = conn_settings['HOST']
        port = conn_settings['PORT']
        db_name = conn_settings['NAME']
        db_url = f'{dialect}://{user}:{password}@{host}:{port}/{db_name}'
        self._engine = create_engine(db_url)


    # Abstract method that must be implemented by all child classes
    @abstractmethod
    def specific_method(self):
        pass

    # You can define more abstract or common methods here
