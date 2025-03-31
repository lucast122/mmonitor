class AppDataRouter:
    """
    A router to control all database operations on models in themmonitor application.
    """
    def db_for_read(self, model, **hints):
        """
        Attempts to read mmonitor models go to mmonitor database.
        """
        if model._meta.app_label == 'mmonitor':
            return 'mmonitor'
        return None

    def db_for_write(self, model, **hints):
        """
        Attempts to write mmonitor models go to mmonitor database.
        """
        if model._meta.app_label == 'mmonitor':
            return 'mmonitor'
        return None

    def allow_relation(self, obj1, obj2, **hints):
        """
        Allow relations if a model in the mmonitor app is involved.
        """
        if obj1._meta.app_label == 'mmonitor' or \
           obj2._meta.app_label == 'mmonitor':
           return True
        return None

    def allow_migrate(self, db, app_label, model_name=None, **hints):
        """
        Make sure the mmonitor app only appears in the 'mmonitor'
        database.
        """
        if app_label == 'mmonitor':
            return db == 'mmonitor'
        return None
