# MMonitor Developer


## Todo

### Desktop

- Release:
	- include all binaries in build script
	- build release on all desired platforms

- Server integration:
	- expand with user / database login form & offline mode
	- integrate MySQL (consider SQAlchemy ORM)

### Server

- Release:
	- Migrate to MySQL (consider Django ORM)

- Desktop integration:
	- Create user in MySQL instance when user registers on the server 


## Add new Dash App

### Desktop

- create new class in `mmonitor.dashapp.app` inheriting from `mmonitor.dashapp.base_app.BaseApp`
- implement layout and callbacks 
- register app in `mmonitor.dashapp.index._init_apps()`:
	- add to `self._apps`: `'/apps/<app>': {'name': '<app>', 'app': <app>.<App>(self._sql)}`

### Server

- copy & paste desktop app file to `dashboard.dashapp`
- adjust code
	- remove `BaseApp` inheritance
	- add `from django_plotly_dash import DjangoDash` to imports
	- add `self.app = DjangoDash('<app>')` to apps' `__init__()`
- register app in `dashboard.urls`:
	- add app to imports `from .dashapp import <app>`
	- create app object `<app>.<App>(db)`
	- add to urlpatterns `path('<app>/', views.load_app, {'name': '<app>'}, name='<app>'),`


## Future Considerations

- Switch from raw queries to an ORM for easier migrations & more robust structure
- Find a way to implement the horizon plot generation & serving natively
- or include an Rscript runtime + necessary libraries in the project
