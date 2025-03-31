from django.urls import path, include
from . import views
from .dashapp.index import Index
from .dashapp.taxonomy import Taxonomy
from .dashapp.diversity import Diversity
from .dashapp.correlations import Correlations
from .dashapp.horizon import Horizon
from .dashapp.qc import QC
from .dashapp.mags import MAGs

app_name = 'dashboard'

urlpatterns = [
    path('', views.load_app, {'name': 'Index'}, name='index'),
        path('get_user_id/', views.get_user_id, name='get_user_id')

]
