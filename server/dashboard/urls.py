from django.urls import path, re_path
from . import views

app_name = 'dashboard'

urlpatterns = [
    path('get_user_id/', views.get_user_id, name='get_user_id'),
    path('', views.spa_view, name='index'),
    re_path(r'^(?P<path>.+)$', views.spa_view, name='spa_catchall'),
]
