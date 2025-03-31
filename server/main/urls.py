from django.urls import path

from . import views

app_name = "main"
urlpatterns = [
    # Homepage
    path("", views.index, name="index"),


]
