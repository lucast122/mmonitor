from django.urls import path, include

from . import views
from .views import UploadMAGAnnotationsView

app_name = "users"
urlpatterns = [
    # Custom logout view that supports both GET and POST
    path("logout/", views.custom_logout, name="logout"),
    
    # Include default auth_urls (our custom logout will override the auth one)
    path("", include("django.contrib.auth.urls")),
    
    # Registration page
    path("register/", views.register, name="register"),
    path('get_user_id/', views.get_user_id, name='get_user_id'),
    path('add_nanopore_record/', views.add_nanopore_record, name='add_nanopore_record'),
    path('overwrite_nanopore_record/', views.overwrite_nanopore_record, name='overwrite_nanopore_record'),
    path('add_sequencing_statistics/', views.add_sequencing_statistics, name='add_sequencing_statistics'),
    path('get_unique_sample_ids/', views.get_unique_sample_ids, name='get_unique_sample_ids'),
    path('feedback/submit/', views.submit_feedback, name='submit_feedback'),
    path('feedback/thank-you/', views.feedback_thank_you, name='feedback_thank_you'),
    path('upload_mag/', views.upload_mag, name='upload_mag'),
    path('mag-viewer/', views.mag_viewer, name='mag_viewer'),
    path('upload_annotations/', UploadMAGAnnotationsView.as_view(), name='upload_annotations'),
]