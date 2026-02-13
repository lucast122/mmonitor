"""MMonitor URL Configuration"""
from django.contrib import admin
from django.urls import path, include, re_path
from django.conf import settings
from django.conf.urls.static import static
from django.views.static import serve
import os

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', include('main.urls')),
    path('users/', include('users.urls')),

    # Serve React frontend assets before the dashboard catch-all
    re_path(r'^dashboard/assets/(?P<path>.*)$', serve, {
        'document_root': os.path.join(settings.BASE_DIR, '..', 'frontend', 'dist', 'assets')
    }),

    path('dashboard/', include('dashboard.urls')),

    # REST API v1
    path('api/v1/', include('api.urls', namespace='api')),

    # Serve static files in development
    re_path(r'^static/(?P<path>.*)$', serve, {'document_root': settings.STATIC_ROOT}),
]

# Serve static files during development and production
if settings.DEBUG:
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
else:
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
