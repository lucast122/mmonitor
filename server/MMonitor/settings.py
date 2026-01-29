# settings.py

import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent

# Ensure DEBUG is set to False in production
DEBUG = True  # Set to True for development

# Static files
STATIC_URL = '/static/'
STATIC_ROOT = os.path.join(BASE_DIR, 'staticfiles')
SITE_URL = 'https://mmonitor.org/'

# Find dash component directories
import dash
import dash_mantine_components
import dash_bootstrap_components
import dash_core_components
import dash_html_components
import dash_table
import dash_bio 
import dash_ag_grid
import dash_iconify

# Add paths to Dash components
venv_path = os.path.join(BASE_DIR, 'py3.11venv', 'lib', 'python3.11', 'site-packages')
STATICFILES_DIRS = [
    os.path.join(BASE_DIR, 'static'),
    # Add Dash component paths with their correct package location
    os.path.dirname(dash.__file__),
    os.path.dirname(dash_mantine_components.__file__),
    os.path.dirname(dash_bootstrap_components.__file__),
    os.path.dirname(dash_core_components.__file__),
    os.path.dirname(dash_html_components.__file__),
    os.path.dirname(dash_table.__file__),
    os.path.dirname(dash_bio.__file__), 
    os.path.dirname(dash_ag_grid.__file__),
    os.path.dirname(dash_iconify.__file__),
]

# Media files
MEDIA_URL = '/media/'
MEDIA_ROOT = BASE_DIR / 'media'

ALLOWED_HOSTS = ['127.0.0.1', 'localhost', 'mmonitor.org', 'www.mmonitor.org']
CSRF_TRUSTED_ORIGINS = ['https://mmonitor.org', 'https://www.mmonitor.org']


INSTALLED_APPS = [
    "django.contrib.admin",
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.messages",
    "django.contrib.staticfiles",

    "django_plotly_dash.apps.DjangoPlotlyDashConfig",
    "dpd_static_support",

    "users",
    "main",
    "dashboard",
    "dashboard.dashapp",

    "channels",
    "django_bootstrap5",
]



MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    'django_plotly_dash.middleware.ExternalRedirectionMiddleware',
    'django_plotly_dash.middleware.BaseMiddleware',
]

PLOTLY_COMPONENTS = [
    "dash_core_components",
    "dash_html_components",
    "dash_renderer",
    "dash_table",
    "dash_bootstrap_components",
    "dash_mantine_components",
    "dash_ag_grid",
    "dash_iconify",
]


PLOTLY_DASH = {
    "serve_locally": False,
    "ws_route": "dpd/ws/channel",
    "http_route": "dpd/views",
    "http_poke_enabled": True,
    "cache_timeout_initial_arguments": 60,
}




# Static files finders
STATICFILES_FINDERS = [
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
]

# Use simple storage backend
STATICFILES_STORAGE = 'django.contrib.staticfiles.storage.StaticFilesStorage'

# Whitenoise settings
WHITENOISE_USE_FINDERS = True
WHITENOISE_ROOT = STATIC_ROOT
WHITENOISE_MAX_AGE = 31536000  # 1 year in seconds
WHITENOISE_SKIP_COMPRESS_EXTENSIONS = []  # Empty list to compress all files

ROOT_URLCONF = 'MMonitor.urls'
DASH_APP_URL = 'https://www.mmonitor.org/dash/'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [BASE_DIR / 'templates'],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

ASGI_APPLICATION = 'MMonitor.asgi.application'

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': BASE_DIR / 'db.sqlite3',
    }
    
}

# Configure the 'default' database dynamically from the DATABASE_URL environment variable
if 'DATABASE_URL' in os.environ:
    import dj_database_url
    DATABASES['default'] = dj_database_url.config(conn_max_age=600, ssl_require=False)

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]

SECRET_KEY = 'django-insecure-9v8c1b(!um4d75#+1gmrg$qunrmtmxk7g7&u)5od-^8vhpk5&^'
LANGUAGE_CODE = 'en-us'
TIME_ZONE = 'UTC'
USE_I18N = True
USE_TZ = True

DEFAULT_AUTO_FIELD = 'django.db.models.BigAutoField'

LOGIN_URL = 'users:login'
LOGOUT_REDIRECT_URL = "main:index"
LOGIN_REDIRECT_URL = '/'

X_FRAME_OPTIONS = 'SAMEORIGIN'
DATA_UPLOAD_MAX_MEMORY_SIZE = 429916160
RSCRIPT = 'RScript'
DATA_UPLOAD_MAX_MEMORY_SIZE = 100*1024*1024


import os

LOG_FILE_PATH = os.path.join(os.path.dirname(__file__), 'django_log.txt')

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'handlers': {
        'file': {
            'level': 'DEBUG',
            'class': 'logging.FileHandler',
            'filename': LOG_FILE_PATH,
        },
    },
    'root': {
        'handlers': ['file'],
        'level': 'DEBUG',
    },
}


# DATABASE_ROUTERS = ['MMonitor.routers.AppDataRouter']  # Commented out - file doesn't exist

# Security settings
X_FRAME_OPTIONS = 'SAMEORIGIN'  # Allow iframe embedding from same origin
SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTO', 'https')
CSRF_COOKIE_SECURE = True
SESSION_COOKIE_SECURE = True
CSRF_TRUSTED_ORIGINS = ['https://mmonitor.org', 'https://www.mmonitor.org']
DPD_DEFAULT_DOWNLOAD = True
