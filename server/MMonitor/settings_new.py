# settings.py

import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent

# Static files
STATIC_URL = '/static/'
STATIC_ROOT = os.path.join(BASE_DIR, 'staticfiles')

STATICFILES_DIRS = [
    os.path.join(BASE_DIR, 'static'),
]


# Media files
MEDIA_URL = '/media/'
MEDIA_ROOT = BASE_DIR / 'media'

# Ensure DEBUG is set to False in production
DEBUG = True

ALLOWED_HOSTS = ['127.0.0.1', 'localhost', 'mmonitor.de', 'www.mmonitor.de', '134.2.78.150', 'mmonitor.org', 'www.mmonitor.org']
CSRF_TRUSTED_ORIGINS = ['https://mmonitor.org', 'https://www.mmonitor.org']

# Application definition
INSTALLED_APPS = [
        'django_plotly_dash',
    'plotly',
    'users',
    'main',
    'dashboard',
    'dashboard.dashapp',
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django_plotly_dash.apps.DjangoPlotlyDashConfig',
    'django_bootstrap5',
    'channels',
    'dash',
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

PLOTLY_DASH = {
    "ws_route": "dpd/ws/channel",
    "http_route": "dpd/views",
    "http_poke_enabled": True,
    "insert_demo_migrations": False,
    "cache_timeout_initial_arguments": 60,
    "view_decorator": None,
    "cache_arguments": False,
    "serve_locally": False,
}

PLOTLY_COMPONENTS = [
    'dash_core_components',
    'dash_renderer',
    'dpd_components',
]

STATICFILES_FINDERS = [
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
    'django_plotly_dash.finders.DashAssetFinder',
    'django_plotly_dash.finders.DashComponentFinder',
    'django_plotly_dash.finders.DashAppDirectoryFinder',
]

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
    },
    'mmonitor': {
        'ENGINE': 'django.db.backends.mysql',
        'NAME': 'mydjangodb',
        'USER': 'mmonitor_remote',
        'PASSWORD': 'asdf!minion',
        'HOST': '134.2.78.150',   # Or an IP Address that your DB is hosted on
        'PORT': '3306',
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

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'handlers': {
        'file': {
            'level': 'DEBUG',
            'class': 'logging.FileHandler',
            'filename': '/Users/timo/Downloads/home/minion-computer/mmonitor_production/MMonitor/server/django_log.txt',
        },
    },
    'root': {
        'handlers': ['file'],
        'level': 'DEBUG',
    },
}

DATABASE_ROUTERS = ['MMonitor.routers.AppDataRouter']
