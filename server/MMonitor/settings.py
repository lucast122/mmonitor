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

STATICFILES_DIRS = [
    os.path.join(BASE_DIR, 'static'),
]

# Media files
MEDIA_URL = '/media/'
MEDIA_ROOT = BASE_DIR / 'media'

ALLOWED_HOSTS = ['127.0.0.1', 'localhost', 'mmonitor.org', 'www.mmonitor.org', '*']  # '*' for development only
CSRF_TRUSTED_ORIGINS = ['https://mmonitor.org', 'https://www.mmonitor.org']


INSTALLED_APPS = [
    "django.contrib.admin",
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.messages",
    "django.contrib.staticfiles",

    # Django REST Framework
    "rest_framework",
    "rest_framework_simplejwt",
    "corsheaders",
    "django_filters",

    # Project apps
    "api",
    "users",
    "main",
    "dashboard",

    "django_bootstrap5",
]



MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'corsheaders.middleware.CorsMiddleware',  # CORS - must be before CommonMiddleware
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]


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

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [
            BASE_DIR / 'templates',
            os.path.join(BASE_DIR, '..', 'frontend', 'dist'),
        ],
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
DATA_UPLOAD_MAX_MEMORY_SIZE = 100*1024*1024
RSCRIPT = 'RScript'


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


# =============================================================================
# Django REST Framework Configuration
# =============================================================================

REST_FRAMEWORK = {
    'DEFAULT_AUTHENTICATION_CLASSES': [
        'rest_framework_simplejwt.authentication.JWTAuthentication',
    ],
    'DEFAULT_PERMISSION_CLASSES': [
        'rest_framework.permissions.IsAuthenticated',
    ],
    'DEFAULT_PAGINATION_CLASS': 'api.pagination.StandardPagination',
    'PAGE_SIZE': 50,
    'DEFAULT_FILTER_BACKENDS': [
        'django_filters.rest_framework.DjangoFilterBackend',
        'rest_framework.filters.SearchFilter',
        'rest_framework.filters.OrderingFilter',
    ],
    'DEFAULT_THROTTLE_CLASSES': [
        'rest_framework.throttling.AnonRateThrottle',
        'rest_framework.throttling.UserRateThrottle'
    ],
    'DEFAULT_THROTTLE_RATES': {
        'anon': '100/hour',
        'user': '1000/hour'
    },
}

# =============================================================================
# Simple JWT Configuration
# =============================================================================

from datetime import timedelta

SIMPLE_JWT = {
    'ACCESS_TOKEN_LIFETIME': timedelta(minutes=60),
    'REFRESH_TOKEN_LIFETIME': timedelta(days=7),
    'ROTATE_REFRESH_TOKENS': True,
    'BLACKLIST_AFTER_ROTATION': True,
    'UPDATE_LAST_LOGIN': True,

    'ALGORITHM': 'HS256',
    'SIGNING_KEY': SECRET_KEY,
    'VERIFYING_KEY': None,
    'AUDIENCE': None,
    'ISSUER': None,

    'AUTH_HEADER_TYPES': ('Bearer',),
    'AUTH_HEADER_NAME': 'HTTP_AUTHORIZATION',
    'USER_ID_FIELD': 'id',
    'USER_ID_CLAIM': 'user_id',

    'AUTH_TOKEN_CLASSES': ('rest_framework_simplejwt.tokens.AccessToken',),
    'TOKEN_TYPE_CLAIM': 'token_type',
}

# =============================================================================
# CORS Configuration
# =============================================================================

# Allow all origins in development, restrict in production
CORS_ALLOWED_ORIGINS = [
    "http://localhost:3000",
    "http://localhost:5173",
    "http://127.0.0.1:3000",
    "http://127.0.0.1:5173",
    "https://mmonitor.org",
    "https://www.mmonitor.org",
]

CORS_ALLOW_CREDENTIALS = True

CORS_ALLOW_HEADERS = [
    'accept',
    'accept-encoding',
    'authorization',
    'content-type',
    'dnt',
    'origin',
    'user-agent',
    'x-csrftoken',
    'x-requested-with',
]

# In development, allow all origins
if DEBUG:
    CORS_ALLOW_ALL_ORIGINS = True
