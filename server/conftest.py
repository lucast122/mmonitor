import os
import django
from django.conf import settings

# Set the Django settings module
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'MMonitor.settings')  # Replace 'MMonitor' with your project name if different

# Configure Django settings
django.setup()
