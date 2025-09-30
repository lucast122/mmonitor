#!/usr/bin/env python
"""
Script to create database tables directly without migrations
"""
import os
import sys
import django

# Add the server directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Setup Django
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'MMonitor.settings')
django.setup()

from django.core.management import execute_from_command_line
from django.db import connection
from users.models import UserQuota, UploadedFile, AnalysisJob, SystemStatus
from django.contrib.auth.models import User

def create_tables():
    """Create tables directly"""
    with connection.schema_editor() as schema_editor:
        # Create UserQuota table
        try:
            schema_editor.create_model(UserQuota)
            print("Created UserQuota table")
        except Exception as e:
            print(f"UserQuota table creation failed: {e}")
        
        # Create UploadedFile table  
        try:
            schema_editor.create_model(UploadedFile)
            print("Created UploadedFile table")
        except Exception as e:
            print(f"UploadedFile table creation failed: {e}")
            
        # Create AnalysisJob table
        try:
            schema_editor.create_model(AnalysisJob)
            print("Created AnalysisJob table")
        except Exception as e:
            print(f"AnalysisJob table creation failed: {e}")
            
        # Create SystemStatus table
        try:
            schema_editor.create_model(SystemStatus)
            print("Created SystemStatus table")
        except Exception as e:
            print(f"SystemStatus table creation failed: {e}")

def create_superuser():
    """Create a superuser for testing"""
    try:
        if not User.objects.filter(username='admin').exists():
            User.objects.create_superuser('admin', 'admin@example.com', 'admin123')
            print("Created superuser: admin/admin123")
        else:
            print("Superuser 'admin' already exists")
    except Exception as e:
        print(f"Superuser creation failed: {e}")

def initialize_system_status():
    """Initialize system status"""
    try:
        SystemStatus.get_current()
        print("Initialized SystemStatus")
    except Exception as e:
        print(f"SystemStatus initialization failed: {e}")

if __name__ == '__main__':
    print("Creating database tables...")
    create_tables()
    print("\nCreating superuser...")
    create_superuser()
    print("\nInitializing system status...")
    initialize_system_status()
    print("\nDone!")
