# MMonitor Web Upload System

## Overview

The MMonitor Web Upload System is a comprehensive web-based interface that allows users to upload FASTQ files and run metagenomics analyses on your server. This system includes user quotas, job queuing, progress tracking, and results management.

## Features Implemented

### 🗄️ Database Models
- **UserQuota**: Daily limits (3 jobs, 2GB storage per user)
- **UploadedFile**: FASTQ file storage with validation and deduplication
- **AnalysisJob**: Job tracking with status, progress, and results
- **SystemStatus**: System-wide monitoring and limits

### 🌐 Web Interface
- **Upload Dashboard**: User quota display and recent jobs overview
- **File Upload**: Drag-and-drop interface with validation
- **Job Management**: List, filter, and monitor analysis jobs
- **Job Status**: Real-time progress tracking with auto-refresh
- **Results Download**: ZIP file download of completed analyses
- **Demo Page**: Showcase system capabilities and statistics

### 🔒 Security & Limits
- Daily quotas: 3 jobs and 2GB storage per user
- File validation: FASTQ format checking
- File deduplication using SHA256 hashes
- System capacity management (5 concurrent jobs, 50 queue slots)
- User authentication required for uploads

### ⚙️ Job Processing
- Background job processor with daemon mode
- Support for two analysis types:
  - Taxonomic Classification (Centrifuge) - 5-10 minutes
  - Full MMonitor Pipeline - 30-60 minutes
- Progress tracking and error handling
- Automatic retry mechanism (up to 3 attempts)
- Job cancellation for queued jobs

## File Structure

```
server/
├── users/
│   ├── models.py              # Database models (updated)
│   ├── forms.py               # Upload and job forms (updated)
│   ├── views.py               # Original views (cleaned up)
│   ├── upload_views.py        # New upload system views
│   ├── urls.py                # URL patterns (updated)
│   ├── signals.py             # Auto-create user quotas (updated)
│   ├── migrations/
│   │   └── 0005_upload_system_models.py  # Database migration
│   ├── management/commands/
│   │   └── process_jobs.py    # Job processor command
│   └── templates/users/
│       ├── upload_dashboard.html
│       ├── upload_files.html
│       ├── job_list.html
│       ├── job_status.html
│       └── demo_page.html
└── templates/
    └── base.html              # Updated navigation menu
```

## Navigation Menu

The system adds the following menu items:
- **Home**: Main landing page
- **Demo**: System showcase and statistics (public)
- **Upload**: Upload dashboard (authenticated users only)
- **Dashboard**: Existing dashboard (authenticated users only)

## Setup Instructions

### 1. Database Migration

Run the migration to create the new database tables:

```bash
python manage.py migrate users
```

### 2. Create System Status

Initialize the system status record:

```python
python manage.py shell
>>> from users.models import SystemStatus
>>> SystemStatus.get_current()
>>> exit()
```

### 3. Start Job Processor

Start the background job processor:

```bash
# One-time processing
python manage.py process_jobs

# Daemon mode (recommended for production)
python manage.py process_jobs --daemon --interval 10
```

### 4. Configure Media Files

Ensure Django is configured to serve uploaded files:

```python
# In settings.py
MEDIA_URL = '/media/'
MEDIA_ROOT = BASE_DIR / 'media'

# In urls.py (for development)
if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
```

### 5. File Storage

Create directories for file storage:

```bash
mkdir -p media/uploads
mkdir -p media/analysis_jobs
```

## Usage Workflow

1. **User Registration/Login**: Users must be authenticated to upload files
2. **Upload Dashboard**: Users see their quota status and recent jobs
3. **File Upload**: Users upload FASTQ files and configure analysis
4. **Job Queuing**: System validates files and queues jobs
5. **Processing**: Background processor runs analyses
6. **Results**: Users can monitor progress and download results

## API Endpoints

The system provides the following user-facing URLs:

- `/users/upload/` - Upload dashboard
- `/users/upload/files/` - File upload form
- `/users/jobs/` - Job list with filtering
- `/users/jobs/<job_id>/` - Job status and details
- `/users/jobs/<job_id>/download/` - Download results
- `/users/demo/` - Demo and statistics page

## System Administration

### Monitoring Jobs

```python
from users.models import AnalysisJob, SystemStatus

# Check system status
status = SystemStatus.get_current()
print(f"Active jobs: {status.active_jobs}")
print(f"Queued jobs: {status.queued_jobs}")

# View recent jobs
recent_jobs = AnalysisJob.objects.filter(status='running')[:10]
```

### Managing Quotas

```python
from users.models import UserQuota

# Reset user quota
user_quota = UserQuota.objects.get(user__username='username')
user_quota.jobs_used_today = 0
user_quota.storage_used_today = 0
user_quota.save()
```

### Maintenance Mode

```python
from users.models import SystemStatus

# Enable maintenance mode
status = SystemStatus.get_current()
status.maintenance_mode = True
status.maintenance_message = "System maintenance in progress"
status.save()
```

## Integration with Existing Pipeline

The job processor includes placeholder functions that need to be integrated with your existing MMonitor pipeline:

- `run_centrifuge()`: Integrate with your Centrifuge setup
- `run_full_pipeline_analysis()`: Connect to your full pipeline
- `run_quality_control()`: Add your QC tools
- `run_assembly()`: Connect to your assembly tools
- `run_annotation()`: Integrate annotation pipeline

## Security Considerations

1. **File Validation**: All uploaded files are validated for FASTQ format
2. **File Size Limits**: 500MB per file, 2GB daily limit per user
3. **User Quotas**: Prevent abuse with daily job and storage limits
4. **System Limits**: Prevent server overload with concurrent job limits
5. **File Deduplication**: Prevent storage waste with hash-based deduplication

## Customization Options

### Adjusting Limits

Modify limits in `users/models.py`:

```python
class UserQuota(models.Model):
    DAILY_JOB_LIMIT = 5  # Change from 3 to 5
    DAILY_STORAGE_LIMIT = 5 * 1024 * 1024 * 1024  # Change to 5GB
```

### Adding Analysis Types

Add new analysis types in `users/models.py`:

```python
ANALYSIS_TYPES = [
    ('centrifuge', 'Taxonomic Classification (Centrifuge)'),
    ('full_pipeline', 'Full MMonitor Pipeline'),
    ('custom_analysis', 'Custom Analysis'),  # Add new type
]
```

## Troubleshooting

### Common Issues

1. **Migration Errors**: Ensure all dependencies are installed
2. **File Upload Errors**: Check media directory permissions
3. **Job Processing Errors**: Verify job processor is running
4. **Template Errors**: Ensure all template files are in place

### Logs

Check Django logs for detailed error information:

```bash
tail -f django_log.txt
```

## Next Steps

1. **Integration**: Connect the job processor with your actual MMonitor pipeline
2. **Testing**: Test with real FASTQ files and validate results
3. **Production**: Deploy with proper web server (nginx/Apache) and database
4. **Monitoring**: Set up system monitoring and alerting
5. **Documentation**: Create user guides and tutorials

## Support

The upload system is designed to be modular and extensible. You can easily:

- Add new analysis types
- Modify quotas and limits
- Customize the user interface
- Integrate with external tools
- Add new file formats

For questions or issues, refer to the code comments and Django documentation.
