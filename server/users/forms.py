from django import forms
from django.core.exceptions import ValidationError
from .models import Feedback, UploadedFile, AnalysisJob
import os


class FeedbackForm(forms.ModelForm):
    class Meta:
        model = Feedback
        fields = ['subject', 'message']


class MultipleFileInput(forms.ClearableFileInput):
    allow_multiple_selected = True


class MultipleFileField(forms.FileField):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("widget", MultipleFileInput())
        super().__init__(*args, **kwargs)

    def clean(self, data, initial=None):
        single_file_clean = super().clean
        if isinstance(data, (list, tuple)):
            result = [single_file_clean(d, initial) for d in data]
        else:
            result = single_file_clean(data, initial)
        return result


class FileUploadForm(forms.Form):
    """Form for uploading FASTQ files"""
    files = MultipleFileField(
        widget=MultipleFileInput(attrs={
            'accept': '.fastq,.fq,.fastq.gz,.fq.gz',
            'class': 'form-control'
        }),
        help_text="Select one or more FASTQ files (.fastq, .fq, .fastq.gz, .fq.gz). Max 500MB per file."
    )
    
    def clean_files(self):
        files = self.files.getlist('files')
        
        if not files:
            raise ValidationError("Please select at least one file.")
        
        total_size = 0
        for uploaded_file in files:
            # Check file extension
            filename = uploaded_file.name.lower()
            if not any(filename.endswith(ext) for ext in UploadedFile.ALLOWED_EXTENSIONS):
                raise ValidationError(
                    f"File '{uploaded_file.name}' must be a FASTQ file with extension: "
                    f"{', '.join(UploadedFile.ALLOWED_EXTENSIONS)}"
                )
            
            # Check file size
            if uploaded_file.size > UploadedFile.MAX_FILE_SIZE:
                max_mb = UploadedFile.MAX_FILE_SIZE // 1024 // 1024
                raise ValidationError(f"File '{uploaded_file.name}' is too large. Max size: {max_mb}MB")
            
            total_size += uploaded_file.size
        
        # Check total size (2GB limit)
        max_total = 2 * 1024 * 1024 * 1024  # 2GB
        if total_size > max_total:
            raise ValidationError(f"Total file size exceeds 2GB limit.")
        
        return files


class AnalysisJobForm(forms.ModelForm):
    """Form for creating analysis jobs"""
    
    class Meta:
        model = AnalysisJob
        fields = ['sample_name', 'analysis_type']
        widgets = {
            'sample_name': forms.TextInput(attrs={
                'class': 'form-control',
                'placeholder': 'Enter a name for your sample (e.g., Sample_001)',
                'maxlength': 255
            }),
            'analysis_type': forms.Select(attrs={
                'class': 'form-select'
            })
        }
    
    def __init__(self, *args, **kwargs):
        self.user = kwargs.pop('user', None)
        super().__init__(*args, **kwargs)
        
        # Add help text
        self.fields['sample_name'].help_text = "A unique identifier for your sample"
        self.fields['analysis_type'].help_text = "Choose the type of analysis to perform"
    
    def clean_sample_name(self):
        sample_name = self.cleaned_data.get('sample_name', '').strip()
        
        if not sample_name:
            raise ValidationError("Sample name is required.")
        
        # Check for valid characters (alphanumeric, underscore, hyphen)
        if not all(c.isalnum() or c in '_-' for c in sample_name.replace(' ', '_')):
            raise ValidationError("Sample name can only contain letters, numbers, spaces, underscores, and hyphens.")
        
        # Replace spaces with underscores
        sample_name = sample_name.replace(' ', '_')
        
        # Check if user already has a job with this sample name
        if self.user and AnalysisJob.objects.filter(user=self.user, sample_name=sample_name).exists():
            raise ValidationError("You already have a job with this sample name. Please choose a different name.")
        
        return sample_name


class JobFilterForm(forms.Form):
    """Form for filtering analysis jobs"""
    
    STATUS_CHOICES = [('', 'All Statuses')] + AnalysisJob.STATUS_CHOICES
    ANALYSIS_TYPE_CHOICES = [('', 'All Types')] + AnalysisJob.ANALYSIS_TYPES
    
    status = forms.ChoiceField(
        choices=STATUS_CHOICES,
        required=False,
        widget=forms.Select(attrs={'class': 'form-select form-select-sm'})
    )
    
    analysis_type = forms.ChoiceField(
        choices=ANALYSIS_TYPE_CHOICES,
        required=False,
        widget=forms.Select(attrs={'class': 'form-select form-select-sm'})
    )
    
    search = forms.CharField(
        required=False,
        widget=forms.TextInput(attrs={
            'class': 'form-control form-control-sm',
            'placeholder': 'Search by sample name...'
        })
    )
