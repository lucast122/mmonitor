from django.contrib import admin

from .models import Feedback
from .models import Metadata
# Register your models here.
from .models import NanoporeRecord
from .models import SequencingStatistics, MAG

admin.site.register(NanoporeRecord)
admin.site.register(SequencingStatistics)
admin.site.register(Metadata)
admin.site.register(Feedback)
admin.site.register(MAG)
