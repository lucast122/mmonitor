from django.core.management.base import BaseCommand
from django.contrib.auth import get_user_model
from django.conf import settings

class Command(BaseCommand):
    help = 'Sets up an offline superuser'

    def handle(self, *args, **kwargs):
        User = get_user_model()
        if not User.objects.filter(username='offlinemode').exists():
            User.objects.create_superuser('offlinemode', 'offline@example.com', 'offlinemode')
            self.stdout.write(self.style.SUCCESS('Successfully created offline superuser'))
