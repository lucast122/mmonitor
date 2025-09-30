from django.db.models.signals import post_save
from django.contrib.auth.models import User
from django.core.mail import send_mail
from django.conf import settings
from django.dispatch import receiver
from .models import UserQuota

@receiver(post_save, sender=User)
def create_user_quota(sender, instance, created, **kwargs):
    """Create UserQuota when a new user is created"""
    if created:
        UserQuota.objects.create(user=instance)

# Email sending temporarily disabled
# @receiver(post_save, sender=User)
# def send_email_to_admin_on_new_user(sender, instance, created, **kwargs):
#     if created:
#         subject = f'MMonitor - New User Created {instance.username}'
#         message = f'A new user has been created: {instance.username}'
#         email_from = settings.EMAIL_HOST_USER
#         recipient_list = ['timo-niklas.lucas@uni-tuebingen.de']  # Replace with your email
#         send_mail(subject, message, email_from, recipient_list)
