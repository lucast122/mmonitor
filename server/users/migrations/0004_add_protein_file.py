# Generated manually to fix protein_file issue

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('users', '0003_remove_mag_annotations_mag_annotations_data_and_more'),
    ]

    operations = [
        migrations.AddField(
            model_name='mag',
            name='protein_file',
            field=models.TextField(blank=True, null=True),
        ),
    ] 