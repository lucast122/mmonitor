import os
import json
from django.core.management.base import BaseCommand
from django.contrib.auth.models import User
from users.models import MAG
from dashboard.models import KeggMap, Annotation

class Command(BaseCommand):
    help = 'Ingest pipeline results into the database'

    def add_arguments(self, parser):
        parser.add_argument('--results-dir', type=str, required=True,
                          help='Directory containing pipeline results')

    def handle(self, *args, **options):
        results_dir = options['results_dir']
        
        # Iterate through MAG directories
        for mag_dir in os.listdir(results_dir):
            mag_path = os.path.join(results_dir, mag_dir)
            if not os.path.isdir(mag_path):
                continue

            # Read checkm2 results
            checkm2_path = os.path.join(mag_path, f"{mag_dir}.checkm2.json")
            if os.path.exists(checkm2_path):
                with open(checkm2_path) as f:
                    checkm2_data = json.load(f)
                
                # Create or update MAG record
                mag, created = MAG.objects.update_or_create(
                    name=mag_dir,
                    defaults={
                        'fasta_path': os.path.join(mag_path, f"{mag_dir}.fasta"),
                        'gff_path': os.path.join(mag_path, f"{mag_dir}.gff"),
                        'completeness': checkm2_data.get('completeness', 0.0),
                        'contamination': checkm2_data.get('contamination', 0.0),
                        'strain_heterogeneity': checkm2_data.get('strain_heterogeneity')
                    }
                )

                # Process annotation summary
                annotation_summary = os.path.join(mag_path, f"{mag_dir}.annotation_summary.json")
                if os.path.exists(annotation_summary):
                    with open(annotation_summary) as f:
                        ann_data = json.load(f)
                    
                    Annotation.objects.update_or_create(
                        mag=mag,
                        defaults={
                            'gene_count': ann_data.get('gene_count', 0),
                            'trna_count': ann_data.get('trna_count', 0),
                            'rrna_count': ann_data.get('rrna_count', 0),
                            'cds_count': ann_data.get('cds_count', 0),
                            'summary_path': annotation_summary
                        }
                    )

                # Process KEGG maps
                kegg_dir = os.path.join(mag_path, 'kegg_maps')
                if os.path.exists(kegg_dir):
                    for map_file in os.listdir(kegg_dir):
                        if map_file.endswith('.png') or map_file.endswith('.svg'):
                            map_id = os.path.splitext(map_file)[0]
                            KeggMap.objects.update_or_create(
                                mag=mag,
                                map_id=map_id,
                                defaults={
                                    'description': f"KEGG map {map_id}",
                                    'image_path': os.path.join(kegg_dir, map_file)
                                }
                            )

                status = 'Created' if created else 'Updated'
                self.stdout.write(self.style.SUCCESS(
                    f'{status} MAG record for {mag_dir}'))
