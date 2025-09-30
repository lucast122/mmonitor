import base64
import json

from django.conf import settings
from django.contrib import messages
from django.contrib.auth import authenticate
from django.contrib.auth import login, logout
from django.contrib.auth.decorators import login_required
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User
from django.core.serializers import serialize
from django.http import JsonResponse, HttpResponseRedirect
from django.shortcuts import render, redirect
from django.views.decorators.csrf import csrf_exempt
from django.urls import reverse

from .forms import FeedbackForm
from .models import NanoporeRecord, SequencingStatistics, UserProfile, MAG
from rest_framework.decorators import api_view
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from rest_framework.authentication import BasicAuthentication
from rest_framework.permissions import IsAuthenticated
from django.db import transaction
from .models import MAGAnnotation
import logging

logger = logging.getLogger(__name__)

@csrf_exempt
def add_sequencing_statistics(request):
    if request.method == "POST":
        # Get the username and password from the Authorization header
        auth_header = request.headers.get('Authorization')
        if not auth_header or not auth_header.startswith('Basic '):
            return JsonResponse({"error": "Missing or invalid Authorization header"}, status=400)
        username, password = base64.b64decode(auth_header.split(' ', 1)[1]).decode().split(':', 1)

        # Authenticate the user
        user = authenticate(request, username=username, password=password)
        if user is None:
            return JsonResponse({"error": "Invalid username or password"}, status=400)

        print("POST request received at /users/add_sequencing_statistics/")
        try:
            data = json.loads(request.body.decode('utf-8'))
            SequencingStatistics.objects.create(
                sample_name=data['sample_name'],
                project_id=data['project_id'],
                subproject_id=data['subproject_id'],
                date=data['date'],
                mean_read_length=data.get('mean_read_length'),
                median_read_length=data.get('median_read_length'),
                mean_quality_score=data.get('mean_quality_score'),
                mean_gc_content=data.get('mean_gc_content'),
                read_lengths=data.get('read_lengths', "[]"),  # Default to empty list if not provided
                avg_qualities=data.get('avg_qualities', "[]"),  # Default to empty list if not provided
                number_of_reads=data.get('number_of_reads'),
                total_bases_sequenced=data.get('total_bases_sequenced'),
                q20_score=data.get('q20_score'),
                q30_score=data.get('q30_score'),
                avg_quality_per_read=json.dumps(data.get('avg_quality_per_read', [])),
                base_quality_avg=json.dumps(data.get('base_quality_avg', {})),
                gc_contents_per_sequence=data.get('gc_contents_per_sequence', "[]"),
                # Default to empty list if not provided

                user=user
            )
            return JsonResponse({"status": "success"}, status=201)
        except Exception as e:
            return JsonResponse({"status": "error", "error": str(e)}, status=400)


@csrf_exempt
def add_nanopore_record(request):
    print(f"Request method: {request.method}")
    print(f"Request path: {request.path}")
    print(f"Request body: {request.body}")
    print(f"Request headers: {request.headers}")

    if request.method == "POST":
        # Get the username and password from the Authorization header
        auth_header = request.headers.get('Authorization')
        if not auth_header or not auth_header.startswith('Basic '):
            return JsonResponse({"error": "Missing or invalid Authorization header"}, status=400)
        username, password = base64.b64decode(auth_header.split(' ', 1)[1]).decode().split(':', 1)

        # Authenticate the user
        user = authenticate(request, username=username, password=password)
        if user is None:
            return JsonResponse({"error": "Invalid username or password"}, status=400)

        print("POST request received at /users/add_nanopore_record/")
        try:
            data = json.loads(request.body)
            print(f"Received data: {data}")
            user_id = data.get('user_id')
            print(f"User ID: {user_id}")
            user = User.objects.get(id=user_id)
            print(f"User: {user}")

            NanoporeRecord.objects.create(
                taxonomy=data.get('taxonomy'),
                tax_genus=data.get('tax_genus'),
                tax_family=data.get('tax_family'),
                tax_order=data.get('tax_order'),
                tax_class=data.get('tax_class'),
                tax_phylum=data.get('tax_phylum'),
                tax_superkingdom=data.get('tax_superkingdom'),
                tax_clade=data.get('tax_clade'),
                tax_subspecies=data.get('tax_subspecies'),
                abundance=data.get('abundance'),
                count=data.get('count'),
                sample_id=data.get('sample_id'),
                project_id=data.get('project_id'),
                subproject=data.get('subproject'),
                date=data.get('date'),
                user=user
            )
            print("Nanopore record created successfully")
            return JsonResponse({"success": True}, status=200)
        except User.DoesNotExist:
            print(f"No user found with ID {user_id}")
            return JsonResponse({"error": "User not found"}, status=400)
        except Exception as e:
            print(f"Error occurred bla: {e}")
            return JsonResponse({"error": str(e), "type": str(type(e))}, status=400)


@csrf_exempt
def get_user_id(request):
    if request.method == 'POST':
        try:
            # Try to get data from JSON body first
            data = json.loads(request.body.decode('utf-8'))
            username = data.get('username')
            password = data.get('password')
        except json.JSONDecodeError:
            # Fall back to form data
            username = request.POST.get('username')
            password = request.POST.get('password')
        
        # Special handling for offline mode
        if username == 'offlinemode' and password == 'offline123':
            # Create or get offline user
            user, created = User.objects.get_or_create(
                username='offlinemode',
                defaults={'is_active': True}
            )
            if created:
                user.set_password('offline123')
                user.save()
            
            # Generate a token for offline mode
            token = base64.b64encode(f"{username}:{password}".encode()).decode()
            return JsonResponse({'user_id': user.id, 'token': token}, status=200)
            
        user = authenticate(request, username=username, password=password)
        if user is not None:
            # Generate a token for authenticated users
            token = base64.b64encode(f"{username}:{password}".encode()).decode()
            return JsonResponse({'user_id': user.id, 'token': token}, status=200)
        else:
            return JsonResponse({'error': 'Invalid username or password.'}, status=400)
    else:
        return JsonResponse({'error': 'Invalid request method.'}, status=400)


@csrf_exempt
def overwrite_nanopore_record(request):
    if request.method == "POST":
        try:
            records = json.loads(request.body)
            if not isinstance(records, list):
                return JsonResponse({"error": "Expected a list of records"}, status=400)

            auth_header = request.headers.get('Authorization')
            if not auth_header or not auth_header.startswith('Basic '):
                return JsonResponse({"error": "Missing or invalid Authorization header"}, status=400)

            username, password = base64.b64decode(auth_header.split(' ', 1)[1]).decode().split(':', 1)
            user = authenticate(request, username=username, password=password)
            if user is None:
                return JsonResponse({"error": "Invalid username or password"}, status=400)

            # Collect unique sample_ids from incoming records
            unique_sample_ids = set(record['sample_id'] for record in records if 'sample_id' in record)

            # Delete existing records for each unique sample_id
            for sample_id in unique_sample_ids:
                NanoporeRecord.objects.filter(sample_id=sample_id, user=user).delete()

            # Process each record
            new_records = []
            for record_data in records:
                sample_id = record_data.get('sample_id')
                # Delete existing records for the given sample_id
                NanoporeRecord.objects.filter(sample_id=sample_id, user=user).delete()

                # Prepare new record
                new_record = NanoporeRecord(
                    taxonomy=record_data.get('taxonomy'),
                    tax_genus=record_data.get('tax_genus'),
                    tax_family=record_data.get('tax_family'),
                    tax_order=record_data.get('tax_order'),
                    tax_class=record_data.get('tax_class'),
                    tax_phylum=record_data.get('tax_phylum'),
                    tax_superkingdom=record_data.get('tax_superkingdom'),
                    tax_clade=record_data.get('tax_clade'),
                    tax_subspecies=record_data.get('tax_subspecies'),
                    abundance=record_data.get('abundance'),
                    count=record_data.get('count'),
                    sample_id=sample_id,
                    project_id=record_data.get('project_id'),
                    subproject=record_data.get('subproject'),
                    date=record_data.get('date'),
                    user=user
                )
                new_records.append(new_record)

            # Bulk insert new records
            NanoporeRecord.objects.bulk_create(new_records)

            return JsonResponse({"message": "Records created successfully"}, status=201)

        except json.JSONDecodeError:
            return JsonResponse({"error": "Invalid JSON format"}, status=400)
        except Exception as e:
            return JsonResponse({"error": str(e)}, status=500)
    else:
        return JsonResponse({"error": "Invalid request method"}, status=405)


@csrf_exempt
def get_unique_sample_ids(request):
    if request.method == 'POST':
        username = request.POST.get('username')
        password = request.POST.get('password')
        user = authenticate(username=username, password=password)

        if user is not None:
            unique_sample_ids = NanoporeRecord.objects.filter(user=user).values_list('sample_id', flat=True).distinct()
            sample_ids_list = list(unique_sample_ids)
            print(sample_ids_list)
            return JsonResponse({'sample_ids': sample_ids_list}, status=200)
        else:
            return JsonResponse({'error': 'Invalid username or password.'}, status=400)
    else:
        return JsonResponse({'error': 'Invalid request method.'}, status=400)


def create_mysql_user(username, password):
    from django.conf import settings
    import mysql.connector

    # Connect to MySQL as root
    cnx = mysql.connector.connect(
        user=settings.DATABASES['mmonitor']['USER'],
        password=settings.DATABASES['mmonitor']['PASSWORD'],
        host=settings.DATABASES['mmonitor']['HOST'],
        database=settings.DATABASES['mmonitor']['NAME']
    )
    cursor = cnx.cursor()

    # Create the MySQL user and grant privileges
    try:
        cursor.execute(f"CREATE USER '{username}'@'%' IDENTIFIED BY '{password}';")
        cursor.execute(f"GRANT SELECT, INSERT, UPDATE, DELETE ON your_database.* TO '{username}'@'%';")

    except mysql.connector.Error as err:
        print(f"Failed creating MySQL user: {err}")
        return False

    cnx.commit()
    cursor.close()
    cnx.close()
    return True


@login_required
def profile_view(request):
    # Fetch the UserProfile for the currently logged-in user
    user_profile = UserProfile.objects.get(user=request.user)
    # Now user_profile contains the data for the currently logged-in user only


def redirect_to_dash(request):
    # Get the Django session ID
    session_id = request.session.session_key

    # Redirect to the Dash app, including the session ID as a URL parameter
    return redirect(f'{settings.DASH_APP_URL}?session_id={session_id}')


def register(request):
    """Register a new user"""
    if request.method != "POST":
        # Display black registration form
        form = UserCreationForm()

    else:
        # Process completed form
        form = UserCreationForm(data=request.POST)

        if form.is_valid():
            new_user = form.save()
            # default to non-active so admin needs to authorize user
            new_user.is_active = False
            new_user.save()
            UserProfile.objects.create(user=new_user, some_field="default value")

            # Login and redirect to homepage
            login(request, new_user)
            return redirect("main:index")

    # Display a blank or invalid form
    context = {"form": form}

    return render(request, "registration/register.html", context)


# USER FEEDBACK
def submit_feedback(request):
    if request.method == 'POST':
        form = FeedbackForm(request.POST)
        if form.is_valid():
            feedback = form.save(commit=False)
            feedback.user = request.user
            feedback.save()
            messages.success(request, 'Your feedback has been submitted!')
            return redirect('users:feedback_thank_you')
    else:
        form = FeedbackForm()
    return render(request, 'submit_feedback.html', {'form': form})


def feedback_thank_you(request):
    return render(request, 'feedback_thank_you.html')


def chart_data_api(request):
    chart_data_queryset = NanoporeRecord.objects.all()
    chart_data_json = serialize('json', chart_data_queryset)
    print(chart_data_json)
    return JsonResponse(chart_data_json, safe=False)

#
# def dashboard_view(request):
#     # Prepare your data for the chart
#     chart_data_queryset = NanoporeRecord.objects.all()
#
#     # Serialize the queryset to JSON
#     chart_data_json = serialize('json', chart_data_queryset)
#
#     print("Serialized chart data: ", chart_data_json)
#
#     request.session['chart_data'] = chart_data_json
#     request.session.save()  # Ensure session data is saved
#
#     # Debug: Log the session ID and data after saving
#     print("Session ID: ", request.session.session_key)
#     print("Session Data: ", request.session['chart_data'])
#
#     # Render your template with the chart data
#     return render(request, 'base.html', {'chart_data': chart_data_json})


from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.contrib.auth import authenticate
import base64
import json

@csrf_exempt
def upload_mag(request):
    if request.method == "POST":
        # Get the username and password from the Authorization header
        auth_header = request.headers.get('Authorization')
        if not auth_header or not auth_header.startswith('Basic '):
            return JsonResponse({"error": "Missing or invalid Authorization header"}, status=400)
        
        try:
            username, password = base64.b64decode(auth_header.split(' ', 1)[1]).decode().split(':', 1)
        except Exception as e:
            return JsonResponse({"error": f"Invalid Authorization header format: {str(e)}"}, status=400)

        # Authenticate the user
        user = authenticate(request, username=username, password=password)
        if user is None:
            return JsonResponse({"error": "Invalid username or password"}, status=400)

        print("POST request received at /users/upload_mag/")
        try:
            # Debug: Print all request information
            print(f"POST parameters: {request.POST}")
            print(f"FILES parameters: {request.FILES}")
            print(f"Content type: {request.content_type}")
            
            # Initialize file content variables
            fasta_content = None
            gff_content = None
            fai_content = None
            completeness = None
            contamination = None
            strain_heterogeneity = None
            kegg_data = None
            
            # Get metadata from POST parameters
            name = request.POST.get('name')
            taxonomy = request.POST.get('taxonomy')
            sample_name = request.POST.get('sample_name')
            
            # Method 2: For standard POST fields (NOT accessing request.body)
            # For JSON payloads, Django won't populate request.POST, so we'll rely on FILES and other params
            if not name:
                name = request.POST.get('mag_name') or request.POST.get('bin_name')
            if not taxonomy:
                taxonomy = request.POST.get('taxonomy')
            if not sample_name:
                sample_name = request.POST.get('sample_id') or request.POST.get('sample')
            
            # Get quality metrics from POST
            completeness = request.POST.get('completeness') or request.POST.get('quality_completeness')
            contamination = request.POST.get('contamination') or request.POST.get('quality_contamination')
            strain_heterogeneity = request.POST.get('strain_heterogeneity') or request.POST.get('quality_strain_heterogeneity')
            
            # Try URL parameters as fallback
            if not name:
                name = request.GET.get('name') or request.GET.get('mag_name') or request.GET.get('bin_name')
            if not taxonomy:
                taxonomy = request.GET.get('taxonomy')
            if not sample_name:
                sample_name = request.GET.get('sample_name') or request.GET.get('sample_id') or request.GET.get('sample')
            
            # Get quality metrics from URL params if not in POST
            if not completeness:
                completeness = request.GET.get('completeness') or request.GET.get('quality_completeness')
            if not contamination:
                contamination = request.GET.get('contamination') or request.GET.get('quality_contamination')
            if not strain_heterogeneity:
                strain_heterogeneity = request.GET.get('strain_heterogeneity') or request.GET.get('quality_strain_heterogeneity')
            
            print(f"Final metadata values:")
            print(f"name: {name}")
            print(f"taxonomy: {taxonomy}")
            print(f"sample_name: {sample_name}")
            print(f"completeness: {completeness}")
            print(f"contamination: {contamination}")
            print(f"strain_heterogeneity: {strain_heterogeneity}")

            if not all([name, taxonomy, sample_name]):
                missing_fields = []
                if not name: missing_fields.append("name")
                if not taxonomy: missing_fields.append("taxonomy")
                if not sample_name: missing_fields.append("sample_name")
                return JsonResponse({
                    "error": f"Missing required metadata fields: {', '.join(missing_fields)}",
                    "received_data": {
                        "from_post": dict(request.POST),
                        "files": [k for k in request.FILES.keys()]
                    }
                }, status=400)

            # Try to get files from multipart form data
            fasta_file = request.FILES.get('fasta_file')
            if not fasta_file:
                fasta_file = request.FILES.get('fasta') or request.FILES.get('genome_file') or request.FILES.get('mag_file')
            if fasta_file:
                fasta_content = fasta_file.read().decode('utf-8')
            
            # GFF is optional for now in multipart uploads
            gff_file = request.FILES.get('gff_file')
            if not gff_file:
                gff_file = request.FILES.get('gff') or request.FILES.get('annotation_file')
            if gff_file:
                gff_content = gff_file.read().decode('utf-8')
                    
            # FAI is optional
            fai_file = request.FILES.get('fai_file')
            if not fai_file:
                fai_file = request.FILES.get('fai') or request.FILES.get('index_file')
            if fai_file:
                fai_content = fai_file.read().decode('utf-8')
            
            print(f"fasta_content: {'Present' if fasta_content else 'Missing'} - {len(fasta_content) if fasta_content else 0} bytes")
            print(f"gff_content: {'Present' if gff_content else 'Missing'} - {len(gff_content) if gff_content else 0} bytes")
            print(f"fai_content: {'Present' if fai_content else 'Missing'} - {len(fai_content) if fai_content else 0} bytes")

            # Only FASTA content is required, GFF is optional for now
            if not fasta_content:
                return JsonResponse({
                    "error": "Missing FASTA content in the request",
                    "received_data": {
                        "files": [k for k in request.FILES.keys()]
                    }
                }, status=400)
            
            # Check if MAG already exists for this user and name
            # If so, update it rather than creating a new one
            try:
                # First try to get the existing MAG
                mag = MAG.objects.get(name=name, user=user)
                print(f"MAG {name} already exists - updating")
                
                # Check if it's a chunk or full upload
                chunk_number = request.POST.get('chunk_number', 0)
                try:
                    chunk_number = int(chunk_number)
                except (ValueError, TypeError):
                    chunk_number = 0
                
                # Update the MAG with new data
                if chunk_number == 0:
                    # First chunk - replace content
                    mag.fasta_file = fasta_content
                    if gff_content:
                        mag.gff_file = gff_content
                    if fai_content:
                        mag.fai_file = fai_content
                else:
                    # Subsequent chunk - append to fasta content
                    mag.fasta_file += fasta_content
                
                # Update other fields if provided
                if completeness is not None:
                    mag.completeness = completeness
                if contamination is not None:
                    mag.contamination = contamination
                if strain_heterogeneity is not None:
                    mag.strain_heterogeneity = strain_heterogeneity
                if taxonomy:
                    mag.taxonomy = taxonomy
                if sample_name:
                    mag.sample_name = sample_name
                
                # Save the updated MAG
                mag.save()
                is_new = False
                
            except MAG.DoesNotExist:
                # Create new MAG if it doesn't exist
                print(f"Creating new MAG {name}")
                mag = MAG.objects.create(
                    name=name,
                    taxonomy=taxonomy,
                    sample_name=sample_name,
                    gff_file=gff_content or "",  # GFF content is optional for now
                    fasta_file=fasta_content,
                    fai_file=fai_content or "",
                    user=user,
                    # Add quality metrics if available
                    completeness=completeness if completeness is not None else None,
                    contamination=contamination if contamination is not None else None,
                    strain_heterogeneity=strain_heterogeneity if strain_heterogeneity is not None else None,
                    # Add annotations if provided
                    annotations_data=None
                )
                is_new = True
            
            print(f"MAG {mag.name} {'created' if is_new else 'updated'} successfully")
            return JsonResponse({"success": True, "mag_id": mag.id, "is_new": is_new}, status=201)
        
        except Exception as e:
            print(f"Error occurred: {e}")
            import traceback
            traceback.print_exc()
            return JsonResponse({"error": str(e), "traceback": traceback.format_exc()}, status=400)
    
    return JsonResponse({"error": "Invalid request method"}, status=405)

@api_view(['POST'])
def offline_login(request):
    """Create an offline user session"""
    try:
        # Create or get offline user
        user, created = User.objects.get_or_create(
            username='offline_user',
            defaults={'is_active': True}
        )
        
        # Create session
        login(request, user)
        
        return Response({'status': 'success'})
    except Exception as e:
        return Response({'error': str(e)}, status=500) 


@login_required
def mag_viewer(request):
    """
    View for the MAG analysis dashboard
    """
    return render(request, 'users/mag_viewer.html', {
        'title': 'MAG Analysis Dashboard'
    })

# Add a custom logout view that supports both GET and POST
def custom_logout(request):
    """
    Custom logout view that supports both GET and POST methods.
    This is necessary because Django's default LogoutView only allows POST.
    """
    logout(request)
    messages.success(request, "You have been successfully logged out.")
    return redirect('main:index')

class UploadMAGAnnotationsView(APIView):
    authentication_classes = [BasicAuthentication]
    permission_classes = [IsAuthenticated]
    
    def post(self, request, *args, **kwargs):
        """Handle POST request for uploading MAG annotations."""
        try:
            # Get data from request
            data = request.data
            sample_id = data.get('sample_id')
            annotations = data.get('annotations', [])
            
            # Validate required fields
            if not sample_id:
                return Response(
                    {"error": "sample_id is required"}, 
                    status=status.HTTP_400_BAD_REQUEST
                )
            
            if not annotations:
                return Response(
                    {"error": "annotations list is empty"}, 
                    status=status.HTTP_400_BAD_REQUEST
                )
            
            # Get the MAG instance
            try:
                mag = MAG.objects.get(name=sample_id, user=request.user)
            except MAG.DoesNotExist:
                return Response(
                    {"error": f"MAG with id {sample_id} not found"}, 
                    status=status.HTTP_404_NOT_FOUND
                )
            
            # Process annotations in a transaction
            with transaction.atomic():
                # Delete existing annotations for this MAG
                MAGAnnotation.objects.filter(sample=mag).delete()
                
                # Create new annotations
                annotations_to_create = []
                for annotation in annotations:
                    # Validate required fields
                    required_fields = ['gene_id', 'contig', 'start', 'stop', 'strand', 'product']
                    missing_fields = [field for field in required_fields if not annotation.get(field)]
                    
                    if missing_fields:
                        return Response(
                            {"error": f"Missing required fields: {', '.join(missing_fields)}"}, 
                            status=status.HTTP_400_BAD_REQUEST
                        )
                    
                    # Create annotation instance
                    mag_annotation = MAGAnnotation(
                        sample=mag,
                        gene_id=annotation['gene_id'],
                        contig=annotation['contig'],
                        start=annotation['start'],
                        stop=annotation['stop'],
                        strand=annotation['strand'],
                        product=annotation['product'],
                        kegg=annotation.get('kegg', ''),
                        ec_numbers=annotation.get('ec_numbers', []),
                        uniref=annotation.get('uniref', []),
                        sequence_nt=annotation.get('sequence_nt', ''),
                        sequence_aa=annotation.get('sequence_aa', '')
                    )
                    annotations_to_create.append(mag_annotation)
                
                # Bulk create annotations
                MAGAnnotation.objects.bulk_create(
                    annotations_to_create,
                    batch_size=1000  # Process in batches of 1000
                )
            
            # Return success response
            return Response({
                "message": f"Successfully uploaded {len(annotations)} annotations for MAG {sample_id}",
                "annotations_count": len(annotations)
            }, status=status.HTTP_200_OK)
            
        except Exception as e:
            logger.error(f"Error uploading annotations: {str(e)}")
            logger.error(traceback.format_exc())
            return Response(
                {"error": f"Server error: {str(e)}"}, 
                status=status.HTTP_500_INTERNAL_SERVER_ERROR
            )

