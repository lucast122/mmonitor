# from django.db import connections
from django.shortcuts import render
# from django.contrib.auth.decorators import login_required
# from django.contrib.admin.views.decorators import staff_member_required


def index(request):
    """Homepage"""
    return render(request, "main/index.html")
