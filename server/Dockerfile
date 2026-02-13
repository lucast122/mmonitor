# Use the official Python image as a base
FROM python:3.12-slim

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Set work directory
WORKDIR /app

# Copy the server application code
COPY server/ /app/

# Copy frontend build output (built externally before docker build)
COPY frontend/dist/ /frontend/dist/

# Install dependencies
RUN apt-get update && apt-get install -y \
    pkg-config \
    gcc \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

# Collect static files
RUN python manage.py collectstatic --noinput

# Expose the port the app runs on
EXPOSE 8025

# Command to run the application
CMD ["gunicorn", "--workers", "3", "--bind", ":8025", "MMonitor.wsgi:application"]
