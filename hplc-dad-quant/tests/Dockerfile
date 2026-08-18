FROM python:3.13-slim-bookworm

# Everything the verifier needs is installed here; test.sh never installs.
RUN pip install --no-cache-dir \
    pytest==9.1.1 \
    pytest-json-ctrf==0.5.2

# Bake the tests and the ground truth into the image - harbor does not upload
# tests/ at verify time.
COPY . /tests/

# The declared artifact lands at /app/results.json, so the directory has to exist.
RUN mkdir -p /app /logs/verifier

WORKDIR /tests
