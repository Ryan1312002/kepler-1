FROM python:3.13-slim-bookworm

# The task is pure data analysis, but a chromatographer would expect a working
# scientific stack to be available without having to fetch one first.
RUN pip install --no-cache-dir numpy==2.5.2 scipy==1.18.0 pandas==3.0.5 matplotlib==3.11.0

WORKDIR /app

# The acquisition data: 23 chromatograms, the run manifest and the analyst's
# method notes.
COPY data /app/data
