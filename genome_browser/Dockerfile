# Genome browser. Built and run by ../deploy/docker-compose.yml.
# Pure Python: FastAPI + stdlib sqlite3, so no compiler and no apt packages are needed.
#
# The image carries no genomes. /data/genomes is a volume: one folder per genome, built from uploads
# (or from GENOMES_IMPORT_DIR) by `python -m genomebrowser.build`. /seed may hold a demo (seed.json).
FROM python:3.13-slim

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY genomebrowser/ genomebrowser/

RUN useradd --uid 10002 --create-home genomes && mkdir -p /data/genomes && chown genomes /data/genomes
USER genomes

ENV PYTHONUNBUFFERED=1 \
    GENOMES_HOST=0.0.0.0 \
    GENOMES_PORT=8766 \
    GENOMES_DIR=/data/genomes \
    SQLITE_TMPDIR=/data/genomes \
    GENOMES_AUTH=header

EXPOSE 8766
HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
  CMD python -c "import urllib.request,sys; sys.exit(0 if urllib.request.urlopen('http://127.0.0.1:8766/healthz', timeout=4).status == 200 else 1)"

CMD ["python", "-m", "genomebrowser.server"]
