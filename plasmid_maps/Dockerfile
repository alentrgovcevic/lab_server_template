# Plasmid map bank: the web app and the scanner in one image.
# Build:  docker build -t plasmidbank .        (from plasmid_maps/; compose does this for you)
# Run:    see ../deploy/docker-compose.yml at the repo root (the app is never published directly; Caddy fronts it).
# Stage 1: build wheels. primer3-py has no prebuilt wheel on some platforms (e.g. arm64 Linux),
# so a compiler is needed here, and only here.
FROM python:3.13-slim AS wheels
RUN apt-get update \
 && apt-get install -y --no-install-recommends build-essential \
 && rm -rf /var/lib/apt/lists/*
COPY requirements.txt .
RUN pip wheel --no-cache-dir --wheel-dir /wheels -r requirements.txt

# Stage 2: runtime, no compiler.
FROM python:3.13-slim

# MAFFT for the Align tab, BLAST+ for "BLAST this"; nothing else from apt.
# Debian 13 ships BLAST+ 2.16; the Mac uses Homebrew's 2.17. Same database format (v5).
RUN apt-get update \
 && apt-get install -y --no-install-recommends mafft ncbi-blast+ \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY requirements.txt .
COPY --from=wheels /wheels /wheels
RUN pip install --no-cache-dir --no-index --find-links /wheels -r requirements.txt && rm -rf /wheels

COPY plasmidbank/ plasmidbank/

# Unprivileged user; /state holds the index, mirror, per-user outboxes, accounts and reports.
RUN useradd --uid 10001 --create-home plasmidbank \
 && mkdir -p /state /data/bank /data/oligos \
 && chown -R plasmidbank /state
USER plasmidbank

ENV PYTHONUNBUFFERED=1 \
    PLASMIDBANK_HOST=0.0.0.0 \
    PLASMIDBANK_PORT=8765 \
    PLASMIDBANK_INDEX=/state/index.sqlite \
    PLASMIDBANK_ROOT=/data/bank \
    PLASMIDBANK_OLIGOS=/data/oligos/oligo_export.tab \
    PLASMIDBANK_OUTBOX=/state/outbox \
    PLASMIDBANK_UPLOADS=/state/uploads \
    PLASMIDBANK_USERS=/state/users.json \
    PLASMIDBANK_STATE=/state \
    PLASMIDBANK_AUTH=local

EXPOSE 8765
HEALTHCHECK --interval=30s --timeout=5s --start-period=120s --retries=3 \
  CMD python -c "import urllib.request,sys; sys.exit(0 if urllib.request.urlopen('http://127.0.0.1:8765/healthz', timeout=4).status == 200 else 1)"

CMD ["python", "-m", "plasmidbank.app.server"]
