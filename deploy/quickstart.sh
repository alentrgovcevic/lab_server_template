#!/bin/bash
# One-command demo: the plasmid map bank (public example plasmids) and the genome browser
# (Lotus japonicus chloroplast, built on first start) behind Caddy on https://localhost:8443.
# Needs Docker and python3.
#
#   deploy/quickstart.sh up          create deploy/.env (random secret), build, start, make an admin
#   deploy/quickstart.sh user NAME   create/reset an account (prints a random password)
#   deploy/quickstart.sh admin NAME  make an account admin (the admin page is /admin)
#   deploy/quickstart.sh down        stop and remove the stack (keeps the state volume)
#
# For your own data, edit deploy/.env (BANK_MOUNT, OLIGO_FILE, SITE_ADDRESS); genomes are uploaded
# on the genome browser's page.
set -euo pipefail
cd "$(dirname "$0")/.."

dc() { docker compose -f deploy/docker-compose.yml --env-file deploy/.env "$@"; }
pw() { python3 -c "import secrets; print(secrets.token_urlsafe(12))"; }
envval() { sed -n "s/^$1=//p" deploy/.env | tail -1; }

case "${1:-}" in
  up)
    command -v docker >/dev/null || { echo "Docker is needed: https://docs.docker.com/get-docker/"; exit 1; }
    if [ ! -f deploy/.env ]; then
      secret=$(python3 -c "import secrets; print(secrets.token_urlsafe(48))")
      sed "s|^PLASMIDBANK_SECRET=.*|PLASMIDBANK_SECRET=$secret|" deploy/.env.example > deploy/.env
      echo "created deploy/.env"
    fi
    mkdir -p deploy/certs
    dc build app genomes
    dc --profile tools run --rm scan >/dev/null   # index the bank once; re-run after adding maps
    dc up -d app genomes caddy
    dc restart app >/dev/null                      # pick up a fresh index if the app was already up
    port=$(envval HTTPS_PORT); site=$(envval SITE_ADDRESS)
    url="https://$site${port:+:$port}"; [ "$port" = 443 ] && url="https://$site"
    printf "waiting for the app"
    for _ in $(seq 1 60); do
      [ "$(curl -sk -o /dev/null -w '%{http_code}' "$url/healthz")" = 200 ] && break
      printf "."; sleep 2
    done; echo
    if dc exec -T app python -m plasmidbank.app.auth list 2>/dev/null | grep -q '(no users)'; then
      p=$(pw)
      dc exec -T -e PLASMIDBANK_NEW_PASSWORD="$p" app python -m plasmidbank.app.auth add admin >/dev/null
      dc exec -T app python -m plasmidbank.app.auth role admin admin >/dev/null
      echo "first account:  admin / $p   (change it after signing in)"
    fi
    echo "plasmid bank:   $url/            (self-signed certificate: accept the browser warning)"
    echo "genome browser: $url/genomes/"
    ;;
  user)
    name="${2:?usage: $0 user NAME}"; p=$(pw)
    dc exec -T -e PLASMIDBANK_NEW_PASSWORD="$p" app python -m plasmidbank.app.auth add "$name"
    echo "user: $name   password: $p"
    ;;
  admin)
    dc exec -T app python -m plasmidbank.app.auth role "${2:?usage: $0 admin NAME}" admin
    ;;
  down)
    dc down
    ;;
  *) sed -n 2,10p "$0"; exit 2 ;;
esac
