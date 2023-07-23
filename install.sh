#!/bin/bash
# takes two paramters, the domain name and the email to be associated with the certificate
DOMAIN=$1
EMAIL=$2

echo DOMAIN=${DOMAIN} >> .env
echo EMAIL=${EMAIL} >> .env

# Phase 1
docker compose -f ./docker-compose-initiate.yaml up -d nginx
docker compose -f ./docker-compose-initiate.yaml up certbot
docker compose -f ./docker-compose-initiate.yaml down

# some configurations for let's encrypt
curl -L --create-dirs -o etc/letsencrypt/options-ssl-nginx.conf https://raw.githubusercontent.com/certbot/certbot/master/certbot-nginx/certbot_nginx/_internal/tls_configs/options-ssl-nginx.conf
if [ ! -f etc/letsencrypt/ssl-dhparams.pem ]; then
  openssl dhparam -out etc/letsencrypt/ssl-dhparams.pem 2048
fi

# Phase 2
crontab ./etc/crontab
docker compose -f ./docker-compose.yaml up -d
