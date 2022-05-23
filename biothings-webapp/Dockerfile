FROM node:18 AS build-webapp

ARG WEBAPP_GIT=https://github.com/biothings/biothings_studio.git
ARG WEBAPP_VERSION=master

WORKDIR /build/src/github.com/biothings/biothings_studio
RUN git clone ${WEBAPP_GIT} .
RUN git checkout ${WEBAPP_VERSION}
WORKDIR /build/src/github.com/biothings/biothings_studio/webapp
RUN npm install && npm run build --legacy-peer-deps

FROM nginx:mainline-alpine
COPY --from=build-webapp --chown=root:www-data /build/src/github.com/biothings/biothings_studio/webapp/dist /srv/www/webapp
COPY webapp.conf /etc/nginx/conf.d/default.conf
RUN mkdir -p /var/log/nginx && \
    ln -svf /dev/stdout /var/log/nginx/access.log && \
    ln -svf /dev/stderr /var/log/nginx/error.log

EXPOSE 80/tcp
