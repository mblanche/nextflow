#!/usr/bin/env bash

docker build -t mblanche/$(basename ${PWD}) . \
    && docker push mblanche/$(basename ${PWD}) 
