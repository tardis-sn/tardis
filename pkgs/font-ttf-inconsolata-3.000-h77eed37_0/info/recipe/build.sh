#!/bin/bash

mkdir -p ${PREFIX}/fonts || true
install -v -m644 ./fonts/ttf/Inconsolata-Regular.ttf ${PREFIX}/fonts/Inconsolata-Regular.ttf
install -v -m644 ./fonts/ttf/Inconsolata-Bold.ttf ${PREFIX}/fonts/Inconsolata-Bold.ttf
