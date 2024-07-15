#!/bin/bash

if ! compgen -G "${PREFIX}/lib/graphviz"/config* > /dev/null; then
    "${PREFIX}/bin"/dot -c || true
fi
