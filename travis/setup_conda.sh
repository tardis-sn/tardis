#!/bin/bash

# Note to the future: keep the conda scripts separate for each OS because many
# packages call ci-helpers with:
#
#   source ci-helpers/travis/setup_conda_$TRAVIS_OS_NAME.sh
#
# The present script was added later.

if [[ $DEBUG == True ]]; then
    set -x
fi

# First check: if the build should be run at all based on the event type

if [[ ! -z $EVENT_TYPE ]]; then
    for event in $EVENT_TYPE; do
        if [[ $TRAVIS_EVENT_TYPE = $event ]]; then
            allow_to_build=True
        fi
    done
    if [[ $allow_to_build != True ]]; then
        travis_terminate 0
    fi
fi

# Second check: if any of the custom tags are used to skip the build

TR_SKIP="\[(skip travis|travis skip)\]"
DOCS_ONLY="\[docs only|build docs\]"

# Travis doesn't provide the commit message of the top of the branch for
# PRs, only the commit message of the merge. Thus this ugly workaround is
# needed for now.

if [[ $TRAVIS_PULL_REQUEST == false ]]; then
    COMMIT_MESSAGE=${TRAVIS_COMMIT_MESSAGE}
else
    COMMIT_MESSAGE=$(git show -s $TRAVIS_COMMIT_RANGE | awk 'BEGIN{count=0}{if ($1=="Author:") count++; if (count==1) print $0}')
fi

# Skip build if the commit message contains [skip travis] or [travis skip]
# Remove workaround once travis has this feature natively
# https://github.com/travis-ci/travis-ci/issues/5032

if [[ ! -z $(echo ${COMMIT_MESSAGE} | grep -E "${TR_SKIP}") ]]; then
    echo "Travis was requested to be skipped by the commit message, exiting."
    travis_terminate 0
elif [[ ! -z $(echo ${COMMIT_MESSAGE} | grep -E "${DOCS_ONLY}") ]]; then
    if [[ $SETUP_CMD != *build_docs* ]] && [[ $SETUP_CMD != *build_sphinx* ]]; then
        echo "Only docs build was requested by the commit message, exiting."
        travis_terminate 0
    fi
fi

echo "==================== Starting executing ci-helpers scripts ====================="

source travis/setup_conda_$TRAVIS_OS_NAME.sh;

echo "================= Returning executing local .travis.yml script ================="

set +x
