#!/bin/sh
github-changes -o tardis-sn -r tardis --only-pulls --use-commit-body -f tmp_changelog.md
pandoc --from markdown --to rst tmp_changelog.md > TMP_CHANGELOG.rst
