#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

# Build the documentation from the SOURCE_BRANCH
# and push it to TARGET_BRANCH.
SOURCE_BRANCH="master"
TARGET_BRANCH="gh-pages"

# Save some useful information
REPO=`git config remote.origin.url`
SSH_REPO=${REPO/https:\/\/github.com\//git@github.com:}
SHA=`git rev-parse --verify HEAD`

# Clone the existing gh-pages for this repo into out/
# Create a new empty branch if gh-pages doesn't exist yet (should only happen on first deply)
git clone $REPO out
cd out

git checkout $TARGET_BRANCH || git checkout --orphan $TARGET_BRANCH
cd ..

# Clean out existing contents
rm -rf out/sphinx_docs/ || exit 0

# Pull from SOURCE_BRANCH again
git pull origin $SOURCE_BRANCH

# Build the Sphinx documentation
cd docs
make html
cd ../  

mkdir -p out/sphinx_docs/
mv -f docs/_build/html/* out/sphinx_docs/
touch out/.nojekyll

# Now let's go have some fun with the cloned repo
cd out
git config --local user.name "Azure Pipelines"
git config --local user.email "azuredevops@microsoft.com"

echo "doing git add/commit/push"

# Commit the "changes", i.e. the new version.
# The delta will show diffs between new and old versions.
git add --all

# Exit if there are no docs changes
if git diff --staged --quiet; then
   echo "exiting with no docs changes"
   exit 0
fi

# Otherwise, commit and push
git commit -m "Deploy to GitHub Pages: ${SHA}"
git push $SSH_REPO $TARGET_BRANCH
cd ..

# See https://github.com/AMReX-Codes/amrex/blob/master/LICENSE for the license
# deploy_docs.sh was modified from the original script in the AMReX code (https://github.com/AMReX-Codes/amrex), 
# which is covered by the following license