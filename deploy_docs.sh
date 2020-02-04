# Modified from https://github.com/AMReX-Codes/amrex/blob/master/build_and_deploy.sh
# and licensed according to BSD-3-Clause-LBNL (https://github.com/AMReX-Codes/amrex/blob/master/LICENSE).

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
git clone $REPO out
cd out

# Create a new empty branch if gh-pages doesn't exist yet (should only happen on first deploy)
git checkout $TARGET_BRANCH || git checkout --orphan $TARGET_BRANCH

# Clean out existing contents
git rm -rf . || exit 0
cd ..

# Pull from SOURCE_BRANCH again
git pull $SSH_REPO $SOURCE_BRANCH

# Build the Sphinx documentation
cd docs
make html
cd ../  

# Move it to the gh-pages branch
mv -f docs/_build/html/* out/
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