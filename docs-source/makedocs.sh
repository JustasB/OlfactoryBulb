# Cause a rebuild of all doc files
touch *.rst

# the sphinx make script assumes 'html' folder instead of 'docs', in which github pages look
# it also expects the doc source files to be under /html/source subfolder
# So:
#   we copy the files under docs-source to the /html/source folder
#   run the make file to generate the html files
#   clean up by restoring the original folder names e.g. docs (needed for github pages)
sudo chmod 777 -R ../docs || true
sudo chmod 777 -R ../html || true

# Nuke old /docs or /html
echo Nuking old [repo]/docs and [repo]/html
rm -rf ../docs || true
rm -rf ../html || true

# Recreate it
echo Creating folder structure for Sphinx
mkdir ../docs
mkdir ../docs/source
cp -a . ../docs/source
mv ../docs/source/CNAME ../docs # This allows docs.olfactorybulb.org to point to the github pages
mv ../docs/source/Makefile ../docs
mv ../docs/source/.nojekyll ../docs # tells github that the html files are pre-generated
mv ../docs ../html
sudo chmod 777 -R ../html
cd ../html # We're now in [repo]/html

echo Now in:
pwd
echo Which contains:
ls -a

# Build the sphinx docker container (skips if have already)
echo Building documentation container...
docker build ../docker/documentation/. --quiet -t ob-sphinx-docs:latest

# Run docker sphix container to build the docs
echo Generating documentation HTML files...
docker run \
    -v $(readlink -f .):/html \
    -v $(readlink -f ../olfactorybulb):/olfactorybulb \
    -v $(readlink -f ../prev_ob_models):/prev_ob_models \
    ob-sphinx-docs:latest


# Cleanup
echo Preparing folder structure for GitHub pages \(i.e. all in [repo]/docs\)
cd ..

# rename /html into /docs, as expected by github pages
mv html docs

# remove redundant doc source files from /docs/source (already in /docs-source)
rm -rf docs/source

# Get back to the original /docs-source
cd docs-source
