# Cause a rebuild of all doc files
touch *.rst

# the make script assumes 'html' folder instead of 'docs'
# it also expects the doc source files to be under /docs/source
# So:
#   we copy the doc source files to the source folder
#   then, rename docs to html
#   run the make file to generate the html files
#   clean up by restoring the original folder names e.g. docs
sudo chmod 777 -R ../docs
rm -Rf ../docs
mkdir ../docs
mkdir ../docs/source
cp -a . ../docs/source
mv ../docs/source/Makefile ../docs
mv ../docs/source/.nojekyll ../docs
mv ../docs ../html
sudo chmod 777 -R ../html
cd ../html # We're now in [repo]/html

#exit 1

# Build the sphinx docker container (skips if have already)
docker build ../docker/documentation/. -t ob-sphinx-docs:1.0

# Run docker sphix container to build the docs
#docker run -v $(readlink -f .):/html -v $(readlink -f ../olfactorybulb):/olfactorybulb ob-sphinx-docs:1.0

## Cleanup
#cd ..
#mv html docs
#rm -rf docs/source
#cd docs-source
