# Create a new release of LRphase and update all repositories. Takes one
# argument to indicate whether this is a major, minor, or patch release.

# Increment version numbers
bump2version --list $1

# Update PyyPI
python3 setup.py sdist bdist_wheel
python3 -m twine upload dist/LRphase-$(cat VERSION)*

# Bioconda is linked to GitHub repo, so no need to initiate updates here.
