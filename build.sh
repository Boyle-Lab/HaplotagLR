# Create a new release of HaplotagLR and update all repositories. Takes one
# argument to indicate whether this is a major, minor, or patch release.

# Increment version numbers
bump2version --verbose --current-version $(cat VERSION) --commit --tag --list $1

# Update PyyPI
python3 setup.py sdist bdist_wheel
python3 -m twine upload dist/HaplotagLR-$(cat VERSION)*

# Bioconda is linked to GitHub repo, so no need to initiate updates here as long as each version is tagged as a release in github.
# See https://widdowquinn.github.io/coding/update-bioconda-package/
