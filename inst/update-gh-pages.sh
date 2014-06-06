#!/bin/sh
DOCS_DIR=web
GIT_DOCS="--git-dir=$DOCS_DIR/.git --work-tree=$DOCS_DIR"
VERSION=$(git rev-parse --short HEAD)

# Set up, first time only:
#
#   cd inst
#   git clone git@github.com:<user>/<repo>.git web
#   cd web
#   git checkout --orphan gh-pages
#   git rm -rf .
#
# After first time will just be a clone branch.
#
# Scripts to support this would be nice, but are more general than
# this project so may get factored out.  Hadley probably has some
# already!
#
# (see wood/.update-gh-pages.sh for a "replace pages" workflow).

# TODO: Deal with removals, or check that the commands below do that.
git $GIT_DOCS add .
git $GIT_DOCS commit -m "Generated static docs @ $SHORT_REF"
git $GIT_DOCS push origin gh-pages
