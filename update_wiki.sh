#!/bin/bash
# Script to update the wiki with converted math notation

set -e

echo "Updating wiki repository with converted math notation..."

# Check if wiki directory exists
if [ ! -d "1D_SPH.wiki" ]; then
    echo "Cloning wiki repository..."
    git clone https://github.com/YuriOku/1D_SPH.wiki.git
fi

# Copy converted files
echo "Copying converted files..."
cp wiki-updates/*.md 1D_SPH.wiki/

# Remove the README that we added for documentation
rm -f 1D_SPH.wiki/README.md

# Commit and push changes
cd 1D_SPH.wiki

echo "Checking for changes..."
if git diff --quiet; then
    echo "No changes to commit"
else
    echo "Committing changes..."
    git add *.md
    git commit -m "Convert math notation to GitHub native LaTeX style

- Changed from old render.githubusercontent.com image URLs
- Now using GitHub's native $...$ (inline) and $$...$$ (display) math notation
- This makes formulas render directly in the wiki without external dependencies"
    
    echo "Pushing changes..."
    git push
    
    echo "Wiki successfully updated!"
fi
