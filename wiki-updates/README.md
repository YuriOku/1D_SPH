# Wiki Math Notation Updates

This directory contains the updated wiki files with converted mathematical notation.

## What was changed

The wiki files have been converted from the old GitHub math rendering style:
```markdown
![latex](https://render.githubusercontent.com/render/math?math=...)
```

To the new GitHub native LaTeX math style:
- Inline math: `$...$`
- Display math: `$$...$$`

## Files converted

- `Formulation.md` - English formulation documentation
- `Kernel-function.md` - English kernel function documentation
- `Kernel-interpolation.md` - English kernel interpolation documentation
- `Volume-element.md` - English volume element documentation
- `定式化.md` - Japanese formulation documentation
- `カーネル関数.md` - Japanese kernel function documentation
- `カーネル補間法.md` - Japanese kernel interpolation documentation
- `体積要素.md` - Japanese volume element documentation

## How to apply these changes to the wiki

Since the GitHub wiki is a separate repository, these files need to be manually pushed to the wiki repository:

```bash
# Clone the wiki repository
git clone https://github.com/YuriOku/1D_SPH.wiki.git

# Copy the updated files
cp wiki-updates/*.md 1D_SPH.wiki/

# Commit and push to the wiki
cd 1D_SPH.wiki
git add *.md
git commit -m "Convert math notation to GitHub native LaTeX style"
git push
```

Alternatively, you can copy and paste the content from these files directly into the wiki editor on GitHub.
