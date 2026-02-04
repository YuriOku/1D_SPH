# Wiki Math Notation Conversion - Summary

## Overview
This PR converts all mathematical notation in the GitHub Wiki from the old image-based rendering style to GitHub's native LaTeX math notation.

## Problem
The wiki was using the old GitHub math rendering format with external image URLs:
```markdown
![equation](https://render.githubusercontent.com/render/math?math=...)
```

This approach has several issues:
- Relies on external rendering service
- Makes formulas harder to read in source
- Doesn't render properly in some contexts
- URL-encoded LaTeX is difficult to maintain

## Solution
Converted all math expressions to GitHub's native LaTeX notation:
- **Inline math**: `$expression$`
- **Display math**: `$$expression$$`

## Files Converted
All 8 wiki documentation files (4 English, 4 Japanese):

### English
1. **Formulation.md** - SPH time evolution equations
2. **Kernel-function.md** - Kernel function documentation
3. **Kernel-interpolation.md** - Kernel interpolation methods
4. **Volume-element.md** - Volume element calculations

### Japanese (日本語)
1. **定式化.md** - SPH時間発展方程式
2. **カーネル関数.md** - カーネル関数の説明
3. **カーネル補間法.md** - カーネル補間法
4. **体積要素.md** - 体積要素の計算

## Examples of Conversion

### Before (Old Style)
```markdown
![v](https://render.githubusercontent.com/render/math?math=%5Ctextstyle+v)
```

### After (New Style)
```markdown
$v$
```

### Before (Display Math)
```markdown
![\begin{align*}
\frac{dv}{dt} &= -\frac{\nabla P}{\rho}
\end{align*}](https://render.githubusercontent.com/render/math?math=...)
```

### After (Display Math)
```markdown
$$
\begin{align*}
\frac{dv}{dt} &= -\frac{\nabla P}{\rho}
\end{align*}
$$
```

## Benefits
- ✅ Cleaner, more maintainable source code
- ✅ Renders natively in GitHub without external dependencies
- ✅ Better performance (no image loading)
- ✅ Works consistently across GitHub interface
- ✅ LaTeX source is directly visible and editable

## How to Apply

The converted files are available in the `wiki-updates/` directory. 

### Option 1: Automatic Script
Run the provided script from the repository root:
```bash
./update_wiki.sh
```

### Option 2: Manual Update
1. Navigate to the wiki repository
2. Copy files from `wiki-updates/` to the wiki
3. Commit and push changes

The wiki repository has already been updated with commit hash: `7537f68`

## Validation
All conversions have been verified to:
- Preserve mathematical accuracy
- Maintain proper LaTeX syntax
- Convert both inline and display math correctly
- Handle complex multi-line equations
- Support both English and Japanese documentation
