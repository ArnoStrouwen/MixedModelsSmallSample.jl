# MixedModelsSmallSample.jl Development Instructions

Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.

## Package Overview

MixedModelsSmallSample.jl provides small sample corrections for confidence intervals and hypothesis tests for fixed effects in mixed models. It implements Kenward-Roger and Satterthwaite adjustments for small sample inference.

## Environment Setup

### Julia Installation

Run these commands to ensure Julia 1.10+ is available:

```bash
julia --version  # Should show Julia 1.10 or higher
```

### Package Dependencies Installation

```bash
# NEVER CANCEL: Initial package installation takes 5-8 minutes due to large dependency tree
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

  - **Expected time: 5-8 minutes**
  - **NEVER CANCEL**: Network timeouts are common, package will eventually install

## Development Workflow

### Code Formatting

Always run code formatting before committing changes:

```bash
# Install JuliaFormatter (once per environment)
julia -e 'using Pkg; Pkg.add("JuliaFormatter")'

# Format all code - FAST operation
julia -e 'using JuliaFormatter; format("."; verbose=true)'
```

  - **Expected time: 2-3 seconds**
  - **NEVER CANCEL**: This is a quick operation

### Testing Limitations

**CRITICAL**: Full package tests cannot run in environments with network restrictions because:

  - MixedModels.jl requires downloading test data artifacts from external sources
  - Network domains like `pkg.julialang.org` and `osf.io` may be blocked
  - This prevents `using MixedModelsSmallSample` from working in restricted environments

**Workaround**:

  - Test individual functions when possible
  - Use CI/CD workflows for full testing
  - Focus on code quality, formatting, and documentation changes

### Documentation Building

```bash
cd docs && julia --project=. -e 'using Pkg; Pkg.instantiate()'
cd docs && julia --project=. make.jl
```

  - **Expected time: 2-4 minutes for setup**

## CI/CD Workflows

### Important Timeout Settings

All workflows use appropriate timeouts:

  - **CI.yml**: 60 minutes timeout for full testing
  - **Format.yml**: Quick formatting validation
  - **Documentation.yml**: Full docs build and deployment
  - **spellcheck.yml**: Typo detection with crate-ci/typos

### Before Committing Changes

Always run these commands before committing:

```bash
# 1. Format code (required)
julia -e 'using JuliaFormatter; format("."; verbose=true)'

# 2. Check repository status
git status

# 3. Commit with descriptive message
git add .
git commit -m "Your descriptive commit message"
```

## Common Tasks
