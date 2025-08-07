# MixedModelsSmallSample.jl Development Instructions

Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.

## Package Overview

MixedModelsSmallSample.jl provides small sample corrections for confidence intervals and hypothesis tests for fixed effects in mixed models. It implements Kenward-Roger and Satterthwaite-Welch adjustments for small sample inference.

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
  - **Known issue**: MixedModels.jl test data artifacts may fail to download due to network restrictions - this does not affect package functionality

### Core Dependencies Verification

Test that core dependencies load (without MixedModels full functionality):

```bash
julia --project=. -e 'using Distributions, LinearAlgebra, StatsAPI, StatsBase; println("Core dependencies loaded successfully")'
```

  - **Expected time: 5-10 seconds**

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
  - Uses Blue style formatting (see `.JuliaFormatter.toml`)

### Code Quality Checks

Install and run code quality tools:

```bash
# Install quality tools (once per environment)
julia -e 'using Pkg; Pkg.add("Aqua")'
julia -e 'using Pkg; Pkg.add("JET")'

# Basic syntax/structure validation
julia --project=. -e 'using JuliaFormatter; format("."); println("Formatting check completed")'
```

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
# NEVER CANCEL: Documentation setup takes 2-4 minutes with many downloads
cd docs && julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Build documentation (when MixedModels works)
cd docs && julia --project=. make.jl
```

  - **Expected time: 2-4 minutes for setup**
  - **Known issue**: Documentation building may fail due to MixedModels dependency issues

## Repository Structure

### Key Files and Directories

```
├── .github/
│   ├── workflows/          # CI/CD pipelines
│   │   ├── CI.yml          # Main testing (Julia 1.10/1.11, Ubuntu/Windows)
│   │   ├── Format.yml      # Code formatting checks
│   │   ├── Documentation.yml # Docs building and deployment
│   │   └── spellcheck.yml  # Typo checking
├── src/
│   ├── MixedModelsSmallSample.jl  # Main module
│   └── show.jl             # Display methods
├── test/
│   ├── runtests.jl         # Test runner
│   ├── blocked experiment.jl
│   ├── split plot experiment.jl
│   ├── strip plot experiment.jl
│   ├── random slope.jl
│   ├── categorical.jl
│   ├── varia.jl
│   └── data/               # Test datasets (CSV files)
├── docs/                   # Documentation source
├── Project.toml            # Package definition
└── .JuliaFormatter.toml    # Code style configuration
```

### Core Package Functions

The package exports these main functions:

  - `adjust_KR()` - Kenward-Roger adjustment
  - `adjust_SW()` - Satterthwaite-Welch adjustment
  - `ftest_KR()` - F-tests with Kenward-Roger correction
  - `ftest_SW()` - F-tests with Satterthwaite-Welch correction
  - `vcov_varpar()` - Variance-covariance of variance parameters

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

### Adding New Functions

 1. Add function to `src/MixedModelsSmallSample.jl`
 2. Export function in module definition
 3. Add display methods to `src/show.jl` if needed
 4. Add tests to appropriate test file
 5. Update documentation in `docs/src/index.md`
 6. Run formatting: `julia -e 'using JuliaFormatter; format("."; verbose=true)'`

### Updating Dependencies

```bash
# Update Project.toml compatibility bounds
# Then regenerate Manifest.toml
julia --project=. -e 'using Pkg; Pkg.update()'
```

### Documentation Updates

 1. Edit `docs/src/index.md` for main documentation
 2. Use Julia docstrings in source code for API documentation
 3. Format markdown: `julia -e 'using JuliaFormatter; format("."; verbose=true)'`

## Validation Scenarios

### Minimum Validation (Network-Restricted Environments)

 1. Verify Julia version compatibility
 2. Test dependency installation (expect some failures)
 3. Run code formatting validation
 4. Check file structure and basic syntax

### Full Validation (Unrestricted Environments)

 1. Complete package installation
 2. Load package: `using MixedModelsSmallSample`
 3. Run test suite: `julia --project=. -e 'using Pkg; Pkg.test()'`
 4. Build documentation
 5. Test example workflows from documentation

## Known Issues and Workarounds

### Network Connectivity Issues

  - **Problem**: Cannot download MixedModels test data artifacts
  - **Symptom**: `Unable to automatically download/install artifact 'TestData'`
  - **Workaround**: Focus on code structure, formatting, and documentation changes
  - **Solution**: Use CI workflows for full testing

### Build Time Expectations

  - **Package installation**: 5-8 minutes (with potential network timeouts)
  - **Code formatting**: 2-3 seconds
  - **Documentation setup**: 2-4 minutes
  - **NEVER CANCEL**: Always wait for operations to complete

### Dependency Management

  - Julia 1.10+ required
  - MixedModels.jl is the primary dependency (may have network issues)
  - Test dependencies: Aqua.jl, JET.jl, SafeTestsets.jl
  - Documentation: Documenter.jl

## Quick Reference Commands

```bash
# Check Julia version
julia --version

# Quick dependency check
julia --project=. -e 'using Distributions, LinearAlgebra; println("Basic deps OK")'

# Format code (always run before commit)
julia -e 'using JuliaFormatter; format("."; verbose=true)'

# Check git status
git status

# View package status
julia --project=. -e 'using Pkg; Pkg.status()'
```

Remember: This package focuses on statistical computations for mixed models. When making changes, prioritize mathematical correctness, proper error handling, and clear documentation of small sample correction methods.
