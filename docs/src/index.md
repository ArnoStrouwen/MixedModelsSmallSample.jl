```@meta
CurrentModule = MixedModelsSmallSample
```

```@docs
MixedModelsSmallSample
```

## Getting Started

We will look at a split plot experiment, described in chapter 10 of:

> Goos, Peter, and Bradley Jones. Optimal design of experiments: a case study approach. John Wiley & Sons, 2011.

Split-plot experiments have hard-to-change factors and easy-to-change factors.
Hard-to-change factors are not properly randomized,
they are kept at the same factor setting for multiple observations.
The dataset can thus be split up into blocks, where the hard-to-change factors remain constant.
These blocks are also called whole plots.
The observations within such a block are more similar to one another than observations coming from a different block.
A random intercept per block is introduced to deal with the correlation between observations.
In the below dataset, FRH and RRH are hard-to-change, while YA and GC are easy-to-change.
The WP column shows that there are 10 blocks in the dataset in which FRH and RRH remain constant.
Efficiency is the output we want to model.

```@raw html
<details><summary>Click here to expand the details loading the data.</summary>
```

```@example split_plot
#! format: off
using DataFrames
df = DataFrame(
  WP = [1,1,1,1,1,
        2,2,2,2,2,
        3,3,3,3,3,
        4,4,4,4,4,
        5,5,5,5,5,
        6,6,6,6,6,
        7,7,7,7,7,
        8,8,8,8,8,
        9,9,9,9,9,
        10,10,10,10,10],
  FRH = [-1,-1,-1,-1,-1,
         0,0,0,0,0,
         1,1,1,1,1,
         -1,-1,-1,-1,-1,
         0,0,0,0,0,
         0,0,0,0,0,
         1,1,1,1,1,
         -1,-1,-1,-1,-1,
         0,0,0,0,0,
         1,1,1,1,1],
  RRH = [-1,-1,-1,-1,-1,
         -1,-1,-1,-1,-1,
         -1,-1,-1,-1,-1,
         0,0,0,0,0,
         0,0,0,0,0,
         0,0,0,0,0,
         0,0,0,0,0,
         1,1,1,1,1,
         1,1,1,1,1,
         1,1,1,1,1],
  YA = [-1,1,0,1,-1,
        1,0,1,0,-1,
        1,1,-1,-1,0,
        -1,0,1,0,-1,
        0,1,0,-1,0,
        1,0,0,0,-1,
        0,1,0,0,-1,
        -1,1,0,1,-1,
        -1,0,0,0,1,
        -1,1,-1,1,0],
  GC = [-1,1,0,-1,1,
        1,0,-1,1,0,
        -1,1,-1,1,0,
        0,1,0,-1,1,
        -1,0,0,1,0,
        0,0,0,0,-1,
        -1, 1,0,1,0,
        -1,-1,0,1,1,
        0,-1,0,1,1,
        1,-1,-1,0,1],
  EFFICIENCY = [
                0.873,0.969,0.959,0.821,1.019,
                0.921,0.923,0.821,0.958,0.914,
                0.705,0.861,0.756,0.861,0.78,
                0.999,1.06,0.91,0.854,1.082,
                0.841,0.885,0.905,1.025,0.936,
                0.844,0.891,0.902,0.871,0.833,
                0.786,0.879,0.86,0.905,0.87,
                0.98,0.907,1.045,1.094,1.159,
                1.0,0.901,0.99,1.059,1.025,
                0.991,0.828,0.853,0.923,0.997
                ]
  )
```

```@raw html
</details>
```

```@example split_plot
df[1:10, :]
```

```@example split_plot
df[45:50, :]
```

We work with a full quadratic model, where the intercept is allowed to vary between blocks.
First, we get the output using only `MixedModels`.
Here, the p-values for the main effects (first-order terms) are all below `1e-10`.

```@example split_plot
using MixedModels
fm = @formula(
    EFFICIENCY ~
        1 +
    (1 | WP) +
    FRH +
    RRH +
    YA +
    GC +
    FRH & RRH +
    FRH & YA +
    FRH & GC +
    RRH & YA +
    RRH & GC +
    YA & GC +
    FRH & FRH +
    RRH & RRH +
    YA & YA +
    GC & GC
)
m = fit(MixedModel, fm, df; REML=true)
```

Now we look at the adjustments made in the Kenward Roger approach.
The p-values for the main effects of the hard-to-change factors are much larger.
This is primarily because the degrees of freedom for these factors are very limited.

Similar results hold for the second-order terms.
Terms involving only the hard-to-change factors have limited degrees of freedom.
The interactions involving both a hard and an easy-to-change factor
have similar degrees of freedom as terms only involving easy-to-change factors.

```@example split_plot
using MixedModelsSmallSample
kr = small_sample_adjust(m, KenwardRoger(; fim=ExpectedFIM()))
```

The standard errors on the estimates are also slightly different.
However, for this specific experiment, the adjustment does not have a large effect.
An example where there is a large difference in standard error will be added later.

## Multiple random effects

```@example slope
using CSV
using DataFrames
using MixedModels

using MixedModelsSmallSample

df = DataFrame(MixedModels.dataset(:sleepstudy))
fm = @formula(reaction ~ 1 + days + (1 + days | subj))
m = fit(MixedModel, fm, df; REML=true)
```

```@example slope
kr = small_sample_adjust(m, KenwardRoger(; fim=ExpectedFIM())) # or ObservedFIM()
```

## More examples

### Scalar Random Effects

  - [Bioequivalence (Homogeneous)](https://github.com/ArnoStrouwen/MixedModelsSmallSample.jl/blob/master/test/bioequivalence%20homogeneous.jl)
  - [Blocked experiment](https://github.com/ArnoStrouwen/MixedModelsSmallSample.jl/blob/master/test/blocked%20experiment.jl)
  - [Split-plot experiment](https://github.com/ArnoStrouwen/MixedModelsSmallSample.jl/blob/master/test/split%20plot%20experiment.jl)
  - [Strip-plot experiment](https://github.com/ArnoStrouwen/MixedModelsSmallSample.jl/blob/master/test/strip%20plot%20experiment.jl)

### Vector Random Effects

  - [Bioequivalence (Intermediate)](https://github.com/ArnoStrouwen/MixedModelsSmallSample.jl/blob/master/test/bioequivalence%20intermediate.jl)
  - [Single random slope](https://github.com/ArnoStrouwen/MixedModelsSmallSample.jl/blob/master/test/random%20slope.jl)

## Similar software

  - The R package [pbkrtest](https://github.com/hojsgaard/pbkrtest)
  - [JMP](https://www.jmp.com/support/help/en/18.1/index.shtml#page/jmp/statistical-details-for-the-kackarharville-correction-2.shtml#)
  - [SAS](https://documentation.sas.com/doc/en/statcdc/14.2/statug/statug_glimmix_details40.htm)

## How It Works

The computations follow the same high-level flow for both Kenward-Roger and Satterthwaite:

 1. Compute ``V(\theta)``, ``V(\theta)^{-1}``, and derivative matrices for either [`Unstructured`](@ref) or [`FactorAnalytic`](@ref).
 2. Compute `W = I(θ)^{-1}` via [`vcov_varpar`](@ref), using either [`ObservedFIM`](@ref) or [`ExpectedFIM`](@ref).
 3. For [`KenwardRoger`](@ref), compute the bias-corrected covariance of `β` (see [`MixedModelsSmallSample.compute_adjusted_covariance`](@ref)).
 4. Compute denominator degrees of freedom and (for KR) the F-scaling factor (see [`MixedModelsSmallSample.compute_dof`](@ref)).

## API Reference

### Main Functions

```@docs
small_sample_adjust
small_sample_ftest
vcov_varpar
```

### Adjustment Methods

```@docs
KenwardRoger
Satterthwaite
```

### Options

```@docs
ObservedFIM
ExpectedFIM
Unstructured
FactorAnalytic
```

### Result Types

```@docs
LinearMixedModelAdjusted
```

## Internals

These types and functions are not exported but may be useful for understanding
or extending the package.

### Variance Decomposition

```@docs
VarianceDecomposition
```

### Covariance Adjustment Computation

```@docs
MixedModelsSmallSample.compute_adjusted_covariance
```

### Degrees of Freedom Computation

```@docs
MixedModelsSmallSample.compute_dof
```
