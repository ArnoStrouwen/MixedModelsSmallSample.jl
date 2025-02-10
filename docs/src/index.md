```@meta
CurrentModule = KenwardRoger
```

# KenwardRoger

The confidence intervals and hypothesis tests for the fixed effects in
[MixedModels.jl](https://juliastats.org/MixedModels.jl/stable/)
are based on large sample approximations.
This package implements small sample corrections for these intervals and tests, as described in:

> Kenward, Michael G., and James H. Roger. "Small sample inference for fixed effects from restricted maximum likelihood." Biometrics (1997): 983-997.

## Getting Started

```@example
using CSV
using DataFrames
using MixedModels

using KenwardRoger

df = DataFrame(MixedModels.dataset(:sleepstudy))
fm = @formula(reaction ~ 1 + days + zerocorr(1 + days | subj)) # zerocorr is necessary
m = fit(MixedModel, fm, df; REML=true)
kr = kenwardroger_matrices(m; FIM_σ²=:expected) # :expected or :observed
estimates = kenwardroger_estimates(m, kr)
```

## Mathematical details

To quantify the uncertainty on the estimates of the variance components,
either the observed or the expected Fisher information matrix can be used,
by passing `:expected` or `:observed` to the keyword argument `FIM_σ²`.
Observed is the default option.

## Limitations

Currently, correlated random effects are not supported.
In the above example the `zerocorr` is necessary.

## More examples

  - [Blocked experiment with random intercept](https://github.com/ArnoStrouwen/KenwardRoger.jl/blob/master/test/blocked%20experiment.jl)
  - [Split plot experiment](https://github.com/ArnoStrouwen/KenwardRoger.jl/blob/master/test/split%20plot%20experiment.jl)
  - [Strip plot experiment](https://github.com/ArnoStrouwen/KenwardRoger.jl/blob/master/test/strip%20plot%20experiment.jl)

## Similar software

  - The R package [pbkrtest](https://github.com/hojsgaard/pbkrtest)
  - [JMP](https://www.jmp.com/support/help/en/18.1/index.shtml#page/jmp/statistical-details-for-the-kackarharville-correction-2.shtml#)
  - [SAS](https://documentation.sas.com/doc/en/statcdc/14.2/statug/statug_glimmix_details40.htm)

## API

```@autodocs
Modules = [KenwardRoger]
```
