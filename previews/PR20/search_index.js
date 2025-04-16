var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = KenwardRoger","category":"page"},{"location":"#KenwardRoger","page":"Home","title":"KenwardRoger","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The confidence intervals and hypothesis tests for the fixed effects in MixedModels.jl are based on large sample approximations. This package implements small sample corrections for these intervals and tests, as described in:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Kenward, Michael G., and James H. Roger. \"Small sample inference for fixed effects from restricted maximum likelihood.\" Biometrics (1997): 983-997.","category":"page"},{"location":"#Getting-Started","page":"Home","title":"Getting Started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using CSV\nusing DataFrames\nusing MixedModels\n\nusing KenwardRoger\n\ndf = DataFrame(MixedModels.dataset(:sleepstudy))\nfm = @formula(reaction ~ 1 + days + zerocorr(1 + days | subj)) # zerocorr is necessary\nm = fit(MixedModel, fm, df; REML=true)\nkr = kenwardroger_matrices(m; FIM_σ²=:expected) # :expected or :observed\nestimates = kenwardroger_estimates(m, kr)","category":"page"},{"location":"#Mathematical-details","page":"Home","title":"Mathematical details","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To quantify the uncertainty on the estimates of the variance components, either the observed or the expected Fisher information matrix can be used, by passing :expected or :observed to the keyword argument FIM_σ². Observed is the default option.","category":"page"},{"location":"#Limitations","page":"Home","title":"Limitations","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Currently, correlated random effects are not supported. In the above example the zerocorr is necessary.","category":"page"},{"location":"#More-examples","page":"Home","title":"More examples","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Blocked experiment with random intercept\nSplit plot experiment\nStrip plot experiment","category":"page"},{"location":"#Similar-software","page":"Home","title":"Similar software","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The R package pbkrtest\nJMP\nSAS","category":"page"},{"location":"#API","page":"Home","title":"API","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modules = [KenwardRoger]","category":"page"}]
}
