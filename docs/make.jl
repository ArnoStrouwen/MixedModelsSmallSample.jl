using MixedModelsSmallSample
using Documenter

DocMeta.setdocmeta!(
    MixedModelsSmallSample, :DocTestSetup, :(using MixedModelsSmallSample); recursive=true
)

makedocs(;
    modules=[MixedModelsSmallSample],
    authors="Arno Strouwen",
    sitename="MixedModelsSmallSample.jl",
    format=Documenter.HTML(;
        canonical="https://ArnoStrouwen.github.io/MixedModelsSmallSample.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
)

deploydocs(;
    repo="github.com/ArnoStrouwen/MixedModelsSmallSample.jl",
    devbranch="master",
    push_preview=true,
)
