# Sainsc.jl

`Sainsc.jl` - Segmentation-free analysis of in situ capture data[^1] - is,
as the name states, a tool to analyse in situ capture-based spatial
transcriptomics data. Actually, it can also be used to analyse imaging-based datasets.
For a more thorough background please refer to the original publication.
or just follow along the provided examplary analysis.

A more _batteries included_[^2] version of this package is available in Python/Rust.
If you are interested have a look at [https://github.com/HiDiHlabs/sainsc](https://github.com/HiDiHlabs/sainsc).

## Citation

If you are using `Sainsc.jl` for your research please consider citing

Müller-Bötticher, N., Tiesmeyer, S., Eils, R., and Ishaque, N.
"Sainsc: a computational tool for segmentation-free analysis of in-situ capture"
bioRxiv (2024) [https://doi.org/10.1101/2024.08.02.603879](https://doi.org/10.1101/2024.08.02.603879)

## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the `Pkg` REPL mode and run:

```
pkg> add Sainsc
```

Or, alternatively, via the `Pkg` API:

```julia
using Pkg
Pkg.add("Sainsc")
```

## Index

```@contents
Pages = [
    "examples/ExampleAnalysis.md",
    "api.md",
]
Depth = 3
```

[^1]: Some people claim it actually stands for _Stupid Acronyms in Science_
[^2]: Pun intended
