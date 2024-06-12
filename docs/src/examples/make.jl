using Literate

ENV["GKSwstype"] = "100"

Literate.markdown("ExampleAnalysis.jl"; credit=false, execute=true, documenter=true)
