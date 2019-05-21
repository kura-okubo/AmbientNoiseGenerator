using PlotlyJS, Random, Distributions

randomseed=100
Random.seed!(randomseed)
rng = Normal(log10(0.1), abs(0.2 * log10(0.1)))

fMrand = 10.0.^(rand(rng, 1000)[:])

trace1 = histogram(x=log10.(fMrand))
layout = Layout(width=500,height=500)
plot(trace1, layout)

