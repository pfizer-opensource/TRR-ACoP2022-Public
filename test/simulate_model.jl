# Step 1: Load the model as @reaction_network from Catalyst with inital values and parameter values from CSV file.
using Catalyst, CSV, DataFrames
veryhitgcsv = "./log/veryhitg.csv"
csvfile = veryhitgcsv
df = CSV.read(csvfile,DataFrame)
include("../src/rxns.jl")
rn = get_rxns(df)

# Step 2: Convert the model to an ODESystem and simplify it.
rn_odesys = structural_simplify(convert(ODESystem,rn;combinatoric_ratelaws=false))

# Step 3: Solve ODE.
using DiffEqCallbacks, OrdinaryDiffEq
tspan = (0.0,10000.0) # [days]
cb = TerminateSteadyState(1e-8,1e-6) # callback to stop at steadystate, more efficient for M-H search
op = ODEProblem(rn_odesys, [], tspan, [], callback=cb)
sol = solve(op, Rodas5())

# Step 4: Visualize solution.
using Plots
p = plot(sol, legend = :outerright)
