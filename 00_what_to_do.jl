include("main.jl")

useDomination = true
params = ExpeParam(useDomination, findAlternativeSolutions = true)
p = 12
maximalTime = 120 # Maximal time in seconds

instancePath = "./data/tsp/fl3795_julia.txt"

println("========== Solve without rounding the distances")
instance = Instance(instancePath)
instance.p = p
resultsC = solveByClusters(instance, params = params, time_limit=maximalTime)

println("\n\n========== Solve with rounded distances")
instance = Instance(instancePath) # The instance must be reinitialized before a new resolution
instance.p = p
resultsM = solveByModuloClusters(instance, params = params, time_limit=maximalTime)

println("\n\n========== RESULTS ")
println("=== Without rounding: ")
println("\tTime: ", round(resultsC["resolutionTime"], digits = 2), "s")
println("\tObjective value: ", resultsC["radius"])
println("\tBest lower bound: ", resultsC["dualBound"])
println("\tOpened sites: ", resultsC["openedSites"])

println("\n=== With rounding: ")
println("\tTime: ", round(resultsM["resolutionTime"], digits = 2), "s")
println("\tObjective value: ", resultsM["radius"])
println("\tBest lower bound: ", resultsM["dualBound"])
println("\tOpened sites: ", resultsM["openedSites"])


