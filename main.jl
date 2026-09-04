using JuMP                  
using CPLEX
using Clustering
using Random
using StatsBase
using Dates

global EPS = 0.0001

Random.seed!(1230);

include("expeParam.jl")
include("instance.jl")
include("localSearch.jl")

"""
Solve a p-center instance by only considering a subset of clients and adding clients as needed.
The clients are clustered and each cluster is divided in quadrants. When clients must be added, we try to find in each quadrant of each cluster one client that would increase the radius of the current solution.
The relaxation of the p-center problem is solved as long as it enables to add clients which would improve the radius. When there are no more, the integer problem is solved once (and then we go back to solving the relaxation).

Input:
- instance: the instance (in which the clients may already include the partition of the clients
- params: structure which define some additional solution parameters
- modulo: all distances computed will be rounded down to 10^(modulo-1)
- initClusters: true if the client clusters must be computed (must be set to true if the problem is solved only using this function)
- initialTimeElapsed: time elapsed in the resolution before calling this function 
- time_limit: remaining time
- useImprovedRounding: true if a fractional solution is rounded by opening the sites by considering them in decreasing order of their y_j variables and if a site is opened only if it improves the radius (without improvement, a site is opened even if it does not improve the radius)
- skipDandOfInit: true if the distances must not be recomputed (must be set to false if the problem is solved only using this function)
- isRelaxation: true if the linear relaxation of the problem is solved
- improveModuloUB: true if this method is called by the solveByModuloClusters method and if its bound is updated within solveByClusters
- computeAlternativeRoundedSolution: true if two integer solutions are computed from a fractional solution (one with the improved rounding, one without) and the best is kept 
"""
function solveByClusters(instance::Instance; params::ExpeParam=ExpeParam(false), modulo::Int=1, initClusters::Bool=true, initialTimeElapsed::Float64=0.0, time_limit::Int=-1, useImprovedRounding::Bool=false, skipDandOfInit::Bool=false, isRelaxation::Bool=false, improveModuloUB::Bool=false, computeAlternativeRoundedSolution::Bool=false)

    startingTime = time()
    
    if initClusters
        createClusters!(instance, params)
        instance.computedDistances = 0
    end

    dominationTime = 0
    if !skipDandOfInit
        initStartingTime = time()
        initializeDandOf!(instance, params, modulo=modulo, time_limit=time_limit)
    end

    isClusterUpdated = true
    
    # True only when the loop stopped because no client could be removed from the clusters any
    # more, i.e. when the radius found on the cluster representatives is proved to be the radius
    # of every client.
    # Used to avoid considering that the result is optimal when the reason to stop is that there is no time left
    hasConverged = false

    radius = nothing
    dichotomyResults = nothing
    dichotomyTime = 0
    clusterUpdateTime = 0
    alternativeSolutionsTime = 0
    
    iterationCount = 0
    radiusOfRepresentatives = nothing

    while isClusterUpdated && !isOver(time_limit, startingTime) # While the number of clusters is increased

        iterationCount += 1

        while isClusterUpdated && !isOver(time_limit, startingTime) && instance.lb < instance.ub
            println("\n-- Continuous resolution (elapsed time ", round(Int, time() - startingTime + initialTimeElapsed), "s)")

            relaxationStartingTime = time()
            ## Relaxation resolution
            # Solve p-center problem associated to the clusters in the instance
            print("Solving with ", length(instance.clusters), " clusters")
            
            if remainingTime(time_limit, startingTime) != -1
                print(", max time: ", remainingTime(time_limit, startingTime))
            end
            print("... ")

            dichotomyResults = solveByDichotomy(instance, isInitialized=true, isRelaxation=true, modulo=modulo, time_limit=remainingTime(time_limit, startingTime), params = params, useImprovedRounding=useImprovedRounding, returnFractionalSolution=isRelaxation, computeAlternativeRoundedSolution=computeAlternativeRoundedSolution)
            dichotomyTime += dichotomyResults["resolutionTime"]
            
            println("... done in ", round(Int, time()-relaxationStartingTime), "s. Radius found: ", dichotomyResults["radius"])

            if dichotomyResults["isOptimal"]
                
                if dichotomyResults["optimalValueFound"]
		    instance.ub = dichotomyResults["radius"]
		    instance.lb = instance.ub
                    isClusterUpdated = false
                else                 
                    lbImproved = instance.lb < dichotomyResults["radius"]
                    if lbImproved
                        println("Set LB to ", dichotomyResults["radius"])
                        setLB!(instance, dichotomyResults["radius"], params)
                    end

                    ## Compute the actual radius of the opened sites
                    radiusOfRepresentatives = typemin(Int64)

                    # Get the distance of each cluster representative to its closest opened site
                    dichotomyOpenedSites = dichotomyResults["openedSites"]
                    for clusterId in 1:length(instance.clusters)
                        cDist = distanceToOpenSites(instance, clusterId, dichotomyOpenedSites, modulo=modulo, stoppingValue = radiusOfRepresentatives, isClusterId = true)

                        if cDist > radiusOfRepresentatives
                            radiusOfRepresentatives = cDist
                        end
                    end

                    openedSitesSolutions = Vector{Vector{Int}}([dichotomyResults["openedSites"]])

                    # If a second solution is obtained from the fractional one at the end of the continuous dichotomy
                    alternativeSolution = dichotomyResults["alternativeRoundedSolution"]
                    if computeAlternativeRoundedSolution && alternativeSolution != nothing
                        push!(openedSitesSolutions, dichotomyResults["alternativeRoundedSolution"])
                    end

                    if params.findAlternativeSolutions
                        alternativeSolutionsStartingTime = time()
                        append!(openedSitesSolutions, findAlternativeOpenedSites(instance, dichotomyResults["openedSites"], radiusOfRepresentatives, modulo, time_limit=remainingTime(time_limit, startingTime)))
                        alternativeSolutionsTime += time() - alternativeSolutionsStartingTime
                    end 

                    ## Update the clusters if necessary
                    clusterUpdateStartingTime = time()

                    if isRelaxation
                        isClusterUpdated, domTime = updateClusters!(instance, radiusOfRepresentatives, openedSitesSolutions, lbImproved, isIntResolution = false, params = params, modulo=modulo, fractionalSitesOpened=dichotomyResults["fractionalSolution"], time_limit=remainingTime(time_limit, startingTime), improveModuloUB=improveModuloUB)
                        dominationTime += domTime
                    else 
                        isClusterUpdated, domTime = updateClusters!(instance, radiusOfRepresentatives, openedSitesSolutions, lbImproved, isIntResolution = false, params = params, modulo=modulo, time_limit=remainingTime(time_limit, startingTime), improveModuloUB=improveModuloUB)
                        dominationTime += domTime
                    end

                    clusterUpdateTime += time() - clusterUpdateStartingTime
                    
                end # if dichotomyResults["optimalValueFound"] (i.e., if the optimal solution is not found)
                
                if modulo == 1
                    println("Bounds: [", instance.lb, ", ", instance.ub, "], instance.of size: ", length(instance.of))
                else
                    displayedUB = string(instance.moduloUB) * "]"

                    if instance.moduloUB == typemax(Int)
                        if instance.ub == typemax(Int)
                            displayedUB = "+oo["
                        else
                            displayedUB = string(instance.ub) * "]"
                        end
                    else
                        if instance.ub == typemax(Int)
                            displayedUB = string(instance.moduloUB) * "]"
                        else 
                            displayedUB = string(min(instance.moduloUB, instance.ub  + modulo)) * "]"
                        end 
                    end 
                    
                    println("Bounds: [", instance.lb, ", ", displayedUB, ", number of distances considered: ", length(instance.of))
                end 
            end # if dichotomyResults["isOptimal"]
            if !dichotomyResults["isOptimal"] # If the relaxation could not be solved (time limit, or a set cover without a conclusion). 
                # Avoid solving it again 
                break
            end
        end # while isClusterUpdated && !isOver(time_limit, startingTime) 

        isClusterUpdated = false
        
         ## Integer resolution
        if !isRelaxation && !isOver(time_limit, startingTime) && instance.lb >= instance.ub

            # The bounds already met: the radius is proved, there is nothing left to solve
            hasConverged = true
        end

        if !isRelaxation && !isOver(time_limit, startingTime) && instance.lb < instance.ub
            println("\n-- Integer resolution (elapsed time ", round(Int, time() - startingTime + initialTimeElapsed), "s)")  

            intStartingTime = time()

            # Solve p-center problem associated to the clusters in the instance
            print("Solving with ", length(instance.clusters), " clusters")
            
            if remainingTime(time_limit, startingTime) != -1
                print(", max time: ", remainingTime(time_limit, startingTime))
            end
            print("... ")
            
            dichotomyResults = solveByDichotomy(instance, isInitialized=true, modulo=modulo, time_limit=remainingTime(time_limit, startingTime), params = params)
            dichotomyTime += dichotomyResults["resolutionTime"]  
            
            radiusOfRepresentatives = dichotomyResults["radius"]
            println("... done in ", round(Int, time()-intStartingTime), "s.")

            if dichotomyResults["isOptimal"]
                if dichotomyResults["optimalValueFound"] 
                    isClusterUpdated = false
                    hasConverged = true
		    instance.lb = radiusOfRepresentatives
		    instance.ub = radiusOfRepresentatives
                else
                lbImproved = instance.lb < radiusOfRepresentatives
                if lbImproved
                    println("Set LB to ", radiusOfRepresentatives)
                    setLB!(instance, radiusOfRepresentatives, params)
                end

                    openedSitesSolutions = Vector{Vector{Int}}([dichotomyResults["openedSites"]])
                    
                    if params.findAlternativeSolutions
                        alternativeSolutionsStartingTime = time()
                        append!(openedSitesSolutions, findAlternativeOpenedSites(instance, dichotomyResults["openedSites"], radiusOfRepresentatives, modulo, time_limit=remainingTime(time_limit, startingTime)))
                        alternativeSolutionsTime += time() - alternativeSolutionsStartingTime
                    end
                    
                    # Update the clusters if necessary
                    clusterUpdateStartingTime = time()
                    
		    try
			isClusterUpdated, domTime = updateClusters!(instance, radiusOfRepresentatives, openedSitesSolutions, lbImproved, isIntResolution = true, params = params, modulo=modulo, time_limit=remainingTime(time_limit, startingTime), improveModuloUB=improveModuloUB)
                        dominationTime += domTime
                        hasConverged = !isClusterUpdated
		    catch e
		        println("Error: ")
			rethrow(e)
		    end

                    clusterUpdateTime += time() - clusterUpdateStartingTime
                end 
                    
                if modulo == 1
                    println("Bounds: [", instance.lb, ", ", instance.ub, "], number of distances considered: ", length(instance.of))
                else
                    displayedUB = string(instance.moduloUB) * "]"

                    if instance.moduloUB == typemax(Int)
                        if instance.ub == typemax(Int)
                            displayedUB = "+oo["
                        else
                            displayedUB = string(instance.ub) * "]"
                        end
                    else
                        if instance.ub == typemax(Int)
                            displayedUB = string(instance.moduloUB) * "]"
                        else 
                            displayedUB = string(min(instance.moduloUB, instance.ub  + modulo)) * "]"
                        end 
                    end 
                    
                    println("Bounds: [", instance.lb, ", ", displayedUB, ", number of distances considered: ", length(instance.of))
                end 
            end
        end 
    end
    totalTime= round(Int, time() - startingTime + initialTimeElapsed)
    results = Dict{String, Any}()

    results["dichotomyTime"] = dichotomyTime
    results["clusterUpdateTime"]= clusterUpdateTime
    results["isOptimal"] = !isRelaxation && !isOver(time_limit, startingTime) && hasConverged
    results["resolutionTime"] = time() - startingTime
    results["alternativeSolutionTime"] = alternativeSolutionsTime
    results["dominationTime"] = dominationTime
    results["n"] = instance.n
    results["p"] = instance.p
    results["radius"] = instance.ub

    if !isRelaxation && !hasConverged && !isOver(time_limit, startingTime)
        println("Warning: the cluster loop stopped before converging although time remained ",
                "(a set cover could not be solved to a conclusion). The radius is NOT proved optimal.")
    end

    if results["isOptimal"] && !isRelaxation && modulo == 1
        results["dualBound"] = results["radius"]
    else
        results["dualBound"] = instance.lb
    end 
    
    results["openedSites"] = instance.ubOpenedSites 
    results["clientsAtTheEnd"] = instance.dRows
    
    return results
end 

"""
Solve a p-center problem using dichotomy

Input
- instance: the instance
- params: additional parameters
- isInitialized: true if the clusters and the distances are already computed
- isRelaxation: true if we only solve the linear relaxation of the problem
- modulo: all distances computed will be rounded down to 10^(modulo-1)
- time_limit: remaining time
- useImprovedRounding: true if a fractional solution is rounded by opening the sites by considering them in decreasing order of their y_j variables and if a site is opened only if it improves the radius (without improvement, a site is opened even if it does
- returnFractionalSolution: true if the fractional solution is returned
- computeAlternativeRoundedSolution: true if two integer solutions are computed from a fractional solution (one with the improved rounding, one without) and the best is kept 
"""
function solveByDichotomy(instance::Instance; params::ExpeParam=ExpeParam(false), isInitialized::Bool=false, isRelaxation::Bool=false, modulo::Int=1, time_limit::Int=-1, useImprovedRounding::Bool=false, returnFractionalSolution::Bool=false, computeAlternativeRoundedSolution::Bool=false)

    startingTime = time()

    if !isInitialized
        
        println("Initializing distances...")
        instance.lb = typemin(Int)
        instance.ub = typemax(Int)
        createUnitaryClusters(instance)
        initializeDandOf!(instance, params, time_limit=time_limit)
        println(" done in ", round(time()-startingTime, digits = 2), "s")
    end
    initializationTime = time()-startingTime
    
    lbId = 1
    ubId = length(instance.of)

    lb = instance.of[lbId]
    ub = instance.of[ubId]
    x = nothing
    isFeasible = nothing
    optimalValueFound = false # True if the optimal value of the problem for ALL clients is found

    isFirstStep = true
    t = 0

    # True if a set cover could not be solved to a conclusion, in which case the dichotomy is
    # stopped and its bounds cannot be claimed to be proved
    isDichotomyAborted = false

    lastFeasibleOpenedSites = instance.ubOpenedSites
    alternativeOpenedSites = nothing # Alternative integer solution obtained from a fractional one (use if computeAlternativeRoundedSolution is true)

    # While:
    # - the bounds are not equal;
    # - the time is not over; and
    # - the UB is not proved to be optimal
    while lbId < ubId && !isOver(time_limit, startingTime)  && !optimalValueFound
            
        if isFirstStep
            testedId = lbId
            isFirstStep = false
        else
            testedId = floor(Int, (lbId+ubId)/2)
        end

        # The optimal value for the problem with all clients is found if:
	# - there is only two possible indexes (lbId and ubId); and
	# - of[lbId] + modulo > ub (in that case, of[lbId] is the optimal value since there is no value rounded down to a multiple of modulo between of[lbId] and ub)
        if lbId+1 == ubId && instance.of[lbId]+modulo > instance.ub
            optimalValueFound = true
            isFeasible = true
        else
            isFeasible,  x, isConclusive= areClientsCoverable(instance, instance.of[testedId], isRelaxation=isRelaxation, time_limit=remainingTime(time_limit, startingTime), params = params)

            if !isConclusive

                # Leave the bounds untouched and stop: treating this as infeasible would raise lbId
                # without proof
                println("Set cover not solved to a conclusion, stopping the dichotomy")
                isDichotomyAborted = true
                break
            end

            if isFeasible
                ubId = testedId
                print("l")
                lastFeasibleOpenedSites = Vector{Int}([])

                if !isRelaxation
                    for j in 1:instance.m
                        if instance.siteDomination[j] == 0 && JuMP.value(x[j]) > 0.9
                            push!(lastFeasibleOpenedSites, j)
                        end 
                    end
                end 
            else
                print("r")
                lbId = testedId + 1
            end
        end             
    end

    radius = instance.of[lbId]
    
    # If the last set cover was infeasible (i.e., if the lower bound was increased)
    if !isDichotomyAborted && !isOver(time_limit, startingTime) && (length(instance.of) == 1 || !isFeasible)
        # Solve the problem for the lower bound
        isFeasible, x, isConclusive = areClientsCoverable(instance, instance.of[lbId], isRelaxation = isRelaxation, time_limit=remainingTime(time_limit, startingTime), params = params)
        if !isConclusive
            isDichotomyAborted = true
        end 
    end

    isOptimal = false

    results = Dict{String, Any}()
    
    if !isDichotomyAborted && !isOver(time_limit, startingTime)

        if time_limit == -1 || remainingTime(time_limit, startingTime) > 5
            isOptimal = true
        end 

        if isRelaxation && !optimalValueFound
            if isOptimal
                
                lastFeasibleOpenedSites = Vector{Int64}([]) 
                sites = collect(1:instance.m)
                xValues = Vector{Float64}(zeros(instance.m))

                for (j, jDomination) in enumerate(instance.siteDomination)
                    if jDomination == 0
                        xValues[j] = JuMP.value(x[j])
                    end
                end 
                
                # Order the site by decreasing values of x
                permutations = sortperm(xValues, rev = true)
                orderedSites = sites[permutations]

                if returnFractionalSolution
                    minimalFractionalValue = 1E-5
                    fractionalSites = Vector{Tuple{Int, Float64}}(undef, instance.m)#maximalNumberOfSites)

                    nextSiteId = 1
                    nextJ = permutations[nextSiteId]
                    nextSiteValue = xValues[nextJ]

                    while nextSiteId <= instance.m && nextSiteValue >= minimalFractionalValue 
                        fractionalSites[nextSiteId] = Tuple{Int, Float64}((nextJ, nextSiteValue))

                        nextSiteId += 1

                        if nextSiteId <= instance.m
                            nextJ = permutations[nextSiteId]
                            nextSiteValue = xValues[nextJ]
                        end 
                    end
                    results["fractionalSolution"] = fractionalSites[1:nextSiteId-1]
                end 

                if !useImprovedRounding || computeAlternativeRoundedSolution
                    lastFeasibleOpenedSites = @view orderedSites[1:min(instance.p, length(orderedSites))] # Open the sites which have the largest value of y variables

                    if computeAlternativeRoundedSolution
                        alternativeOpenedSites, testRadius = getIntegerSolution(instance, orderedSites, xValues, permutations)

                        if alternativeOpenedSites == lastFeasibleOpenedSites
                            alternativeOpenedSites = nothing
                            println("Identical alternative solution found")
                        end 
                    end 
                else 
                    lastFeasibleOpenedSites, testRadius =  getIntegerSolution(instance, orderedSites, xValues, permutations)
                end
            end 
        end
        results["radius"] = radius 
        results["dualBound"] = radius
    else
        results["radius"] = min(instance.ub, instance.of[ubId]) # instance.of[ubId] can be > instance.ub when ubId is not updated
        results["dualBound"] = instance.of[lbId]
    end

    if optimalValueFound 
        openedSites = instance.ubOpenedSites
    else
        openedSites = lastFeasibleOpenedSites
    end

    results["initializationTime"] = initializationTime
    results["openedSites"] = openedSites
    results["resolutionTime"] = time() - startingTime
    results["isOptimal"] = isOptimal
    results["n"] = instance.n
    results["optimalValueFound"] = optimalValueFound
    results["alternativeRoundedSolution"] = alternativeOpenedSites

    return results
end 

"""
Solve a set cover to determine if the clients can be covered within a given distance

Input
- instance: the instance
- distance: the distance within which we test if the clients can be covered
- isRelaxation: true if we only solve the linear relaxation of the problem
- time_limit: remaining time
- params: additional parameters
"""
function areClientsCoverable(instance, distance::Int; isRelaxation::Bool=false, time_limit::Int=-1, params::ExpeParam=ExpeParam(false))

    m = Model(CPLEX.Optimizer)
    set_silent(m)

    if isRelaxation
        @variable(m, 0 <= x[j in 1:instance.m;  instance.siteDomination[j] == 0] <= 1)
    else
        @variable(m, x[j in 1:instance.m;  instance.siteDomination[j] == 0], Bin)
    end
    
    if time_limit <= 5 && time_limit != -1
        return false, x, false
    else
        if time_limit != -1
            set_optimizer_attribute(m, "CPX_PARAM_TILIM", time_limit)
        end 
    end 

    if params.useNullObjective
        @objective(m, Min, 0)
    else
        @objective(m, Min, sum(x[j] for j in 1:instance.m if instance.siteDomination[j] == 0))
    end

    if params.useCutOff && !isRelaxation
        set_optimizer_attribute(m, "CPXPARAM_MIP_Tolerances_UpperCutoff", instance.p)
        set_optimizer_attribute(m, "CPX_PARAM_INTSOLLIM", 1)
    end 

    @constraint(m, client[i in 1:instance.dRows], sum(x[j] for j in 1:instance.m if instance.siteDomination[j] == 0 && instance.d[i, j] <= distance) >= 1)

    if params.usePMaxConstraint || params.useNullObjective
        @constraint(m, sum(x[j] for j in 1:instance.m if instance.siteDomination[j] == 0) <= instance.p)
    end 

    isCoverable = false

    # True if the answer is a proof. False if the set cover was interrupted (e.g., time limit)
    # Reporting it as "not coverable" would raise the lower bound of the
    # dichotomy without any proof
    isConclusive = true

    try
        optimize!(m)
    catch e
        println("error while solving the set cover: ", sprint(showerror, e))
        return false, x, false
    end    

    status = termination_status(m)

     if primal_status(m) == MOI.FEASIBLE_POINT
         isCoverable = JuMP.objective_value(m) <= instance.p + 10^-3

        # A solution using more than p sites only proves that the clients are not coverable if the
        # search finished; if CPLEX stopped early, a cover with <= p sites may still exist
        if !isCoverable && status != MOI.OPTIMAL && status != MOI.OBJECTIVE_LIMIT
            isConclusive = false
        end
    else

        # No solution found: only CPLEX proving infeasibility (or cutting everything off with the
        # cutoff set to p) shows that the clients are not coverable
        isConclusive = status == MOI.INFEASIBLE || status == MOI.INFEASIBLE_OR_UNBOUNDED || status == MOI.OBJECTIVE_LIMIT
    end

    return isCoverable, x, isConclusive

end  

"""
Solve a p-center instance by clustering the clients. The distances are rounded down at each iteration providing upper and lower bounds. At each iteartion the rounded p-center problem is solved using the solveByCluster method.
 At the next iteration the rounding is more precise by one digit until the optimal solution is reached.

Input
- instance: the instance (in which the clients may already include the partition of the clients
- params: structure which define some additional solution parameters
- time_limit: remaining time
- useImprovedRounding: true if a fractional solution is rounded by opening the sites by considering them in decreasing order of their y_j variables and if a site is opened only if it improves the radius (without improvement, a site is opened even if it does
- isRelaxation: true if the linear relaxation of the problem is solved
- useKmeansA: compute a k-means solution with p clusters and use it to create an initial solution
- improveModuloUB: true if at each iteration, the global upper bound is updated within the solveByCluster method
- base: base used to round the distances (i.e., the distances will be rounded down to their closest base^i value with i decreasing during the iterations)
- computeAlternativeRoundedSolution: true if two integer solutions are computed from a fractional solution (one with the improved rounding, one without) and the best is kept
- limitIterationsTime: true if the duration of each iteration (except the last one) is limited
"""
function solveByModuloClusters(instance::Instance; params::ExpeParam=ExpeParam(false), time_limit::Int=-1, useImprovedRounding::Bool=false, isRelaxation::Bool=false, useKmeansA::Bool=false, improveModuloUB::Bool=false, base::Int=10, computeAlternativeRoundedSolution::Bool=false, limitIterationsTime::Bool=false)

    maxModuloExponent = 1

    results = Dict{String, Any}()
    results["dichotomyTime"] = 0
    results["dominationTime"] = 0
    results["alternativeSolutionTime"] = 0
    results["clusterUpdateTime"]= 0
    
    previousLB = typemin(Int)
    previousUB = typemax(Int)

    # Initialize the clusters        
    println("Creating the clusters...")
    createClusters!(instance, params, useKmeansA=useKmeansA)

    if useKmeansA
        results["KMAUB"] = instance.ub
    end

    instance.computedDistances = 0

    # Get the average distance between the cluster representatives and the sites to set the first modulo value
    averageDistance = Float64(0.0)
    for (clusterId, cluster) in enumerate(instance.clusters)
        for siteId in 1:instance.m
            averageDistance += distance(instance, cluster.representativeId, siteId)
        end
    end

    averageDistance /= (length(instance.clusters)*instance.m)
    maxModuloExponent = num_digits(averageDistance, base) -1

    startingTime = time()
    instance.moduloOpenedSites = Vector{Int}([])
    
    moduloExp = maxModuloExponent
    clusterResults = nothing

    resolutionTimeByIterations = Vector{Float64}([])
    
    while moduloExp >= 0 && !isOver(time_limit, startingTime) && instance.lb < instance.ub

        iterationStartingTime = time()

        # The maximal solution time for this iteration is the remaining by default
        # If the time is limited only allow 25% of this time (but no less than 5 minutes)
        maxIterationTime = remainingTime(time_limit, startingTime)
        if limitIterationsTime && moduloExp > 0
            maxIterationTime = max(round(Int, maxIterationTime/4), 300)
        end

        println("\n\n=== Solving with modulo ", moduloExp, " (elapsed time: ", round(Int, time() - startingTime), "s, bounds [", instance.lb, ", ", instance.ub, "])")
        clusterResults = solveByClusters(instance, params=params, modulo=base^moduloExp, initClusters=false, initialTimeElapsed = time() - startingTime, time_limit=maxIterationTime, useImprovedRounding=useImprovedRounding, isRelaxation=isRelaxation, improveModuloUB=improveModuloUB, computeAlternativeRoundedSolution=computeAlternativeRoundedSolution)
        results["dichotomyTime"] += clusterResults["dichotomyTime"]
        results["alternativeSolutionTime"] += clusterResults["alternativeSolutionTime"]
        results["dominationTime"] += clusterResults["dominationTime"]
        results["clusterUpdateTime"] += clusterResults["clusterUpdateTime"]
        instance.ub = instance.moduloUB
        
        # Set the bounds accordingly
        if isRelaxation || moduloExp > 0 || !clusterResults["isOptimal"]
            instance.lb = clusterResults["dualBound"]
        else
            instance.lb = clusterResults["radius"]
        end

        boundOfOpenedSites = typemin(Int)

	boundStartingTime = time()
	boundComputed = true
        if length(clusterResults["openedSites"]) > 0
            for clientId in 1:size(instance.clientsCoordinates, 2)
                clientDistance = distanceToOpenSites(instance, clientId, clusterResults["openedSites"], modulo=1, stoppingValue = boundOfOpenedSites)

                if clientDistance > boundOfOpenedSites
                    boundOfOpenedSites = clientDistance
                end 
		if rem(clientId, 100000) == 1 && isOver(time_limit, startingTime)
		    boundComputed = false
		    break
                end
            end 
        else 
            boundOfOpenedSites = typemax(Int)   
        end

        if boundComputed
            if boundOfOpenedSites < instance.moduloUB
                instance.moduloUB = boundOfOpenedSites
                instance.moduloOpenedSites = clusterResults["openedSites"]
            end
        end
        instance.ub = instance.moduloUB

        if !isOver(time_limit, startingTime)
            computeLB(instance, params=params, useImprovedRounding=useImprovedRounding)
        end 

        moduloExp -= 1  

        push!(resolutionTimeByIterations, time() - iterationStartingTime)
    end 

    results["isOptimal"] = instance.lb >= instance.ub
    results["radius"] = instance.moduloUB
    results["resolutionTime"] = time()-startingTime
    results["resolutionTimeByIterations"] = resolutionTimeByIterations
    results["openedSites"] = instance.moduloOpenedSites
    results["n"] = instance.n
    results["clientsAtTheEnd"] = instance.dRows
    results["dualBound"] = instance.lb
    
    return results
end

"""
Returns the remaining time
"""
function remainingTime(time_limit::Int, startingTime)
    if time_limit == -1
        return -1
    else
        return round(Int, time_limit - (time() - startingTime))
    end 
end 

"""
Test if the time is over
"""
function isOver(time_limit::Int, startingTime)
    return time_limit != -1 && time() - startingTime > time_limit - 5
end 

"""
Cluster the clients of the instance

Input
- instance: the instance (in which the clients may already include the partition of the clients
- params: structure which define some additional solution parameters
- useKmeansA: compute a k-means solution with p clusters and use it to create an initial solution
"""
function createClusters!(instance::Instance, params::ExpeParam; useKmeansA::Bool=false)
    
    # Cluster the clients
    startingTime = time()

    instance.lb = typemin(Int)
    instance.ub = typemax(Int)

    # If we compute an initial solution of the p-center with the adapted k-means
    if useKmeansA
        result = kmeans(instance.clientsCoordinates, instance.p) # Compute the clusters and their barycenters
        instance.ub = typemin(Int)
        instance.ubOpenedSites = Vector{Int64}(undef, instance.p)

        # For each cluster
        for clusterId in 1:instance.p
            c = Cluster()
            c.clientsId = findall(result.assignments .== clusterId)
            c.representative = Vector{Float64}(zeros(2))

            updateRepresentative!(c, instance) # Get its representative id

            instance.ubOpenedSites[clusterId] = c.representativeId

            repCoordinates = @view instance.clientsCoordinates[:, c.representativeId]

            # Update UB if necessary
            for (clientClusterId, clientId) in enumerate(c.clientsId)
                d = distance((@view instance.clientsCoordinates[:, clientId]), repCoordinates)

                if d > instance.ub
                    instance.ub = d
                end
            end 
        end
    end 
    
    k = min(instance.n, instance.p + params.clusterCountIncrease)
    instance.initialClusterCount = k

    result = kmeans(instance.clientsCoordinates, k)
    kmeansTime = round(Int, time() - startingTime)

    # Create the clusters
    instance.clusters = Vector{Cluster}(undef, k)

    # Add the clusters found by k-means
    for clusterId in 1:k
        instance.clusters[clusterId] = Cluster(instance, findall(result.assignments .== clusterId))
    end

    instance.areClientsIdenticalToSites = false
end 


"""
Select p sites to open from a fractional value of the vector y.
A site is opened if it enables to reduce the radius.

Input:
- instance: the instance
- orderedSites: list of sites indexes ordered by decreasing of value of y (dominated sites are considered to have y = 0)
- yValues: value of the y variables
- permutations: permutation of the sites from 1:n to their decreasing order of value of variables y
"""
function getIntegerSolution(instance::Instance, orderedSites::Vector{Int}, yValues::Vector{Float64}, permutations::Vector{Int}; percentageOfRadiusReduction=1)

    # Open the first site
    openedSites = Vector{Int}(undef, instance.p)
    openedSites[1] = orderedSites[1]
    addedSites = 1

    # Compute the radius of each cluster (which is currently the distance to the first site)
    clusterRadius = Vector{Int}(undef, length(instance.clusters))

    # Radius with the current sites
    radius = -1

    # Clusters which allocation distance is equal to the radius
    clustersAtRadius = Vector{Int}()

    # For each cluster
    for clusterId in 1:length(instance.clusters)

         # Its allocation distance is equal to the distance with the unique current site
         cr = instance.d[clusterId, openedSites[1]]
        clusterRadius[clusterId] = cr

        # Update radius and clustersAtRadius if necessary
        if cr >= radius
            if cr == radius
                push!(clustersAtRadius, clusterId)
            else
                clustersAtRadius = Vector{Int}([clusterId])
                radius = cr
            end
        end 
    end
    
    # Minimal percentage in [0, 1] of clusters (among the clusters which allocation distance is maximal) which allocation distance must be reduced to enable the addition of a site. Initially equal to 1, its value will decrease if no suitable site is found.
    percentageOfRadiusReduction = 1.0

    # Number of clients at maximal radius that a site must cover to be added
    cARPercentage1 = percentageOfRadiusReduction * length(clustersAtRadius)

    # Also used to avoid checking unecessary clusters
    cARPercentage2 = length(clustersAtRadius) - cARPercentage1

    # Get the index of the last site that could be considered
    maxSiteIndex = min(10*instance.p, size(instance.sitesCoordinates, 2))

    # Find through binary search the first index of a site which y value is equal to 0
    lb = 1
    ub = maxSiteIndex

    while lb < ub

        testedSiteIndex = ceil(Int, lb+ (ub-lb)/2)

        if yValues[permutations[testedSiteIndex]] <= EPS
            ub = testedSiteIndex - 1
        else
            lb = testedSiteIndex
        end 
    end

    maxSiteIndex = max(instance.p, lb)
    siteTestedId = 2
    
    while addedSites < instance.p

        ## Test if adding the next site reduces the radius
        ## (i.e., if it enables to cover all the  which allocation distance is equal to the radius)
        clusterTestedId = 1
        clusterReduced = 0

        # If the site must cover new clusters at the maximal allocation distance
        if cARPercentage1 > 0

            # While a sufficient number of reduced clusters are not found and can still be found
            # i.e., clusters which allocation distance is reduced + clusters at the radius - id of the currently tested cluster + 1 >= minimal percentage of clusters at the radius
            while clusterReduced < cARPercentage1 && cARPercentage2 + clusterReduced - clusterTestedId + 1 >= 0 && clusterTestedId <= length(clustersAtRadius)
                if instance.d[clustersAtRadius[clusterTestedId], orderedSites[siteTestedId]] < radius
                    clusterReduced += 1
                end 
                clusterTestedId += 1 
            end
        elseif  orderedSites[siteTestedId] in @view openedSites[1:addedSites] # If the site must not cover any cluster but  it not already in the solution do not add it
            clusterReduced = -1
        end  

        # If the site enables to cover a sufficient number of clusters in clusterAtRadius, add it
        if clusterReduced >= cARPercentage1
            addedSites += 1
            openedSites[addedSites] = orderedSites[siteTestedId]

            isRadiusReduced = clusterReduced == length(clustersAtRadius)
            # If the radius is not reduced
            # i.e., if the allocation distance of all the clusters in clustersAtRadius are not reduced
            if !isRadiusReduced

                # Remove from clusterAtRadius the clusters which allocation distance is reduced
                for (clusterAtRadiusId, clusterId) in Iterators.reverse(enumerate(clustersAtRadius))

                    #clusterId = clustersAtRadius[clusterAtRadiusId]
                    if instance.d[clusterId, openedSites[addedSites]] < radius
                        deleteat!(clustersAtRadius, clusterAtRadiusId)
                    end
                end
            else # If the radius is reduced
                
                ## Update radius and clustersAtRadius
                radius = -1
                clustersAtRadius = Vector{Int}()
            end

            # For each cluster
            for (clusterId, cRadius) in enumerate(clusterRadius)

                # If the new site enables to improve its allocation distance
                if instance.d[clusterId, openedSites[addedSites]] < cRadius
                    clusterRadius[clusterId] = instance.d[clusterId, openedSites[addedSites]]

                    # Update radius and clustersAtRadius if necessary
                    if isRadiusReduced && clusterRadius[clusterId] >= radius
                        if clusterRadius[clusterId] == radius
                            push!(clustersAtRadius, clusterId)
                        else
                            clustersAtRadius = Vector{Int}([clusterId])
                            radius = clusterRadius[clusterId]
                        end
                    end 
                end 
            end

            siteTestedId = 2
            percentageOfRadiusReduction = 1.0
            cARPercentage1 = percentageOfRadiusReduction * length(clustersAtRadius)
            cARPercentage2 = length(clustersAtRadius) - cARPercentage1

        end 

        if siteTestedId < maxSiteIndex
            siteTestedId += 1
        else
            # If no suitable site is found, divide by 2 the percentage of clusters which allocation distance must be reduced to allow the addition of a site
            percentageOfRadiusReduction /= 2.0

            # If no site covering a client at the radius has been found, add any site
            if cARPercentage1 < 1
                cARPercentage1 = 0
            else
                cARPercentage1 = percentageOfRadiusReduction * length(clustersAtRadius)  
            end 
            
            cARPercentage2 = length(clustersAtRadius) - cARPercentage1
            siteTestedId = 2
        end 
    end

    return openedSites[1:addedSites], radius
end 

"""
Return the radius of an instance for a set of opened sites

Input
- instance: the instance
- sites: the index of the opened sites
"""
function getRadius(instance, sites)
    radius = -1

    for i in 1:size(instance.clientsCoordinates, 2)
         d = distanceToOpenSites(instance, i, sites)
        if d > radius
             radius = d
         end
    end

    return radius
end 


"""
Compute the linear relaxation of the p-center by considering only the clients already considered in the instance

- instance: the instance (in which the clients may already include the partition of the clients
- params: structure which define some additional solution parameters
- useImprovedRounding: true if a fractional solution is rounded by opening the sites by considering them in decreasing order of their y_j variables and if a site is opened only if it improves the radius (without improvement, a site is opened even if it does
- initialTimeElapsed: time elapsed in the resolution before calling this function 
"""
function computeLB(instance::Instance; params::ExpeParam=ExpeParam(false), initialTimeElapsed::Float64=0.0, useImprovedRounding::Bool=false)

    useDomination = params.useDomination
    params.useDomination = false
    
    startingTime = time()
    
    initializeDandOf!(instance, params)

    isClusterUpdated = true
    radius = nothing
    dichotomyResults = nothing
    dichotomyTime = 0
    clusterUpdateTime = 0
    alternativeSolutionsTime = 0
    radiusOfRepresentatives = nothing

    if instance.lb < instance.ub
        println("\n-- Computing LB")

        relaxationStartingTime = time()
        ## Relaxation resolution
        # Solve p-center problem associated to the clusters in the instance
        print("Solving with ", length(instance.clusters), " clusters... ")

        dichotomyResults = solveByDichotomy(instance, isInitialized=true, isRelaxation=true, params = params, useImprovedRounding=useImprovedRounding)
        dichotomyTime += dichotomyResults["resolutionTime"]
        
        println("... done in ", round(Int, time()-relaxationStartingTime), "s. Radius found: ", dichotomyResults["radius"])

        lbImproved = instance.lb < dichotomyResults["radius"]
        if lbImproved
            instance.lb = dichotomyResults["radius"]
            println("Set LB to ", dichotomyResults["radius"])
        end
    end
    
    totalTime= round(Int, time() - startingTime + initialTimeElapsed)

    results = Dict{String, Any}()

    results["resolutionTime"] = time() - startingTime

    params.useDomination = useDomination
    return results
end 


"""
Count the number of digits of a number

- Input
x: the number
b: the base considered
"""
function num_digits(x, b)
    return floor(Int, log(x) / log(b)) + 1
end
