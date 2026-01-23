"""
Find alternative solutions of p opened sites which do not increase the radius and would be satsfied by more clients which are not yet considered

Input
- instance: the instance
- openedSites: index of the sites opened in the solution
- radius: radius of the solution
- modulo: all distances computed will be rounded down to 10^(modulo-1)
- time_limit: remaining time
"""
function findAlternativeOpenedSites(instance::Instance, openedSites::AbstractArray{Int64}, radius::Int64, modulo::Int; time_limit::Int=-1)
    
    alternativeOpenedSites = Vector{Vector{Int}}()
    startingTime = time()

    # Initialize the structure used to efficiently perform and evaluate the moves
    # (a move consists in exchanging a site in the solution with another non-dominated one)
    # (a move is invalid if it deteriorates the radius of the cluster representatives)
    status = LocalSearchStatus(instance, openedSites, radius, modulo)

    # If all the clients are covered at the radius
    if length(status.setC) == 0
        return alternativeOpenedSites
    end
    
    movesCount = length(status.currentSites) * length(status.nonDominatedSites)
    hasValidMoves = true
    coverAllSetC = false
    
    iterationWithoutImprovement = 0
    
    while iterationWithoutImprovement < status.nIter && hasValidMoves && !coverAllSetC && !isOver(time_limit, startingTime)

        ## Perturbation: perform alpha random moves
        for i in 1:status.alpha

            # Generate a random move
            removedSite = rand(1:length(status.currentSites))
            addedSite = rand(1:length(status.nonDominatedSites))
            attempt = 1

            # Find the next valid move
            while !status.isMoveValid[removedSite, addedSite] && attempt <= movesCount && !isOver(time_limit, startingTime)
                attempt += 1
                if addedSite < length(status.nonDominatedSites)
                    addedSite += 1
                else
                    addedSite = 1
                    if removedSite < length(status.currentSites)
                        removedSite += 1
                    else
                        removedSite = 1
                    end
                end 
            end

            if attempt > movesCount
                hasValidMoves = false
            else 
                applyMove!(status, instance, removedSite, addedSite, radius, modulo)
            end 
        end

        ## Local search: perform the best move until no more move improves the objective of covering setC
        solutionImproved = true

        while solutionImproved && hasValidMoves && !isOver(time_limit, startingTime)

            bestMoves = Vector{Tuple{Int, Int}}()
            bestMoveValue = nothing

            solutionImproved = false

            # For each valid move
            for currentSiteId in 1:length(status.currentSites)
                for nonDomSiteId in 1:length(status.nonDominatedSites)
                    if status.isMoveValid[currentSiteId, nonDomSiteId]

                        # If the move is among the best known and improves the coverage of setC
                        if ((bestMoveValue == nothing || status.movesValue[currentSiteId, nonDomSiteId] >= bestMoveValue) && status.movesValue[currentSiteId, nonDomSiteId] > 0)

                            # If it is the first move with this value found
                            if (bestMoveValue == nothing || status.movesValue[currentSiteId, nonDomSiteId] > bestMoveValue)
                                bestMoveValue = status.movesValue[currentSiteId, nonDomSiteId]
                                bestMoves = Vector{Tuple{Int, Int}}([(currentSiteId, nonDomSiteId)])
                            else
                                push!(bestMoves, (currentSiteId, nonDomSiteId))
                            end 
                        end 
                    end 
                end
            end

            if length(bestMoves) > 0
                # Randomly select a move among the best ones and apply it
                (removedSite, addedSite) = bestMoves[rand(1:length(bestMoves))]
                
                applyMove!(status, instance, removedSite, addedSite, radius, modulo)
                solutionImproved = true 
            end 
        end # while solutionImproved

        # If a better solution is obtained, update the best solution
        if status.currentValue > status.bestSiteValue 
            status.bestSites = copy(status.currentSites)
            status.bestSiteValue = status.currentValue

            # If all the elements in setC are covered by the new solution
            if status.bestSiteValue == length(status.setC)
                previousSize = length(status.setC)
                addBetaClientsToSetC!(status, instance, radius, modulo)

                if length(status.setC) == previousSize
                    coverAllSetC = true
                end 
            end

            push!(alternativeOpenedSites, status.bestSites)
            iterationWithoutImprovement = 0
        else
            iterationWithoutImprovement += 1
        end 
    end

    return alternativeOpenedSites
end 

"""
Structure that represents the status of the local search
"""
mutable struct LocalSearchStatus

    # Opened sites which cover the most clients of setC
    bestSites::Vector{Int}

    # Number of sites of setC covered by bestSites
    bestSiteValue::Int
    
    # Sites for which we try to find a better solution
    setC::Vector{Int}

    # Sites in the current solution
    currentSites::Vector{Int}

    # Number of elements of setC which are covered by currentSites
    currentValue::Int

    # clusterCoveredBy[c] contains the list of sites in the current solution which distance to the representative of cluster c is <= radius
    # (it is used to know if a move deteriorate the radius of the representatives)
    # (a site is represented by its index in currentSites)
    # (the sites associated to a cluster are in increasing order)
    clusterCoveredBy::Vector{Vector{Int}}

    # setCCoveredBy[i] contains the list of sites in the current solution which distance to setC[i] is <= radius
    # (it is used to know the value of a move)
    # (a site is represented by its index in currentSites)
    # (the sites associated to an element in setC are in increasing order)
    setCCoveredBy::Vector{Vector{Int}}

    # movesValue[i, j] = difference of the number of clients from setC covered when exchanging currentSites[i] and nonDominatedSites[j]
    movesValue::Matrix{Int64}

    # isMoveValid[i, j] = true iff exchanging sites currentSites[i] and nonDominatedSites[j] in the current solution still covers all the representative with a radius of at most radius and if site nonDominatedSites[j] is not already in the solution
    isMoveValid::Matrix{Bool}
    
    # Number of clients not covered considered
    beta::Int

    # Number of steps of the perturbation step
    alpha::Int

    # Maximal number of steps without improvement of the best solution
    nIter::Int

    candidates::Vector{Int}
    
    # isSetCCandidate[i] = 1 iff client i is not the representative of its cluster and if is it not already in setC
    setCCandidateWeights::Vector{Float64}

    # Id of the sites which are not dominated (i.e., for which instance.siteDomination == 0)
    nonDominatedSites::Vector{Int}
    
    function LocalSearchStatus()
        return new()
    end
end 

"""
Initialization of the local search status

Input
- instance: the instance
- openedSites: index of the sites opened in the solution
- radius: radius of the solution
- modulo: all distances computed will be rounded down to 10^(modulo-1)
"""
function LocalSearchStatus(instance::Instance, openedSites::AbstractArray{Int64}, radius::Int, modulo::Int)

    this = LocalSearchStatus()

    this.beta = 5
    this.alpha = 3
    this.nIter = 3

    this.currentSites = copy(openedSites)
    this.bestSites = copy(openedSites)
    
    ## Find the clients which are not representative of their clusters
    # (required since clients in setC must not be representative of their cluster)
    candidateCount = 0
    
    for (clusterId, cluster) in enumerate(@view instance.clusters[1:instance.initialClusterCount])
        cluster = instance.clusters[clusterId]
        for quadrant in cluster.quadrants
            candidateCount += length(quadrant)
        end
    end
    
    this.candidates = Vector{Int}(undef, candidateCount)

    candidateId = 1
    
    for (clusterId, cluster) in enumerate(@view instance.clusters[1:instance.initialClusterCount])
        for quadrant in cluster.quadrants
            for clientClusterId in quadrant
                this.candidates[candidateId] = cluster.clientsId[clientClusterId]
                candidateId += 1
            end 
        end
    end 

    # Initialize setC
    this.setC = Vector{Int64}()
    this.setCCoveredBy = Vector{Vector{Int}}()

    # Moves are only evaluated for non dominated sites
    this.nonDominatedSites = findall(instance.siteDomination .== 0)
    
    # Initialize moves
    this.movesValue = Matrix{Int64}(zeros(length(this.currentSites), length(this.nonDominatedSites)))
    this.isMoveValid = Matrix{Bool}(ones(length(this.currentSites), length(this.nonDominatedSites)))

    addBetaClientsToSetC!(this, instance, radius, modulo)

    # If all clients are covered at the radius
    if length(this.setC) == 0
        return this
    end 
    
    # None of the clients added to setC will be covered by the initial solution
    this.bestSiteValue = 0
    this.currentValue = this.bestSiteValue
    
    # Initialize clusterCoveredBy
    this.clusterCoveredBy = Vector{Vector{Int}}(undef, length(instance.clusters))

    for cId in 1:length(instance.clusters)
        cCenters = Vector{Int}()
        for (siteId, currentSite) in enumerate(this.currentSites)
            if instance.d[cId, currentSite] <= radius
                push!(cCenters, siteId)
            end
        end 
        this.clusterCoveredBy[cId] = cCenters
    end

    # Initialize isMoveValid (movesValue has already been initialized in addBetaClientsToSetC)
    # For each non-dominated site
    for (nonDomSiteId, nonDominatedSite) in enumerate(this.nonDominatedSites)

        isInCurrentSite = nonDominatedSite in this.currentSites
        
        # For each site in the current solution
        for currentSiteId in 1:length(this.currentSites)

            if isInCurrentSite
                this.isMoveValid[currentSiteId, nonDomSiteId] = false
            else
                nonDomSiteInM = nonDominatedSite
                this.isMoveValid[currentSiteId, nonDomSiteId] = isMoveValid(this, instance, currentSiteId, nonDomSiteId, radius)
            end 
        end
    end

    return this
end 

"""
Exchange two sites in the current solution.

Input
- status: current status of the local search
- instance: the instance
- removedSiteInP: index in {1, ..., p} of the site removed
- addedSite: index in status.nonDominatedSites of the site added
- radius: radius of the solution
- modulo: all distances computed will be rounded down to 10^(modulo-1)
"""
function applyMove!(status, instance, removedSiteInP, addedSite, radius, modulo)

    removedSiteInM = status.currentSites[removedSiteInP]
    addedSiteInM = status.nonDominatedSites[addedSite]
    status.currentSites[removedSiteInP] = addedSiteInM

    ## Update status.isMoveValid
    # testMovesValidityOfSite[j] is true if a nonDominatedSite[j] could have invalid moves which become valid
    testMovesValidityOfSite = Vector{Bool}(zeros(length(status.nonDominatedSites)))

    # The removed site could now have valid moves
    removedSiteInNDS = findfirst(status.nonDominatedSites .== removedSiteInM)

    if removedSiteInNDS != nothing
        testMovesValidityOfSite[removedSiteInNDS] = true
    end 

    for (clusterId, currentClusterCoveredBy) in enumerate(status.clusterCoveredBy)

        removedIndex = searchsortedfirst(currentClusterCoveredBy, removedSiteInP)
        isCoveredByRemovedSite =  removedIndex <= length(currentClusterCoveredBy) && currentClusterCoveredBy[removedIndex] == removedSiteInP
        isCoveredByAddedSite = instance.d[clusterId, addedSiteInM] <= radius
        
        if isCoveredByRemovedSite

            # If the cluster is covered by one less site 
            if !isCoveredByAddedSite

                deleteat!(currentClusterCoveredBy, removedIndex)
                
                # If there was two sites that were covering this cluster
                # (the number of current sites covering the cluster goes from 2 to 1)
                if length(currentClusterCoveredBy) == 1

                    coveringSiteInP = currentClusterCoveredBy[1]
                    
                    # Each valid move involving the only remaining site covering the cluster...
                    for (nonDomSiteId, nonDominatedSite) in enumerate(status.nonDominatedSites)
                        
                        # ... becomes invalid if it does not cover the cluster
                        if status.isMoveValid[coveringSiteInP, nonDomSiteId] && instance.d[clusterId, nonDominatedSite] > radius
                            status.isMoveValid[coveringSiteInP, nonDomSiteId] = false
                        end 
                    end
                end 
            end
        # If the cluster is covered by one more site
        elseif isCoveredByAddedSite
            insert!(currentClusterCoveredBy, removedIndex, removedSiteInP)

            # If there is at most another site in the solution that covers the cluster
            # (the number of current sites covering the cluster goes from 0 to 1 or from 1 to 2) 
            if length(currentClusterCoveredBy) <= 2

                # Each invalid move which adds a site not covering the cluster could now be valid
                for (nonDomSiteId,  nonDominatedSite) in enumerate(status.nonDominatedSites)
                    if instance.d[clusterId, nonDominatedSite] > radius
                        testMovesValidityOfSite[nonDomSiteId] = true
                    end
                end 
            end 
        end
    end # for clusterId in 1:length(instance.clusters)

    # For each non dominated sites which could have invalid moves that become valid
    for (nonDomSiteId, testMoves) in enumerate(testMovesValidityOfSite)
        if testMoves && !(status.nonDominatedSites[nonDomSiteId] in status.currentSites)
            for currentSiteInP in 1:length(status.currentSites)
                if !status.isMoveValid[currentSiteInP, nonDomSiteId]
                    status.isMoveValid[currentSiteInP, nonDomSiteId] = isMoveValid(status, instance, currentSiteInP, nonDomSiteId, radius)
                end 
            end 
        end 
    end

    # The moves involving the added site are not valid anymore
    for currentSiteInP in 1:length(status.currentSites)
        status.isMoveValid[currentSiteInP, addedSite] = false
    end 
    
    ## Update status.movesValue
    for (setCId, coveredBy) in enumerate(status.setCCoveredBy)

        removedIndex = searchsortedfirst(coveredBy, removedSiteInP)
        isCoveredByRemovedSite =  removedIndex <= length(coveredBy) && coveredBy[removedIndex] == removedSiteInP
        isCoveredByAddedSite = distance(instance, status.setC[setCId], addedSiteInM, modulo=modulo) <= radius

        if isCoveredByRemovedSite

            # If the element of setC is covered by one less site
            if !isCoveredByAddedSite
                
                deleteat!(coveredBy, removedIndex)
                
                # If it was the only site which was covering the element of setC
                # (the number of current sites covering the element goes from 1 to 0)
                if length(coveredBy) == 0

                    status.currentValue -= 1

                    for (nonDomSiteId, ndSites) in enumerate(status.nonDominatedSites)

                        coverSetCElement = distance(instance, status.setC[setCId], ndSites, modulo=modulo) <= radius
                        if coverSetCElement
                            # Increment the value of each move leading to the cover of the element of setC  
                            for currentSiteInP in 1:length(status.currentSites)
                                status.movesValue[currentSiteInP, nonDomSiteId] += 1
                            end
                        else
                            # Increment the value of each move removing removedSiteInP and not leading to the cover of the element of setC  
                            status.movesValue[removedSiteInP, nonDomSiteId] += 1
                        end
                    end 

                # If there was two sites that were covering this cluster
                # (the number of current sites covering the element goes from 2 to 1)
                elseif length(coveredBy) == 1

                    coveringSiteInP = coveredBy[1]
                    
                    # Each move involving the only remaining site covering the cluster...
                    for (nonDomSiteId, ndSites) in enumerate(status.nonDominatedSites)
                        
                        # ... is decremented if the added site does not cover the element of setC
                        if distance(instance, status.setC[setCId], ndSites, modulo=modulo) > radius
                            status.movesValue[coveringSiteInP, nonDomSiteId] -= 1
                        end 
                    end
                end 
            end 

            # If the element of setC is covered by one more site        
        elseif isCoveredByAddedSite

            # If there was previously only one site covering this element of setC
            # (the number of current sites covering the element goes from 1 to 2) 
            if length(coveredBy) == 1

                onlySiteInP = coveredBy[1]

                # Each move involving the only site previously covering the element of setC...
                for (nonDomSiteId, ndSites) in enumerate(status.nonDominatedSites)
                    # ... and in which the added site does not cover the element of setC are incremented
                    if distance(instance, status.setC[setCId], ndSites, modulo=modulo) > radius
                        status.movesValue[onlySiteInP, nonDomSiteId] += 1
                    end
                end 
            end
                        
            # If there was previously no site covering this element of setC
            # (the number of current sites covering the element goes from 0 to 1)
            if length(coveredBy) == 0

                status.currentValue += 1

                for (nonDomSiteId, ndSites) in enumerate(status.nonDominatedSites)
                    
                    coverSetCElement = distance(instance, status.setC[setCId], ndSites, modulo=modulo) <= radius

                    if coverSetCElement

                        # Decrement the value of each move leading to the cover of the element of setC   
                        for currentSiteInP in 1:length(status.currentSites) 
                            status.movesValue[currentSiteInP, nonDomSiteId] -= 1
                        end
                    else
                        # Decrement the value of each move removing removedSiteInP and not leading to the cover of the element of setC
                        status.movesValue[removedSiteInP, nonDomSiteId] -= 1
                    end
                end 
            end
            insert!(coveredBy, removedIndex, removedSiteInP)
        end                          
    end

    ## Recompute movesValue for the removed site

    # Get the number of elements of setC covered by removedSiteInNDS but not by the current sites
    positiveDelta = 0

    for (setCId, cElement) in enumerate(status.setC)
        if length(status.setCCoveredBy[setCId]) == 0 && distance(instance, cElement, removedSiteInM, modulo=modulo) <= radius
            positiveDelta += 1
        end 
    end

    negativeDeltas = Vector{Int}(zeros(length(status.currentSites)))

    for (setCId, coveredBy) in enumerate(status.setCCoveredBy)

        # If the element of setC is only covered by one site in the current solution and not covered by the removed site
        if length(coveredBy) == 1 && distance(instance, status.setC[setCId], removedSiteInM, modulo=modulo) > radius
            # Exchanging the removed site with this one site would remove the covering of setC[setCId]
            negativeDeltas[coveredBy[1]] -= 1
        end 
    end

    if removedSiteInNDS != nothing
        for (currentSiteId, negativeDelta) in enumerate(negativeDeltas)
            status.movesValue[currentSiteId, removedSiteInNDS] = positiveDelta + negativeDelta
        end
    end 
end

"""
Test if a move is valid
"""
function isMoveValid(status, instance, removedSiteInP, nonDomSiteId, radius)

    clusterId = 1
    isValid = true

    #while clusterId <= length(instance.clusters) && isValid
    for (clusterId, coveredBy) in enumerate(status.clusterCoveredBy)
        
        # If cluster clusterId is only covered by the removed site
        if length(coveredBy) == 1 && coveredBy[1] == removedSiteInP

            # If adding site nomDomSiteId would not cover it
            if instance.d[clusterId, status.nonDominatedSites[nonDomSiteId]] > radius
                isValid = false
                break
            end
        end
    end
    return isValid
end 

"""
Add clients to set C
"""
function addBetaClientsToSetC!(status::LocalSearchStatus, instance::Instance, radius::Int, modulo::Int)

    startingTime = time()
    desiredLength = length(status.setC) + status.beta

    testedCandidates = Vector{Bool}(zeros(length(status.candidates)))
    selectedIds = Vector{Int}()
    testedCount = 0
    candidatesIdAddedToSetC = Vector{Int}()
    addToSetC = Vector{Bool}(zeros(status.beta))

    # While enough clients have not been found and there are still candidates
    # Add a time condition since this loop is sometime very long for an unknown reason...
    while length(status.setC) < desiredLength && testedCount < length(status.candidates)

        # Test if a random set of candidates have a larger radius than the representatives
        # (considering a set of clients enables parallelism which would not be possible on the parent while loop)
        candidatesIdToTest = Vector{Int}()

        # For each remaining client to test
        for clientToTest in 1:desiredLength-length(status.setC)

            if testedCount < length(status.candidates)
                # Get a random client from the candidates which is not already in candidatesIdToTest
                candidateId = rand(1:length(status.candidates))

                while testedCandidates[candidateId]
                    candidateId += 1
                    if candidateId > length(status.candidates)
                        candidateId = 1
                    end
                end

                testedCount += 1
                testedCandidates[candidateId] = true
                push!(candidatesIdToTest, candidateId)
            end 
        end

        # Reinitialize addToSetC
        addToSetC = Vector{Bool}(zeros(status.beta)) 
        
        for (candidateId2, candidateId1)  in enumerate(candidatesIdToTest)
            #candidateId1 = candidatesIdToTest[candidateId2]
            
            if distanceToOpenSites(instance, status.candidates[candidateId1], status.currentSites, modulo = modulo, stoppingValue = radius) > radius
                addToSetC[candidateId2] = true
            end
        end

        for (candidateId2, addCandidate2) in Iterators.reverse(enumerate(addToSetC))
            if addCandidate2

                candidateId1 = candidatesIdToTest[candidateId2]
                clientId = status.candidates[candidateId1]

                push!(candidatesIdAddedToSetC, candidateId1)
                
                # If yes, add it to setC
                push!(status.setC, clientId)

                # This client is covered by no site in the current solution
                push!(status.setCCoveredBy, Vector{Int}())

                # Increment the move value of any site covering this client
                for (nonDomSiteId, ndSites)  in enumerate(status.nonDominatedSites)
                    if distance(instance, clientId, ndSites, modulo=modulo) <= radius
                        for currentSiteInP in 1:length(status.currentSites)
                            status.movesValue[currentSiteInP, nonDomSiteId] += 1
                        end
                    end 
                end
            end 
        end
    end

    # Remove the nodes added to setC from the future candidates
    sort!(candidatesIdAddedToSetC)

    for (candidateId2, addedToSetC) in Iterators.reverse(enumerate(candidatesIdAddedToSetC))
         deleteat!(status.candidates, addedToSetC)
    end
end 

