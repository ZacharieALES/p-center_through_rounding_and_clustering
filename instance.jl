"""
Structure which represents an instance of the p-center problem
"""
mutable struct Cluster

    # instanceId of the clients in the cluster. This list is ordered in increasing order.
    # (the instanceId of a client corresponds to row in instance.clientsCoordinates which contains the coordinates of the client)
    clientsId::Vector{Int64}

    # Coordinates of the representative of its cluster
    representative::Vector{Float64}

    # Id of the representative in the rows of the array Instance.clientsCoordinates
    representativeId::Int64

    # Partition of clientsId according to their quadrants
    # Each client is represented by its clusterId
    # (the clusterId of a client corresponds to its index in the vector clientsId of the cluster)
    # quadrants[q][i] indicates the cluster id of the ith client in the qth quadrant of the cluster
    quadrants::Vector{Vector{Int64}}

    function Cluster()
        return new()
    end
end

"""
Structure which represents an instance of the p-center problem
"""
mutable struct Instance
    n::Int64 # Number of clients
    m::Int64 # Number of candidate centers
    p::Int64 # Number of centers to open
    sitesCoordinates::Matrix{Float64}
    clientsCoordinates::Matrix{Float64} # clientsCoordinates[:, i] are the coordinates of the ith client
    clusters::Vector{Cluster}
    d::Matrix{Int64} # Distance matrix d[i, j] = distance between cluster i and site j
    dRows::Int # Number of rows of d in use (i.e., number of clusters whose distances are in d)    
    of::Vector{Int64} # Ordered distances between clients and sites

    computedDistances::Int
    lb::Int64 # lower bound on the optimal radius
    ub::Int64 # upper bound on the optimal radius
    ubOpenedSites::Vector{Int64} # the p sites opened in the solution of value ub

    moduloUB::Int64 # Upper bound on the optimal radius without rounding the distances (only used in solveByModulo)
    moduloOpenedSites::Vector{Int64} # Sites opened in the solution that provides moduloUB (only used in solveByModulo)
    areClientsIdenticalToSites::Bool

    siteDomination::Vector{Int64} # siteDomination[j] is the id of the first site which dominates site j or 0 if it is not dominated

    initialClusterCount::Int

    function Instance()
        return new()
    end
end

"""
Constructor
"""
function Cluster(instance::Instance, clientsId::Vector{Int64})

    this = Cluster()
    this.clientsId = clientsId
    this.representative = Vector{Float64}(zeros(2))

    distToBarycenters = updateRepresentative!(this, instance)

    # Distances between the points of each quadrant and the barycenter
    quadrantsDistances = Vector{Vector{Float64}}([])
    this.quadrants = Vector{Vector{Int64}}([])

    for i in 1:4
        push!(this.quadrants, Vector{Int64}([]))
        push!(quadrantsDistances, Vector{Float64}([]))
    end

    for (clientClusterId, clientInstanceId) in enumerate(this.clientsId)

        if instance.clientsCoordinates[1, clientInstanceId] > this.representative[1]
            if instance.clientsCoordinates[2, clientInstanceId] > this.representative[2]
                push!(this.quadrants[1], clientClusterId)
                push!(quadrantsDistances[1], distToBarycenters[clientClusterId])
            else
                push!(this.quadrants[2], clientClusterId)
                push!(quadrantsDistances[2], distToBarycenters[clientClusterId])
            end
        elseif instance.clientsCoordinates[2, clientInstanceId] > this.representative[2]
            push!(this.quadrants[3], clientClusterId)
            push!(quadrantsDistances[3], distToBarycenters[clientClusterId])
        else
            if (@view instance.clientsCoordinates[:, clientInstanceId]) != this.representative
                push!(this.quadrants[4], clientClusterId)
                push!(quadrantsDistances[4], distToBarycenters[clientClusterId])
            end
        end
    end

    # Sort the quadrants so that the first clients are the furthest from the barycenter of the cluster
    # (to consider them first when removing clients)
    for quadrantId in 1:4
        this.quadrants[quadrantId] = this.quadrants[quadrantId][sortperm(quadrantsDistances[quadrantId], rev=true)]
    end

    return this
end



"""
Constructor from an input file
"""
function Instance(path::String) 

    clientsCoordinates, sitesCoordinates, p = readInstanceFile(path)
    this = Instance()
    this.p = p
    this.computedDistances = 0
    this.dRows = 0
    this.clientsCoordinates = transpose(clientsCoordinates)
    this.sitesCoordinates = transpose(sitesCoordinates)
    this.n = size(this.clientsCoordinates, 2)
    this.m = size(this.sitesCoordinates, 2)
    this.lb = typemin(Int)
    this.ub = typemax(Int)
    this.moduloUB = typemax(Int)
    this.ubOpenedSites = Vector{Int}([])
    this.moduloOpenedSites = Vector{Int}([])
    
    createUnitaryClusters(this)

    # Test if the clients are identical to the sites (in that case there are no dominated sites)
    this.areClientsIdenticalToSites = false

    if this.m == length(this.clusters)
        this.areClientsIdenticalToSites = true
        siteId = 1

        while siteId <= this.m && this.areClientsIdenticalToSites
            if (@view this.clientsCoordinates[:, siteId]) != (@view this.sitesCoordinates[:, siteId])
                this.areClientsIdenticalToSites = false
            end
            siteId += 1
        end
    end

    return this
end

"""
Append to `values` the coordinates found in `content`, a piece of a coordinate block such as
"12661.6667 86437.2222; 12162.5000 86755.0000". Returns the number of rows appended.
"""
function appendCoordinates!(values::Vector{Float64}, content::AbstractString)

    rowCount = 0

    for rowString in split(content, ';')
        row = strip(rowString)
        if isempty(row)
            continue
        end 

        valueCountBefore = length(values)
        for token in split(row)
            push!(values, parse(Float64, token))
        end

        if length(values) - valueCountBefore != 2
            error("expected 2 coordinates per row, got \"", row, "\"")
        end
        rowCount += 1
    end

    return rowCount
end

"""
Resolve the coordinates named `name`, following concatenations and aliases.
"""
function resolveCoordinates(name::AbstractString, blocks::Dict{String, Matrix{Float64}},
                            concatenations::Dict{String, Vector{String}},
                            aliases::Dict{String, String}, depth::Int = 0)

    if depth > 10
        error("cycle while resolving \"", name, "\" in the instance file")
    end 

    if haskey(blocks, name)
        return blocks[name]
    end 

    if haskey(concatenations, name)
        return reduce(vcat, [resolveCoordinates(part, blocks, concatenations, aliases, depth + 1)
                             for part in concatenations[name]])
    end

    if haskey(aliases, name)
        return resolveCoordinates(aliases[name], blocks, concatenations, aliases, depth + 1)
    end 

    error("\"", name, "\" is not defined in the instance file")
end

"""
Read an instance file and return (clientsCoordinates, sitesCoordinates, p), the coordinates being
one row per point as they appear in the file.

Supported forms:
    n = 3496                          scalar
    clientsCoordinates = [            coordinate block, possibly spanning many lines
       12661.6667 86437.2222;
       ... ]
    sitesCoordinates = clientsCoordinates              alias
    clientsCoordinates = [part1; part2; part3]         concatenation of blocks defined above
"""
function readInstanceFile(path::String)

    blocks = Dict{String, Matrix{Float64}}()
    concatenations = Dict{String, Vector{String}}()
    aliases = Dict{String, String}()
    scalars = Dict{String, Int}()

    # Block currently being read ("" when we are not inside a block)
    currentName = ""
    values = Vector{Float64}()
    rowCount = 0

    function closeBlock!()
        blocks[currentName] = permutedims(reshape(values, 2, rowCount))
        currentName = ""
        values = Vector{Float64}()
        rowCount = 0
    end

    for rawLine in eachline(path)
        line = strip(rawLine)
        if (isempty(line) || startswith(line, "#"))
            continue
        end 

        if currentName != ""

            # Inside a coordinate block: it ends with the line holding the closing bracket
            content = line
            isClosed = endswith(content, "]")
            if isClosed
                content = content[1:prevind(content, lastindex(content))]
            end 
            rowCount += appendCoordinates!(values, content)
            if isClosed
                closeBlock!()
            end 
            continue
        end

        assignment = match(r"^(\w+)\s*=\s*(.*)$", line)
        if assignment === nothing
            continue
        end 
        name = String(assignment.captures[1])
        rhs = strip(assignment.captures[2])

        if startswith(rhs, "[")
            content = rhs[nextind(rhs, firstindex(rhs)):end]
            isClosed = endswith(content, "]")
            if isClosed
                content = content[1:prevind(content, lastindex(content))]
            end 

            parts = filter(!isempty, strip.(split(content, r"[;,]")))

            # [part1; part2] : a concatenation of blocks defined earlier.
            # Every part must be a bare identifier: testing for the mere presence of a letter would
            # misread a coordinate written in scientific notation ("0.00000e+00") as a block name.
            if !isempty(parts) && all(part -> occursin(r"^[A-Za-z_]\w*$", part), parts)

                concatenations[name] = String.(parts)
            else
                currentName = name
                rowCount += appendCoordinates!(values, content)
                if isClosed
                    closeBlock!()
                end 
            end
        elseif occursin(r"^-?\d+$", rhs)
            scalars[name] = parse(Int, rhs)
        elseif occursin(r"^\w+$", rhs)
            aliases[name] = String(rhs)
        end
    end

    if currentName != ""
        error("unterminated block \"", currentName, "\" in ", path)
    end 

    clientsCoordinates = resolveCoordinates("clientsCoordinates", blocks, concatenations, aliases)
    sitesCoordinates = resolveCoordinates("sitesCoordinates", blocks, concatenations, aliases)

    # p is optional in the instance files
    p = get(scalars, "p", 2)

    return clientsCoordinates, sitesCoordinates, p
end

"""
Create one cluster for each client
"""
function createUnitaryClusters(instance)
    instance.clusters = Vector{Cluster}(undef, instance.n)
    
    for clientId in 1:instance.n
        instance.clusters[clientId] = Cluster(instance, [clientId])
    end
end

"""
Maximal amount of spare capacity allocated at once when the distance matrix grows.

The matrix gains rows as clusters are disaggregated. Reallocating it at every disaggregation
copies the whole matrix each time, which is quadratic in the number of rounds, so capacity is
allocated in advance instead. Plain doubling cannot be used: one row is 8*m bytes, i.e. about
6 MB when m = 750000, so doubling from 2000 rows would ask for 12 GB on top of the 12 GB already
allocated. Growth is therefore geometric only while the spare capacity stays within this budget,
and linear (by this many bytes) beyond it.
"""
const MAX_SPARE_CAPACITY_BYTES = 256 * 1024 * 1024

"""
Set to `rowCount` the number of rows of instance.d which are in use, allocating more capacity if
needed. The values already in the matrix are preserved and the new rows are set to zero.
"""
function setDRows!(instance::Instance, rowCount::Int)

    capacity = size(instance.d, 1)

    if rowCount > capacity

        rowBytes = 8 * instance.m
        maxSpareRows = max(1, div(MAX_SPARE_CAPACITY_BYTES, max(rowBytes, 1)))

        # Geometric growth, capped so that at most MAX_SPARE_CAPACITY_BYTES is allocated at once
        newCapacity = max(rowCount, min(2 * capacity, capacity + maxSpareRows))

        newD = Matrix{Int64}(undef, newCapacity, instance.m)
        copyto!(view(newD, 1:instance.dRows, :), view(instance.d, 1:instance.dRows, :))
        instance.d = newD
    end

    if rowCount > instance.dRows
        fill!(view(instance.d, instance.dRows+1:rowCount, :), 0)
    end

    instance.dRows = rowCount
end

"""
Initialize matrix "d" and vector "of" according to the clusters of the instance
"""
function initializeDandOf!(instance::Instance, params::ExpeParam; modulo::Int=1, time_limit::Int=-1)

    # Initialize instance.d, reusing the capacity already allocated if there is enough of it
    # (this function is called once per modulo iteration, on clusters which keep growing)
    clusterCount = length(instance.clusters)

    if isdefined(instance, :d) && size(instance.d, 2) == instance.m && size(instance.d, 1) >= clusterCount
        instance.dRows = 0
        setDRows!(instance, clusterCount)
    else
        instance.d = Matrix{Int64}(zeros(clusterCount, instance.m))
        instance.dRows = clusterCount
    end

    startingUpdateTime = time()

    # Compute all values of instance.d
    for cId in 1:length(instance.clusters)
        for siteId in 1:instance.m
            updateD!(instance, cId, siteId, modulo=modulo)
        end
    end
    
    startingTime = time()
    instance.siteDomination = Vector{Int}(zeros(instance.m))

    startingDominationTime = time()
    if params.useDomination  && !instance.areClientsIdenticalToSites
        initializeDominations(instance, time_limit=remainingTime(time_limit, startingTime))
    end
    dominationTime = time() - startingDominationTime

    # Initialize instance.of
    ofInitStartTime = time()
    distances = Set{Int}()
    if !isOver(time_limit,startingUpdateTime) 
        for siteId in 1:instance.m

 
            # Add the distance between the site and the cluster only if the cluster is not dominated
            if instance.siteDomination[siteId] == 0
                for cId in 1:length(instance.clusters)
                    push!(distances, instance.d[cId, siteId])
                end 
            end 
        end
    end

    # of must never be empty: solveByDichotomy indexes of[1] before doing anything else.
    # It would be empty if the time limit interrupted the loop above.
    if isempty(distances)
        push!(distances, max(instance.lb, 0))
    end
    instance.of = Vector{Int64}(undef, length(distances))

    for (id, dist) in enumerate(distances)
        instance.of[id] = dist
    end

    sort!(instance.of)

    return dominationTime
end

"""
Test if a site is dominated by another

Input
- instance: the instance
- idSite1: site for which we test if it is dominated
- idSite2: site for which we test if it dominates site idSite1
- startingCluster: we only test the domination for cluster of id at least startingCluster
"""
function isSiteDominated(instance::Instance, idSite1::Int64, idSite2::Int64; startingCluster::Int64=1)
    isS1Dominated = true
    clusterId = startingCluster

    # While all the clients have not been tested and site 1 is dominated by site 2
    # (use instance.dRows rather than length(instance.clusters) since when new clusters are created their distances have not yet been added in instance.d. It is not useful to test the distances to the new clusters in that case since we only test if non dominated sites become dominated and adding client cannot make them dominated)
    while clusterId <= instance.dRows && isS1Dominated

        if instance.d[clusterId, idSite2] != instance.lb && instance.d[clusterId, idSite1] < instance.d[clusterId, idSite2]
            isS1Dominated = false
        end

        clusterId += 1
    end

    return isS1Dominated
end 

"""
After computing an optimal p-center solution (radius, openedSites) for the clusters of the instance:
- remove from each cluster the clients which are at a distance > radius of their closest opened site;
(each of them becomes a new cluster)
- update the representatives of the clusters

Input
- openedSitesSolutions: each element of this vector is a set of p opened sites. Its size is 1 if alternative solutions are not found by local search
"""
function updateClusters!(instance::Instance, radiusOfRepresentatives::Int64, openedSitesSolutions::Vector{Vector{Int}}, lbImproved::Bool; isIntResolution::Bool, params::ExpeParam=ExpeParam(false), modulo::Int=1, fractionalSitesOpened::Vector{Tuple{Int, Float64}}=Vector{Tuple{Int, Float64}}([]), time_limit::Int=-1, improveModuloUB::Bool=false)

    startingTime = time()
    solutionsCount = length(openedSitesSolutions)
    initialClusterCount = min(instance.n, instance.p + params.clusterCountIncrease)

    minimalRadiusOnRepresentatives = radiusOfRepresentatives

    # If the openedSitesSolutions was rounded from a fractional solution, it does not necessarily give the optimal radius on the representatives
    if !isIntResolution
        minimalRadiusOnRepresentatives = 0
    end 
    
    radiusForAllClients = Vector{Int}(minimalRadiusOnRepresentatives .* ones(solutionsCount))

    # If the radius on the original distances must be computed
    moduloRadiusForAllClients = Vector{Int}()
    if improveModuloUB && modulo > 1
        moduloRadiusForAllClients = Vector{Int}(minimalRadiusOnRepresentatives .* ones(solutionsCount))
    end 

    previousClusterCount = length(instance.clusters)
    # For each initial cluster (i.e., clusters which are not reduced to one client)
    for (clusterId, cluster) in enumerate(@view instance.clusters[1:initialClusterCount])

        if length(cluster.clientsId) > 1

            # Get the id of the clients to remove from the clusters
            # - in instance.clientsCoordinates (clientsInstanceIdToRemove); and
            # - in cluster.clientsId (clientsClusterIdToRemove)
            clientsInstanceIdToRemove, clientsClusterIdToRemove = removeQuadrantClientsFromCluster(cluster, params, instance, openedSitesSolutions, radiusOfRepresentatives, radiusForAllClients, moduloRadiusForAllClients, isIntResolution, modulo=modulo)

            # If there is a fractional solution provided and it is fractional (i.e., more than p variables are positives)
            if length(fractionalSitesOpened) > instance.p
                fractionalClientsInstanceIdToRemove, fractionalClientsClusterIdToRemove = removeClientsFromFractionalSites(cluster, instance, fractionalSitesOpened, modulo=modulo)

                for (clientRemoveId, clientInstanceId) in enumerate(fractionalClientsInstanceIdToRemove)
                    if !(clientInstanceId in clientsInstanceIdToRemove) # If it is not already planned to remove this client
                        push!(clientsInstanceIdToRemove, clientInstanceId)
                        push!(clientsClusterIdToRemove, fractionalClientsClusterIdToRemove[clientRemoveId])
                    end
                end 
            end 
            
            # If clients must be removed
            if length(clientsInstanceIdToRemove) > 0

                # Order the clients according to their cluster id
                permutation = sortperm(clientsClusterIdToRemove)
                clientsInstanceIdToRemove = clientsInstanceIdToRemove[permutation]
                clientsClusterIdToRemove = clientsClusterIdToRemove[permutation]
                
                addedClusters = Vector{Cluster}(undef, length(clientsClusterIdToRemove))
                
                # For each client to remove
                for (clientId, clientClusterId) in Iterators.reverse(enumerate(clientsClusterIdToRemove))

                    # Get its ids in the cluster list and in the instance coordinates
                    clientInstanceId = clientsInstanceIdToRemove[clientId]

                    # Remove from the cluster 
                    deleteat!(cluster.clientsId, clientClusterId)

                    # All clients cluster id > to the cluster id of the client removed must be decremented
                    for quadrant in cluster.quadrants
                        for id in eachindex(quadrant)
                            if quadrant[id] > clientClusterId
                                quadrant[id] -= 1
                            end
                        end
                    end

                    addedClusters[clientId] = Cluster(instance, [clientInstanceId])

                end
                # Create a new cluster reduced to one client
                append!(instance.clusters, addedClusters)
            end 
        end         
    end
    
    # Get the radius of the first solution
    openedSites = openedSitesSolutions[1]
    
    # It is the worst radius for this solution over all the clusters
    bestRadiusForAllClients = radiusForAllClients[1]

    for (viewSolutionId, solutionRadius) in enumerate(@view radiusForAllClients[2:end])
        
        if solutionRadius < bestRadiusForAllClients
            bestRadiusForAllClients = solutionRadius
            openedSites = openedSitesSolutions[viewSolutionId + 1]
        end 
    end

    # If we update moduloUB
    if improveModuloUB && modulo > 1
        # Use the representative of each cluster to update the radius of the solutions
        for (clusterId, cluster) in enumerate(instance.clusters)
            for solutionId in 1:length(moduloRadiusForAllClients)

                # If the radius without rounding of this solution can still be better than instance.moduloB
                if moduloRadiusForAllClients[solutionId] < instance.moduloUB 
                    cModuloDist = distanceToOpenSitesUB(instance, cluster.representativeId, openedSitesSolutions[solutionId], stoppingValue = moduloRadiusForAllClients[solutionId])
                    if cModuloDist > moduloRadiusForAllClients[solutionId]
                        moduloRadiusForAllClients[solutionId] = cModuloDist
                    end 
                end
            end 
        end

        # Get the solution that returns the best radius
        bestModuloRadius = moduloRadiusForAllClients[1]
        bestIndex = 1
        for (viewSolutionId, solutionRadius) in enumerate(@view moduloRadiusForAllClients[2:end])
            
            if solutionRadius < bestModuloRadius
                bestModuloRadius = solutionRadius
                bestIndex = viewSolutionId + 1
            end 
        end

        if bestModuloRadius < instance.moduloUB
            instance.moduloUB = bestModuloRadius
            instance.moduloOpenedSites = openedSitesSolutions[bestIndex]
            println("moduloUB improved in update cluster to: ", instance.moduloUB)
        end
    end
    newClustersCreated = previousClusterCount != length(instance.clusters)
    isOptFound = isIntResolution && !newClustersCreated

    # If the radius is better than the UB and the optimal value is not reached
    ubImproved = bestRadiusForAllClients < instance.ub

    if ubImproved
        println("Set UB to ", bestRadiusForAllClients)
        setUB!(instance, bestRadiusForAllClients, openedSites, params)
    end 

    startDomTime = time()
    if !isOptFound
        if params.useDomination && (lbImproved || ubImproved)
            testBoundDominations(instance, time_limit=remainingTime(time_limit, startingTime))
        end 
    end
    dominationTime = time() - startDomTime

    # If new clusters have been created
    if newClustersCreated

        # Update d and of for the new clusters
        # (rows are taken from the spare capacity of d, so most of the time nothing is copied)
        setDRows!(instance, length(instance.clusters))

        for newClusterId in previousClusterCount+1:length(instance.clusters)

            # Update d for this new cluster
            for siteId in 1:instance.m
            updateD!(instance, newClusterId, siteId, modulo=modulo)

            # Only update vector of for non dominated sites  
            if instance.siteDomination[siteId] == 0
                
                # clusterIdIsLast is true since we know that all the clusters associated to distances in instance.of have an id < newClusterId
                updateOf!(instance, newClusterId, siteId)
            end 
        end
    end 
        
   startDomTime = time()
   # Test if the existing dominations still apply after adding the new clusters
   if params.useDomination && !instance.areClientsIdenticalToSites && !isOver(time_limit, startingTime)

       mayBecomeDominated = Vector{Bool}(zeros(instance.m))
       
       for (idSite1, site1Domination) in enumerate(instance.siteDomination) # Update the domination according to the added clusters  
           if isOver(time_limit, startingTime)
               break
           end

           # If the site was dominated
           if site1Domination != 0

               # If the addition of the new clusters make that site idSite1 is not dominated anymore by the same site
               if !isSiteDominated(instance, idSite1, site1Domination, startingCluster = previousClusterCount+1)
                   
                   mayBecomeDominated[idSite1] = true

                   # Update the bounds of this site for clients 1:previousClusterCount
                   # (for the other clients the bounds have already been taken into account in updateD!)
                   for clusterId in 1:previousClusterCount
                       instance.d[clusterId, idSite1] = boundedDistance(instance, instance.d[clusterId, idSite1])
                   end
               end 
           end 
       end

       becameNonDominated = Vector{Bool}(zeros(instance.m))

       for (idSite1, site1MayBecomeDominated) in enumerate(mayBecomeDominated) # Update the domination according to the added clusters  
           if isOver(time_limit, startingTime)
               break
           end

           # If the site was dominated
           if site1MayBecomeDominated

               # Test the next sites for domination
               
               instance.siteDomination[idSite1] = 0
               firstSiteToTest = 1

               # While all the sites have not been tested and site idSite1 is not dominated 
               for (viewIdSite2, site2Domination) in enumerate(@view instance.siteDomination[firstSiteToTest:end])
                   idSite2 = viewIdSite2 + firstSiteToTest - 1
                   if idSite1 != idSite2 && site2Domination == 0 && isSiteDominated(instance, idSite1, idSite2)
                       instance.siteDomination[idSite1] = idSite2
                       break
                   end
               end
               
               # If the site is not dominated anymore
               if instance.siteDomination[idSite1] == 0
                   becameNonDominated[idSite1] = true
               end 
           end 
       end

       # Test if each site idSite1 which became non dominated is not dominated by a site which became non dominated after idSite1 (i.e., a site which has an id > idSite1)
       for (idSite1, site1BecameNonDominated) in Iterators.reverse(enumerate(@view becameNonDominated[1:instance.m-1]))
           if isOver(time_limit, startingTime)
               break
           end
           if site1BecameNonDominated
               
               for (viewIdSite2, site2BecameNonDominated) in enumerate(@view becameNonDominated[idSite1+1:end])
                   idSite2 = viewIdSite2+idSite1
                   if site2BecameNonDominated && isSiteDominated(instance, idSite1, idSite2)
                       instance.siteDomination[idSite1] = idSite2
                       becameNonDominated[idSite1] = false
                       break
                   end
               end
           end 
       end
       
       # Update of the site that are not dominated anymore
       for (siteId, siteBecameNonDominated) in enumerate(becameNonDominated)
           if isOver(time_limit, startingTime)
               break
           end
           if siteBecameNonDominated
               for clusterId in 1:length(instance.clusters)
                   updateOf!(instance, clusterId, siteId)
               end
           end 
       end
   end
   dominationTime = time() - startDomTime
   end

   return previousClusterCount != length(instance.clusters), dominationTime
end



"""
Remove at most one client in each quadrant
"""
function removeQuadrantClientsFromCluster(cluster::Cluster, params::ExpeParam, instance::Instance, openedSitesSolutions::Vector{Vector{Int64}}, radiusOfRepresentatives::Int64, radiusForAllClients::Vector{Int64}, moduloRadiusForAllClients::Vector{Int64}, isIntResolution::Bool; modulo::Int=1)

    # Id of the clients to remove from the clusters
    # ... in instance.clientsCoordinates
    clientsInstanceIdToRemove = Vector{Int64}([])

    # ... in cluster.clientsId
    clientsClusterIdToRemove = Vector{Int64}([])

    solutionsCount = length(openedSitesSolutions)

    improveModuloUB = length(moduloRadiusForAllClients) > 0

    # Remove at most one client in each quadrant for each solution
    for (quadrantId, quadrant) in enumerate(cluster.quadrants)

        invalidClientId = Vector{Int}(.-ones(solutionsCount))

        # For each client of the quadrant
        for (clientQuadrantId, clientClusterId) in enumerate(quadrant)

            clientInstanceId = cluster.clientsId[clientClusterId]

            # For each solution
            for (solutionId, invalidId) in enumerate(invalidClientId)

                # If there are currently no invalid client in this quadrant for this solution, test if the distance of the client to its closest site is <= radiusOfRepresentatives
                # Otherwise, only test if it is <= radiusForAllClients[solutionId]
                stoppingValue = invalidId == -1 ? radiusOfRepresentatives : radiusForAllClients[solutionId]

                # Get the distance to its closest opened site
                cDist = distanceToOpenSites(instance, clientInstanceId, openedSitesSolutions[solutionId], modulo=modulo, stoppingValue = stoppingValue)

                # If its radius is larger than the radius of the representatives
                # (or if the solution is fractional since in that case the radius of an alternative solution can be < radiusOfRepresentatives)
                if !isIntResolution || cDist > radiusOfRepresentatives

                    # If it is even larger than the radius of all the previously tested clients for this solution
                    if cDist > radiusForAllClients[solutionId]
                        radiusForAllClients[solutionId] = cDist
                    end

                    # If it is the first invalid client found for this radius
                    if invalidId == -1 && cDist > radiusOfRepresentatives
                        if !(clientQuadrantId in invalidClientId)

                            # If this client has not already been added to the list of removed clients for a previous solution
                            if length(clientsInstanceIdToRemove) == 0 || clientsInstanceIdToRemove[end] != clientInstanceId
                                push!(clientsClusterIdToRemove, clientClusterId)
                                push!(clientsInstanceIdToRemove, clientInstanceId)
                            end 
                        end
                        invalidClientId[solutionId] = clientQuadrantId
                    end
                end
                
                if improveModuloUB
                    # If the radius without rounding of this solution can still be better than instance.moduloB
                    if moduloRadiusForAllClients[solutionId] < instance.moduloUB
                        cModuloDist = distanceToOpenSitesUB(instance, clientInstanceId, openedSitesSolutions[solutionId], stoppingValue = moduloRadiusForAllClients[solutionId])

                        if cModuloDist > moduloRadiusForAllClients[solutionId]
                            moduloRadiusForAllClients[solutionId] = cModuloDist
                        end 
                    end 
                end 
            end 
        end # for clientQuadrantId in 1:length(quadrant)

        # If invalid clients are found, remove them from the quadrant
        uniqueSortedClientId = sort(unique(invalidClientId))

        for (invalidId, invalidClient) in Iterators.reverse(enumerate(uniqueSortedClientId))
            if invalidClient != -1
                deleteat!(quadrant, invalidClient)
            end 
        end
    end # for quadrantId in 1:length(cluster.quadrants)

    return clientsInstanceIdToRemove, clientsClusterIdToRemove
end

"""
Compute the distance of a client to its closest open site

Input
- (optional) stoppingValue: stop testing the opened sites if an opened site at distance <= stoppingValue is found (used to stop early when testing if the client does satisfy the current radius of the cluster representatives)
"""
function distanceToOpenSites(instance::Instance, clientId::Int64, openedSites::AbstractArray{Int64}; modulo::Int=1, stoppingValue::Int64=typemin(Int64), isClusterId::Bool=false)

    bestDist = -1

    if isClusterId
        bestDist = instance.d[clientId, openedSites[1]]
    else
        bestDist = distance(instance, clientId, openedSites[1], modulo=modulo)
    end 

    stoppingValue = max(instance.lb, stoppingValue)
    siteId = 2

    while siteId <= length(openedSites) && bestDist > stoppingValue

        newDist = -1

        if isClusterId
            newDist = instance.d[clientId, openedSites[siteId]]
        else
            newDist = distance(instance, clientId, openedSites[siteId], modulo=modulo)
        end
        
        if newDist < bestDist
            bestDist = newDist
        end 

        siteId += 1
    end

    return bestDist
end 


"""
Compute the distance between a cluster representative and a site.
Add this distance to instance.d
"""
function updateD!(instance, clusterId, siteId; modulo::Int=1)

    c = instance.clusters[clusterId]
    dist = distance(instance, c.representativeId, siteId, modulo=modulo)

    if dist < instance.lb
        dist = instance.lb
        
    elseif dist > instance.ub + 1 && instance.ub < typemax(Int)
        dist = instance.ub + 1
    end

    instance.d[clusterId, siteId] = dist
    
end


"""
Add the distance between a cluster and a site to instance.of 
"""
function updateOf!(instance, clusterId, siteId)
    
    dist = instance.d[clusterId, siteId]

    # Get the id of the first value in "of" which is >= dist
    distInsertId = searchsortedfirst(instance.of, dist)

    # If dist is not included in "of"
    if distInsertId > length(instance.of) || instance.of[distInsertId] != dist

        # Add it
        insert!(instance.of, distInsertId, dist)
    end 
end 

"""
Set the lower bound of the instance to a given value (and update the distances accordingly)
"""
function setLB!(instance::Instance, lb::Int64, params::ExpeParam)
    instance.lb = lb

    # For each non-dominated site 
    # for j  in 1:size(instance.d, 2)
    for (j, jDomination) in enumerate(instance.siteDomination)

        if jDomination == 0
            
            # For each cluster
            for c in 1:instance.dRows

                if instance.d[c, j] < lb
                    instance.d[c, j] = lb
                end
            end 
        end
    end

    # Add the lb to instance.of if it not already in it and return its id in instance.of
    lbId = addToSortedVector!(instance.of, lb)

    # Remove the previous distances
    instance.of = instance.of[lbId:end]

end

"""
Set the upper bound of the instance to a given value (and update the distances accordingly)
"""
function setUB!(instance::Instance, ub::Int64, openedSites::Vector{Int64}, params::ExpeParam)

    instance.ub = ub
    instance.ubOpenedSites = copy(openedSites)
    
    # For each non-dominated site
    for (j, jDomination) in enumerate(instance.siteDomination)
        if jDomination == 0
            
            # For each cluster
            for c in 1:instance.dRows
                if instance.d[c, j] > ub + 1
                    instance.d[c, j] = ub + 1
                end
            end 
        end
    end
    
    # Add ub+1 to instance.of if it is not already in it and return its id in instance.of
    ubId = addToSortedVector!(instance.of, ub + 1)

    # Remove the previous distances
    instance.of = instance.of[1:ubId]        
end

"""
Add a value to a vector which is sorted if it is not already in it
"""
function addToSortedVector!(vector::Vector{Int64}, value::Int64)
    id = searchsortedfirst(vector, value)

    if id > length(vector) || vector[id] != value
        insert!(vector, id, value)
    end

    return id
end 

"""
Add the values of a sorted vector into another sorted vector
"""
function addToSortedVector!(vector::Vector{Int64}, vectorValues::Vector{Int64})
    for value in vectorValues
        addToSortedVector!(vector, value)
    end 
end 

"""
Test if after moving bounds, non-dominated site become dominated

Input
- instance: the instance
"""
function testBoundDominations(instance::Instance; time_limit::Int=-1)

    startingTime = time()
    
    if !instance.areClientsIdenticalToSites
        
        for (idSite1, site1Domination) in enumerate(instance.siteDomination)
            if isOver(time_limit, startingTime)
                break
            end

            # If the site was previously non-dominated
            if site1Domination == 0


                # While all the sites have not been tested and site idSite1 is not dominated 
                for (idSite2, site2Domination) in enumerate(instance.siteDomination)

                    id2Dominated = true

                    # Only test if idSite1 is dominated by idSite2 if idSite2 is not already dominated
                    # Otherwise, if the distances of these two sites are all equal, the could both dominate each other and none of them would be considered
                    if idSite1 != idSite2 && site2Domination == 0 && isSiteDominated(instance, idSite1, idSite2)
                        instance.siteDomination[idSite1] = idSite2
                        break
                    end
                end
            end
        end 
    end
end 


"""
Test if each site is dominated by another for the set of clusters in the instance

Input
- instance: the instance
- computeUndominatedDistances: true if the distance of undominated sites must be computed (only done once at the initialization; not done afterwards when the bounds on the radius are modified since the distance have already been computed at the initialization)
"""
function initializeDominations(instance::Instance; time_limit::Int=-1)

    startingTime = time()

    for (idSite1, site1Domination) in enumerate(instance.siteDomination)

        # Otherwise test all the potential dominators
        if isOver(time_limit, startingTime)
            break
        end 

        # While all the sites have not been tested and site idSite1 is not dominated 
        for (idSite2, site2Domination) in enumerate(instance.siteDomination)

            # Only test if idSite1 is dominated by idSite2 if idSite2 is not already dominated
            # Otherwise, if the distances of these two sites are all equal, the could both dominate each other and none of them would be considered
            if idSite1 != idSite2 && site2Domination == 0 && isSiteDominated(instance, idSite1, idSite2)
                instance.siteDomination[idSite1] = idSite2
                break
            end
        end
    end
end 

"""
Compute the distance between a client and a site
"""
function distance(instance::Instance, clientId::Int, siteId::Int; modulo::Int=1)

    instance.computedDistances += 1
    d = distance((@view instance.clientsCoordinates[:, clientId]), (@view instance.sitesCoordinates[:, siteId]), modulo = modulo)
    return boundedDistance(instance, d)
end

"""
Clip a distance in [instance.lb, instance.ub+1]
"""
function boundedDistance(instance::Instance, dist::Int)
    
    if dist < instance.lb
        dist = instance.lb
    elseif instance.ub != typemax(Int) && dist > instance.ub + 1
        dist = instance.ub + 1
    end
    return dist
    
end 

"""
Compute a distance between two points
"""
function distance(pointA::AbstractArray{Float64}, pointB::AbstractArray{Float64}; modulo::Int=1)
    
    distance = round(Int, sqrt((pointA[1] - pointB[1])^2 + (pointA[2] - pointB[2])^2))

    if modulo == 1
        return distance
    else
        return distance - rem(distance, modulo)
    end 
end 

"""
Update the representative of a cluster
"""
function updateRepresentative!(c::Cluster, instance::Instance)

    # Distance of each client in the cluster to the barycenter
    distToBarycenters = Vector{Float64}(undef, length(c.clientsId))
    
    if length(c.clientsId) == 1
        c.representative = @view(instance.clientsCoordinates[:, c.clientsId[1]])
        c.representativeId = c.clientsId[1]
    else
        
        barycentreX = round(Int, sum(instance.clientsCoordinates[1, c.clientsId]) / Float64(length(c.clientsId)))
        barycentreY = round(Int, sum(instance.clientsCoordinates[2, c.clientsId]) / Float64(length(c.clientsId)))

        closestClientClusterId = 1
        closestDistance = abs(barycentreX - instance.clientsCoordinates[1, c.clientsId[1]])^2 + abs(barycentreY - instance.clientsCoordinates[2, c.clientsId[1]])^2

        distToBarycenters[1] = closestDistance

        # enumerate starts at 1 on the view, while the client it yields sits at position viewId+1
        # in c.clientsId: distToBarycenters must be indexed by the position, not by the view index
        for (viewId, clientId) in enumerate(@view c.clientsId[2:end])
            clientClusterId = viewId + 1
            dist = abs(barycentreX - instance.clientsCoordinates[1, clientId])^2 + abs(barycentreY - instance.clientsCoordinates[2, clientId])^2
            
            distToBarycenters[clientClusterId] = dist

            if dist < closestDistance
                closestDistance = dist
                closestClientClusterId = clientClusterId
            end 
        end

        c.representativeId = c.clientsId[closestClientClusterId]
        c.representative = @view instance.clientsCoordinates[:, c.representativeId]
    end

    return distToBarycenters
end 

"""
Remove clients from the fractional solution returned by the dichotomy

Input:
- sitesOpened: tuples (j, yj) with j a site index and yj the value of yj in the fractional solution. Sorted by decreasing values of yj
"""
function removeClientsFromFractionalSites(cluster::Cluster, instance::Instance, sitesOpened::Vector{Tuple{Int, Float64}}; modulo::Int=1)

    # Id of the clients to remove from the clusters
    # ... in instance.clientsCoordinates
    clientsInstanceIdToRemove = Vector{Int64}([])

    # ... in cluster.clientsId
    clientsClusterIdToRemove = Vector{Int64}([])

    # Remove at most one client in each quadrant for each solution 
    for (quadrantId, quadrant) in enumerate(cluster.quadrants)

        invalidClientId = -1

        # For each client of the quadrant (while an invalid client has not been found)
        for (clientQuadrantId, clientClusterId) in enumerate(quadrant)

            clientInstanceId = cluster.clientsId[clientClusterId]

            # Distance to the client and value in the fractional solution of the closest sites of the client in increasing value of the distance
            # closestSites[4][1]: distance between the client and the fourth closest opened sites
            # closestSites[4][2]: value in the fractional solution of this site
            closestSites = Array{Tuple{Int, Float64}}([])

            # Sum of the fractional values of the variables of the sites in closestSites
            currentCovering = 0

            # Sum of the distances to the closest sites weighted by their values
            currentRadius = 0

            # distance to the furthest site in closestSites
            worstDistance = typemax(Int)

            #println("Client ", clientInstanceId)
            
            # For each opened site
            for tupleJ in sitesOpened
                j = tupleJ[1]
                yj = tupleJ[2]

                dij = distance(instance, clientInstanceId, j, modulo=modulo)

                # If the site can be involved in the covering of client i
                if dij < worstDistance || currentCovering < 1

                    # Add j to the covering of client i
                    currentRadius += dij * yj 
                    currentCovering += yj 
                    id = searchsortedfirst(closestSites, dij, by = v -> v[1])

                    # Insert a new entry in closestSites if the distance is not already in closest sites
                    if id > length(closestSites) || closestSites[id][1] != dij
                        insert!(closestSites, id, (dij, yj))
                    else # Otherwise increase the covering of the tuple of distance dij
                        closestSites[id] = Tuple{Int, Float64}((dij, closestSites[id][2] + yj))
                    end  

                    if id == length(closestSites)
                        worstDistance = dij
                    end

                    # If the client is too much covered, remove the furthest covering site(s)
                    if currentCovering > 1
                        additionalCovering = currentCovering-1 # quantity of covering that must be removed
                        lastIndex = length(closestSites) 
                        
                        while additionalCovering > 0

                            furthestTuple = closestSites[lastIndex]
                            
                            # If the furthest site covers i more than the additional covering, reduce its covering of i
                            if furthestTuple[2] > additionalCovering 
                                currentCovering = 1
                                closestSites[lastIndex] = Tuple{Int, Float64}((furthestTuple[1], furthestTuple[2] - additionalCovering))
                                currentRadius -= furthestTuple[1] * additionalCovering
                                additionalCovering = 0

                            else # Otherwise, remove it completely from the covering of i
                                additionalCovering -= furthestTuple[2]
                                lastIndex -= 1
                                currentRadius -= furthestTuple[1] * furthestTuple[2]
                                currentCovering -= furthestTuple[2]
                            end
                        end

                        closestSites = closestSites[1:lastIndex]
                        worstDistance = closestSites[lastIndex][1]
                    end

                    if currentCovering > 1-1E-4

                        # Stop trying to consider closer sites if the radius is lower than the LB
                        # (in that case the client can not be invalid)
                        if currentRadius <= instance.lb + 1E-4
                            break
                        end
                    end 
                end # if dij < worstDistance || currentCovering < 1
            end # for each site opened

            if currentRadius > instance.lb + 1E-4   

                # Add this client in the list of clients to remove
                push!(clientsClusterIdToRemove, clientClusterId)
                push!(clientsInstanceIdToRemove, clientInstanceId)
                invalidClientId = clientQuadrantId
                break
            end 
        end # end for each client in the quadrant

        # If an invalid client is found, remove it from the quadrant
        if invalidClientId != -1
            deleteat!(quadrant, invalidClientId) 
        end
    end # end for each quadrant

    return clientsInstanceIdToRemove, clientsClusterIdToRemove
end 


"""
Compute the distance of a client to its closest open site

Input
- (optional) stoppingValue: stop testing the opened sites if an opened site at distance <= stoppingValue is found (used to stop early when testing if the client does satisfy the current radius of the cluster representatives)
"""
function distanceToOpenSitesUB(instance::Instance, clientId::Int64, openedSites::AbstractArray{Int64}; stoppingValue::Int64=typemin(Int64))

    bestDist = -1
    bestDist = distanceUB((@view instance.clientsCoordinates[:, clientId]), (@view instance.sitesCoordinates[:, 1]))

    stoppingValue = max(instance.lb, stoppingValue)
    siteId = 2

    while siteId <= length(openedSites) && bestDist > stoppingValue

        newDist = distanceUB((@view instance.clientsCoordinates[:, clientId]), (@view instance.sitesCoordinates[:, siteId]))
        
        if newDist < bestDist
            bestDist = newDist
        end 

        siteId += 1
    end

    return bestDist
end 


function distanceUB(pointA::AbstractArray{Float64}, pointB::AbstractArray{Float64})
    return round(Int, sqrt((pointA[1] - pointB[1])^2 + (pointA[2] - pointB[2])^2)) 
end 
