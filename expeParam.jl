mutable struct ExpeParam

    clusterCountIncrease::Int # The number of clusters will be equal to p + clusterCountIncrease
    useDomination::Bool

    usePMaxConstraint::Bool # True if constraint sum xi <= p is used in the set cover problem
    useNullObjective::Bool # True if the objective of the set cover is null (in that case the constraint sum xi <= p is used)
    useCutOff::Bool # True if the cutoff parameter of the set cover is set to p and the resolution stops after the first integer solution

    findAlternativeSolutions::Bool # True if the local search algorithm is used to find alternative solutions which have the same radius for more non-considered clients
    function ExpeParam()
        return new()
    end
end

# Constructeur de la structure
function ExpeParam(useDomination::Bool; clusterCountIncrease::Int=2, usePMaxConstraint::Bool=false, useNullObjective::Bool=false, useCutOff::Bool=false, findAlternativeSolutions::Bool=false)
               
    this = ExpeParam()
    this.useDomination = useDomination
    this.clusterCountIncrease = clusterCountIncrease
    this.usePMaxConstraint = usePMaxConstraint
    this.useNullObjective = useNullObjective
    this.useCutOff = useCutOff
    this.findAlternativeSolutions = findAlternativeSolutions
    
    return this                                                                                                                         
end 
