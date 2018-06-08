# fitness and constraint functions of the creatures

@everywhere function fitness(ct::Creature)
    # it's the fitness function of the corresponding creature
    # ---> return a value that descripe the fitness
    # ---> to the environment of the creature

    vars = ct.vars
    fit = vars[1]+2*vars[2]+3*vars[3]+4*vars[4]+5*vars[5]
    fit = abs(42-fit)
    return fit
end

@everywhere function constraint(vars::Array{Float64})
    # it's the constraint function consist of constraint conditions
    # ---> if satisfied, then return true
    # ---> if unsatisfied, return fasle

    return true
end
