# it is a genetic algorithm moudle

module GAMD
using Plots

export Population, Creature, garun

@everywhere struct Population
    # it's the type of the whole population of the creatures
    # parameters

    ncre :: Int64 # population numbers
    nmax :: Int64 # max generation
    nmin :: Int64 # min generation
    npar :: Int64 # parents numbers
    # nchl :: Int64 # children numbers
    mr   :: Float64 # initial mutation probability
    rpr  :: Float64 # reproduction parobability of good parents
                    # while the bad parents reproduction rate is (1-rpr)
end

@everywhere mutable struct Creature
    # it's the type of the individual creature

    # parameters
    vars :: Array{Float64}

    # functions
    fitness   :: Function

    # judge if the constraint condition is satisfied
    Creature(vars, fitness) = constraint(vars) ? new(vars, fitness) :
        new(ctvars(varmin, varmax), fitness)
end

@everywhere function ctvars(varmin::Array{Float64}, varmax::Array{Float64})
    # it's the function to initilize the single individual vars
    # ---> return the new random creature's vars

    dvar = varmax-varmin
    vars = zeros(length(dvar))
    for i in 1:length(dvar)
        vars[i] = varmin[i] + dvar[i]*rand() # random initilize Creature's vars
    end

    return vars
end


@everywhere function creator(p::Population, varmin::Array{Float64}, varmax::Array{Float64})
    # it's the function to initilize the population
    # ---> return the initilization of population

    # the creature paremeters range should be ensured here
    cts = Creature[] # innitilize the creature array

    for i in 1:p.ncre
        vars = ctvars(varmin, varmax) # random initilize Creatures' vars
        push!(cts,Creature(vars, fitness))
    end

    return cts
end

@everywhere function murate(gen, ngen, mutrate)
    # creatures mutate function
    return mutrate/2*cos(π/gen*ngen)+mutrate/2+0.05
end

@everywhere function calfit(ctso::Array{Creature})
    # it's the fitness calculation function
    # ---> return the fitness array

    cts = copy(ctso) # store the creature of population data into "cts"
    # fits = Float64[] # fitness value initilization
    fittmp = pmap(fitness, cts)

    # @sync @parallel for i in 1:length(cts)
    #     push!(fits, cts[i].fitness(cts[i]))
    # end
    #
    # fittmp = fetch(fits)

    return fittmp
end

@everywhere function selection(ctso::Array{Creature}, fitss::Array{Float64}, npare::Int64)
    # it's the selection function
    # ---> return the selected first level parents creatures for next generation
    # ---> return the others as second level parents creatures
    # ---> return the best fitness creature and its fitness values
    # ---> return the average fitness value of the current generation

    cts = copy(ctso)
    fits = copy(fitss)
    mincre = cts[indmin(fits)]
    min = minimum(fits)
    ave = mean(fits)

    pare = Creature[]
    godind = Int[] # record the index of parents
    for i in 1:npare # pull out the parents creatures
        ind = indmin(fits)
        push!(pare, cts[ind])
        append!(godind, ind)
        deleteat!(fits, ind)
    end

    for i in 1:npare
        deleteat!(cts, godind[i])
    end

    return pare, cts, mincre, min, ave
end

@everywhere function crossover(pare::Array{Creature}, cts::Array{Creature}, mr::Float64, rpr::Float64, λ = 0.6)
    # make the better creatures get more access to reproduce better children
    # make the worse creatures get less chance to reproduce children

    newcres = Creature[]
    npare = length(pare)

    # the good parents creatures reproduction
    for i in 1:npare
        if rand() > rpr
            vars = λ*pare[i].vars + (1-λ)*cts[rand(1:length(cts))].vars
            push!(newcres, Creature(vars, fitness))
        else
            for j in 1:npare
                if i==j
                    continue
                end

                if rand() < mr
                    push!(newcres, Creature(ctvars(varmin,varmax), fitness))
                else
                    vars = λ*pare[i].vars+(1-λ)*pare[j].vars
                    push!(newcres, Creature(vars, fitness))
                end

            end
        end
    end

    # the bad ones reproduction
    for i in 1:length(cts)
        if rand() > rpr
            vars = λ*cts[i].vars+(1-λ)*pare[rand(1:npare)].vars
            push!(newcres, Creature(vars, fitness))
        end
    end

    newcres = vcat(newcres, pare)
    return newcres
end

@everywhere function rms(minval::Float64, meanay::Array{Float64}, cfit)
    # it's the residual error function
    # ---> return the difference between the current fitness result with
    # ---> the mean fitness result of the current 10 generations

    if minval <=cfit
        return 0.
    else
        ave = mean(meanay)
        rmsv = abs(minval/ave-1)
        return rmsv
    end
end

@everywhere function convergence(gen::Int64, nmax::Int64, rmsv::Float64, judgement=1e-1)
    # it's the convergence function
    # the generation numbers and the rms is the judgement conditions
    # ---> if convergence succeed, return true
    # ---> if convergence fail, return false

    if gen>=nmax || rmsv<=judgement
        return true
    else
        return false
    end
end

@everywhere function gaplot(gen::Int64, nmin::Int64, min::Array{Float64}, ave::Array{Float64},
    rms::Array{Float64})
    # it's the plot function
    # ---> display the evolution curve real-time

    if gen >= 2
        x = 1:gen
        minplot = plot(x, min, label="minimum fitness", lw=6)
        maplot  = plot!(minplot, ave, seriestype=:scatter, label="average fitness",
            title="Envolution Curve", xlabel="generations", ylabel = "fitness", lw=3)
        if gen <= nmin
            display(maplot)
            savefig(pwd()*"//output//plots//envolution_curve_$(gen)")
        else
            y = nmin+1:gen
            rmsplot = plot(y, rms, xlabel = "generations",
                ylabel = "rms", label = "RMS", lw = 4)
            plottmp = plot(maplot, rmsplot, layout=(2,1))
            display(plottmp)
            savefig(pwd()*"//output//plots//envolution_curve_$(gen)")
        end
    end
end

@everywhere function garun(p::Population, varmin::Array{Float64},
    varmax::Array{Float64}, cfit = 1., λ=0.6, judgement = 1e-1)
    # it's the main run function
    # ---> run the whole GA project and give all the output you need


    # create output files
    logfile = open(pwd()*"//output//ga_log.txt","w");

    print("=======================================\n")
    print("====== GA PROCESS IS GOING AHEAD ======\n\n")
    write(logfile, "=======================================\n")
    write(logfile, "====== GA PROCESS IS GOING AHEAD ======\n\n")

    tic();
    gen = 1 # initilize the generation
    cts = creator(p, varmin, varmax) # initilize the creatures
    print("+++++++++++++++++ cts is $(cts) +++++++++++++++\n\n")
    p = Population(ncre, nmax, nmin, npar, mr, rpr)

    # parameters initilization
    minfit = Float64[]
    avefit = Float64[]
    mincre = Creature[]
    rmsrec = Float64[]

    while true
        # generation goes ahead

        # # cores calculation distribution
        # npop = length(cts)
        # npop_ = Int64(npop%np)
        # ncore_ = Int64((npop - npop_)/np)
        #
        # # ====== fitness function parallel computation ======
        # sel = Array{Any}(np)
        # for i in 1:np-1
        #     sel[i] = remotecall_fetch(calfit, i, cts[ncore_*(i-1)+1:ncore_*i])
        # end
        #
        # # sel1 = remotecall_fetch(calfit, 1, cts[ncore_*(1-1)+1:ncore_*1])
        # # sel2 = remotecall_fetch(calfit, 2, cts[ncore_*(2-1)+1:ncore_*2])
        # # sel3 = remotecall_fetch(calfit, 3, cts[ncore_*(3-1)+1:ncore_*3])
        #
        # sel[np] = remotecall_fetch(calfit, np, cts[ncore_*(np-1)+1:end])
        #
        # # sel4 = remotecall_fetch(calfit, np, cts[ncore_*(np-1)+1:end])
        # # sel[1] = sel1
        # # sel[2] = sel2
        # # sel[3] = sel3
        # # sel[4] = sel4
        # # ====== parallel computation complete ======
        #
        # fittmp = Float64[]
        # for i in 1:np # combine sel arraies into one array
        #     append!(fittmp, sel[i])
        # end

        fittmp = calfit(cts)
        print("++++++++++fittmp is $(fittmp) ++++++++++++\n\n")
        seltmp = selection(cts, fittmp, p.npar) # select the suitable parents

        pare = seltmp[1]
        chil = seltmp[2]
        push!(mincre, seltmp[3])
        push!(minfit, seltmp[4])
        push!(avefit, seltmp[5])

        cts = crossover(pare, chil, p.mr, p.rpr, λ) # mutate and crossover and generate new generation
        print("+++ No. $(gen) generation is generated successful +++\n")
        print("+++ The average value of the generation is $(avefit[gen]) +++\n")
        print("+++ The min value of the generation is $(minfit[gen]) +++\n")
        print("+++ The correspond creature of the minfit is $(mincre[gen]) +++\n\n")
        write(logfile, "+++ No. $(gen) generation is generated successful +++\n
        +++ The average value of the generation is $(avefit[gen]) +++\n
        +++ The min value of the generation is $(minfit[gen]) +++\n
        +++ The correspond creature of the minfit is $(mincre[gen]) +++\n\n")

        if gen > p.nmin
            rmstmp = rms(minfit[gen], minfit[end-p.nmin:end], cfit)
            push!(rmsrec, rmstmp)
            print("+++ The current generation RMS is $(rmstmp) +++\n\n")
            write(logfile,"+++ The current generation RMS is $(rmstmp) +++\n\n")
            if convergence(gen, p.nmax, rmstmp, judgement)
                print("====== GA PROCESS GET CONVERGED ======\n\n")
                write(logfile, "====== GA PROCESS GET CONVERGED ======\n\n")
                break
            end
        end

        gaplot(gen, p.nmin, minfit, avefit, rmsrec) # display the curves

        gen += 1
    end

    close(logfile);

    toc()
    return  minfit, avefit, mincre[gen], cts
end


end
