"""
Sets the tolerance for the iLQR solver to the intermediate tolerances, until the last 
outer loop iteration.
"""
function set_tolerances!(solver::ALSolver{T}, i::Int, cost_tol=solver.opts.cost_tolerance, 
                         grad_tol=solver.opts.gradient_tolerance) where T
    opts = options(solver)
    if i != solver.opts.iterations_outer
        opts.cost_tolerance = opts.cost_tolerance_intermediate
        opts.gradient_tolerance = opts.gradient_tolerance_intermediate
    else
        opts.cost_tolerance = cost_tol 
        opts.gradient_tolerance = grad_tol 
    end

    return nothing
end

function solve!(solver::ALSolver)
    reset!(solver)

    conset = get_constraints(solver)
    if !is_constrained(solver)
        ilqr = get_ilqr(solver)
        solve!(ilqr)
        terminate!(solver)
        return solver
    end

    # Terminal tolerances
    cost_tol = solver.opts.cost_tolerance
    grad_tol = solver.opts.gradient_tolerance

    Z̄ = solver.ilqr.Z̄
    for al_iter = 1:solver.opts.iterations_outer
        # Set potentially looser tolerances for inner iLQR solve 
        set_tolerances!(solver, al_iter, cost_tol, grad_tol)

        # Solve problem with iLQR
        solve!(solver.ilqr)

        # Check solver status
        status(solver) > SOLVE_SUCCEEDED && break

        # Evaluate the constraints and the cost for the current trajectory 
        J = TO.cost(solver, Z̄)  # NOTE: this evaluates the constraints
        c_max = max_violation(conset)
        μ_max = max_penalty(conset)
        record_iteration!(solver, J, c_max, μ_max)

        # Check if it's converged
        isconverged = evaluate_convergence(solver)
        if isconverged
            break
        end

        # Outer loop updates
        dualupdate!(conset)
        # lin_cons_update!(conset)
        penaltyupdate!(conset)

        # Reset iLQR solver
        # TODO: is this necessary? it gets reset at the beginning of the iLQR solve method 
        reset!(solver.ilqr)
    end
    # Reset tolerances to their original values
    solver.opts.cost_tolerance = cost_tol
    solver.opts.gradient_tolerance = grad_tol
    terminate!(solver)
    return solver
end

function lin_cons_update!(conset::ALConstraintSet)
    #need to get all polytopes as A,b matrices
    println("going through cons")
    #save the linear constraints into another array
    lincons_A = []
    lincons_b = []
    # large_λ_arr = Vector{Bool}()
    # idx_needs_change = Vector{Int}()
    
    for alcon in conset.constraints
        if isa(alcon.con, TO.LinearConstraint{})
            # println("linear constraint")
            # check if constraint unique, if yes, add to lincons
            if !any([A == alcon.con.A for A in lincons_A])
                # display(alcon.con.A)
                push!(lincons_A, alcon.con.A)
                push!(lincons_b, alcon.con.b)
            end
        end
    end
    # splice!(conset.constraints, 2)
    for alcon_idx in 1:length(conset.constraints)
        alcon = conset.constraints[alcon_idx]
        if isa(alcon.con, TO.LinearConstraint{})
            max_λ = maximum(alcon.λ[1])
            # println("Max lamda = ",max_λ)
            if max_λ > 30.0
                dist_to_obs = zeros(length(lincons_A))
                cur_state = state(alcon.Z[1][alcon.inds[1]])
                for i in 1:length(lincons_A)
                    if lincons_A[i] == alcon.con.A
                        dist_to_obs[i] = -1000.0
                    else
                        dist_to_obs[i] = minimum(lincons_A[i]*cur_state[1:3] .- lincons_b[i])
                    end
                end
                # print(dist_to_obs)
                max_val, max_idx = findmax(dist_to_obs)#find the most comfortable constraint
                if max_val > -1.0 #if the constraint is violated a little bit
                    println("Found replacing constraint")
                    stagecon = TO.LinearConstraint(13,4,lincons_A[max_idx], lincons_b[max_idx],Inequality(),1:3)
                    conset.constraints[alcon_idx] = ALConstraint{Float64}(alcon.Z[1], stagecon, alcon.inds, zeros(100)) 
                    # println("Found replacing constraint")
                    # @set alcon.con = TO.LinearConstraint(13,4,lincons_A[i], lincons_b[i],Inequality(),1:3)
                    # @set alcon.con.b = lincons_b[i]
                    # replace_con_ineq = 
                    # replace_con = ALConstraint{Float64}(Z̄, replace_con_ineq, alcon.inds)
                    # println("switching constraint")
                    # # simply remove this constraint and add a new one in the same indexed position
                    # #!!!![FATAL] [1687809539.259629]: Error: setfield! immutable struct of type ALConstraint cannot be changed

                    # # alcon.con.A = lincons_A[i]
                    # # alcon.con.b = lincons_b[i]
                    # splice!(conset.constraints, alcon_idx, replace_con)
                    break
                end
                
            end
        end
    end
    
                

    #             max_λ = maximum(alcon.λ[1])
    #             idx_needs_change = push!(idx_needs_change, alcon.inds[1])
    #             if max_λ > 1e5
    #                 push!(large_λ_arr, true)
    #             else
    #                 push!(large_λ_arr, false)
    #             end
    #         end
    #         push!(lincons, alcon.con)
    #         max_λ = maximum(alcon.λ[1])
    #         idx_needs_change = push!(idx_needs_change, alcon.inds[1])
    #         if max_λ > 1e5
    #             push!(large_λ_arr, true)
    #         else
    #             push!(large_λ_arr, false)
    #         end
    #     end
    # end
    # for  
        # RD.evaluate(con::LinearConstraint, z::AbstractKnotPoint) = con.A*RD.getdata(z)[con.inds] .- con.b
    # check if trajectory current position is in any other polytope
                
                #if yes, then switch the constraint to that polytope
        # println(typeof(alcon))
        # println(alcon.inds)
        # display(alcon.λ)

            


        # for i in eachindex(alcon.inds)
        #     λ, μ, c = alcon.λ[i], alcon.μ[i], alcon.vals[i]
            # conval = alcon.con[i]
            # println(typeof(conval))
            # if isa(conval, LinearConstraint{}) 
            #     println("Linear constraint max value is: ", maximum(conval.λ))
            #     println(typeof(conval))
            # end
        # end
        # c = conval.vals
        # λ = conval.λ
        # μ = conval.μ
	    # λ_max = conval.params.λ_max
	    # cone = TO.sense(conval.con)
	# λ_min = TO.sense(conval.con) == Equality() ? -λ_max : zero(λ_max)
        # for i in eachindex(conval.inds)
        #     λ[i] .= dual_update(cone, SVector(λ[i]), SVector(c[i]), SVector(μ[i]), λ_max) 
        # end
    
end

function record_iteration!(solver::ALSolver, J, c_max, μ_max)
    stats = solver.stats 
    record_iteration!(stats, c_max=c_max, penalty_max=μ_max, is_outer=true)
    lg = getlogger(solver)
    @log lg "iter" stats.iterations
    @log lg "AL iter" stats.iterations_outer
    @log lg "cost" J
    @log lg "||v||" c_max
    @log lg μ_max
end

function evaluate_convergence(solver::ALSolver)
    lg = getlogger(solver)
    iter = solver.stats.iterations
    isconverged = false
    if solver.stats.c_max[iter] < solver.opts.constraint_tolerance
        @log lg "info" "Constraint tolerance met."
        solver.stats.status = SOLVE_SUCCEEDED
        isconverged = true
    end
    if solver.opts.kickout_max_penalty && solver.stats.penalty_max[i] >= solver.opts.penalty_max
        @log lg "info" "Hit max penalty."
        isconverged = true
    end
    if iter >= solver.opts.iterations
        @log lg "info" "Hit max iterations."
        solver.stats.status = MAX_ITERATIONS
        isconverged = true
    end
    if solver.stats.iterations_outer >= solver.opts.iterations_outer
        @log lg "info" "Hit max AL iterations."
        solver.stats.status = MAX_ITERATIONS_OUTER
        isconverged = true
    end
    return isconverged
end