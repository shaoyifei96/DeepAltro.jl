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
    # @warn "AL: going into AL iterations"
    # @warn "print something else"
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
            # @warn "AL:converged"
            break
        end
        # @warn "AL:updating duals"
        # Outer loop updates
        dualupdate!(conset)
        penaltyupdate!(conset)
        lin_cons_update!(conset)

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
    debug_flag = false
    if debug_flag
        @warn "going through cons"
    end
    #save the linear constraints into another array
    lincons_A = []
    lincons_b = []
    lincon_idx = []
    if debug_flag
    idx_start_end = zeros(Int, 0, 2)
    end
    assign_matrix = zeros(Int, 100)
    # large_λ_arr = Vector{Bool}()
    # idx_needs_change = Vector{Int}()
    cons_counter_first_loop = 1
    for (idx,alcon) in enumerate(conset.constraints)
        if isa(alcon.con, TO.LinearConstraint{})
            if !isempty(alcon.inds) 
            #using idx here is wrong!! results in constraint idx (including other contraints) to be counted as polybox            # check if constraint unique, if yes, add to lincons
                push!(lincons_A, alcon.con.A)
                push!(lincons_b, alcon.con.b)
                push!(lincon_idx, idx)
                if debug_flag
                    idx_start_end = [idx_start_end; [alcon.inds[1], alcon.inds[end]]']
                end
            # else
            #     idx_start_end = [idx_start_end; [-1, -1]']
            # end
                assign_matrix[alcon.inds] .= cons_counter_first_loop
                cons_counter_first_loop += 1
            end
            
        end
    end
    if debug_flag
        println("reconstructed linear constraints, con x [start end]")
        display(idx_start_end)
        println("Og Assignment", assign_matrix)
        assign_matrix_bkup = copy(assign_matrix)
    end

    # splice!(conset.constraints, 2)
    
    for alcon_idx in 1:length(conset.constraints) # polytope alcon_idx
        alcon = conset.constraints[alcon_idx]
        if isa(alcon.con, TO.LinearConstraint{}) && !isempty(alcon.inds)
            # display(alcon.λ) #[collocation pts, plane]
            # display(alcon.λ)
            λmax = map(maximum, alcon.λ)#[collocation pts, 1]
            if debug_flag
                printstyled("Found Linear Cons \n"; color = :red)
                og_polytope_idx = findfirst(x -> x == alcon.con.A, lincons_A)
                println("Initial Allocation:", "poly:",og_polytope_idx, conset.constraints[alcon_idx].inds)
            end
            for collocation_idx_1_n in Iterators.reverse(eachindex(λmax)) #collocation pt
                max_λ = λmax[collocation_idx_1_n]
                if max_λ > 20.0
                    og_con_val = 0
                    collocation_idx = alcon.inds[collocation_idx_1_n]

                    dist_to_obs = zeros(length(lincons_A))
                    cur_state = state(alcon.Z[1][collocation_idx])
                    for i in 1:length(lincons_A)
                        if lincons_A[i] == alcon.con.A
                            dist_to_obs[i] = 1000#
                            og_con_val = maximum(lincons_A[i]*cur_state[1:3] .- lincons_b[i])
                        else
                            dist_to_obs[i] = maximum(lincons_A[i]*cur_state[1:3] .- lincons_b[i])
                        end
                    end
                    
                    min_val, min_polytope_idx_1_n = findmin(dist_to_obs)#find the most comfortable constraint
                    min_polytope_idx = lincon_idx[min_polytope_idx_1_n]#find the index in the constraint set
                    if debug_flag
                        printstyled("knot", collocation_idx," lambda ", max_λ, "poly:", og_polytope_idx, "\n"; color = :green)
                        println(dist_to_obs)
                        println("violation = " ,min_val, " at polytope ", min_polytope_idx_1_n)
                    end
                    if min_val < og_con_val #if the constraint is violated a little bit
                        # requiring it to be bigger than 0 is unreasonable, they are not going to be satisfied
                        if debug_flag
                            println("Found replacing constraint")
                        end
                        printstyled("adding  polytope ",min_polytope_idx_1_n, " to collocation ", collocation_idx, "\n", color = :blue)
                            # stagecon = TO.LinearConstraint(13,4,lincons_A[min_idx], lincons_b[min_idx],Inequality(),1:3)
                        # conset.constraints[alcon_idx] = ALConstraint{Float64}(alcon.Z[1], stagecon, alcon.inds, zeros(100)) 
                        assign_matrix[collocation_idx] = min_polytope_idx_1_n #this is different, for printing
                        #add new constraint inds and remove from current constraint inds
                        # println(conset.constraints[min_polytope_idx].nx)
                        # println(conset.constraints[min_polytope_idx].nu)
                        push!(conset.constraints[min_polytope_idx].nx,conset.constraints[min_polytope_idx].nx[end])
                        push!(conset.constraints[min_polytope_idx].nu,conset.constraints[min_polytope_idx].nu[end])
                        push!(conset.constraints[min_polytope_idx].inds, collocation_idx)#add constraint to 
                        push!(conset.constraints[min_polytope_idx].vals, zero(conset.constraints[min_polytope_idx].vals[end]))
                        # display(conset.constraints[min_idx].jac[end])
                        # println("p",p,"w",w)
                        push!(conset.constraints[min_polytope_idx].jac,  zero(conset.constraints[min_polytope_idx].jac[end]))
                        push!(conset.constraints[min_polytope_idx].jac_scaled, zero(conset.constraints[min_polytope_idx].jac_scaled[end]))
                        push!(conset.constraints[min_polytope_idx].λ, zero(conset.constraints[min_polytope_idx].λ[end]))
                        push!(conset.constraints[min_polytope_idx].μ, zero(conset.constraints[min_polytope_idx].μ[end]).+1)
                        push!(conset.constraints[min_polytope_idx].μinv, zero(conset.constraints[min_polytope_idx].μinv[end]).+1)
                        push!(conset.constraints[min_polytope_idx].λbar, zero(conset.constraints[min_polytope_idx].λbar[end]))
                        push!(conset.constraints[min_polytope_idx].λproj, zero(conset.constraints[min_polytope_idx].λproj[end]))
                        push!(conset.constraints[min_polytope_idx].λscaled, zero(conset.constraints[min_polytope_idx].λscaled[end]))
                        push!(conset.constraints[min_polytope_idx].viol, zero(conset.constraints[min_polytope_idx].viol[end]))
                        #this we can calculate
                        push!(conset.constraints[min_polytope_idx].∇proj, zero(conset.constraints[min_polytope_idx].∇proj[end]))
                        push!(conset.constraints[min_polytope_idx].∇proj_scaled, zero(conset.constraints[min_polytope_idx].∇proj_scaled[end]))
                        push!(conset.constraints[min_polytope_idx].∇²proj, zero(conset.constraints[min_polytope_idx].∇²proj[end]))
                        push!(conset.constraints[min_polytope_idx].grad, zero(conset.constraints[min_polytope_idx].grad[end]))
                        push!(conset.constraints[min_polytope_idx].hess, zero(conset.constraints[min_polytope_idx].hess[end]))
                        # using lagrange multipiler of a member of the new constraint, this may be not good, try zeros first

                        @assert conset.constraints[alcon_idx].inds[collocation_idx_1_n] == collocation_idx
                        # push!(marked_for_deletion, collocation_idx_1_n)

                        deleteat!(conset.constraints[alcon_idx].nx, collocation_idx_1_n)
                        deleteat!(conset.constraints[alcon_idx].nu, collocation_idx_1_n)
                        deleteat!(conset.constraints[alcon_idx].inds, collocation_idx_1_n)
                        deleteat!(conset.constraints[alcon_idx].vals, collocation_idx_1_n)
                        deleteat!(conset.constraints[alcon_idx].jac, collocation_idx_1_n)
                        deleteat!(conset.constraints[alcon_idx].jac_scaled, collocation_idx_1_n)
                        deleteat!(conset.constraints[alcon_idx].λ, collocation_idx_1_n)
                        deleteat!(conset.constraints[alcon_idx].μ, collocation_idx_1_n)
                        deleteat!(conset.constraints[alcon_idx].μinv, collocation_idx_1_n)
                        deleteat!(conset.constraints[alcon_idx].λbar, collocation_idx_1_n)
                        deleteat!(conset.constraints[alcon_idx].λproj, collocation_idx_1_n)
                        deleteat!(conset.constraints[alcon_idx].λscaled, collocation_idx_1_n)
                        deleteat!(conset.constraints[alcon_idx].viol, collocation_idx_1_n)
                        deleteat!(conset.constraints[alcon_idx].∇proj, collocation_idx_1_n)

                        deleteat!(conset.constraints[alcon_idx].∇proj_scaled, collocation_idx_1_n)
                        deleteat!(conset.constraints[alcon_idx].∇²proj, collocation_idx_1_n)
                        deleteat!(conset.constraints[alcon_idx].grad, collocation_idx_1_n)
                        deleteat!(conset.constraints[alcon_idx].hess, collocation_idx_1_n)


                        if debug_flag
                            printstyled("removing  polytope ",og_polytope_idx, "on collocation ", collocation_idx, "\n",color = :blue)
                            println("Change Allocation:", "box:", min_polytope_idx_1_n, conset.constraints[min_polytope_idx].inds)
                        end

                    else
                        if debug_flag
                            println("No other constraint is better")
                        end
                    end
                else

                end
            end
            if debug_flag
                println("Final Allocation:", "box:", og_polytope_idx, conset.constraints[alcon_idx].inds)
                println("=====================================")
            end
        end
        
    end
    if debug_flag
        println("Diff Assignment", assign_matrix - assign_matrix_bkup)
        println("New Assignment", assign_matrix)
    end
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