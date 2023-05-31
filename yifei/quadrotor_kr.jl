function Quadrotor_kr(Rot=UnitQuaternion{Float64}; traj_ref, obstacles, time_total,
        costfun=:Quadratic, normcon=false)
        model = RobotZoo.Quadrotor{Rot}()
        n,m = RD.dims(model)

        opts = SolverOptions(
            penalty_scaling=100.,
            penalty_initial=1.0,
            projected_newton=false,
        )

        # discretization
        N = 101 # number of knot points
        u0 = @SVector fill(0.5*9.81/4, m) #TODO: Change to actual vehicle mass

        tf = convert(Float64, time_total) # Assigned from SplineTrajectory segment 0 total_time
        # what is a reasonable longest final time!?
        dt = tf/(N-1) # total time

        # Initial condition
        x0_pos = traj_ref[1] #this is planning in local coordinates
        x0 = RobotDynamics.build_state(model, x0_pos, UnitQuaternion(I), zeros(3), zeros(3))

        # cost
        costfun == :QuatLQR ? sq = 0 : sq = 1
        rm_quat = @SVector [1,2,3,4,5,6,8,9,10,11,12,13]
        Q_diag = Dynamics.fill_state(model, 1e-5, 0.0, 0.0, 0.0)
                                    #       x         q:1e-5*sq      v 1e-3    w
        Q = Diagonal(Q_diag)
        R = Diagonal(@SVector fill(1e-2,m))
        q_nom = UnitQuaternion(I)
        v_nom, ω_nom = zeros(3), zeros(3)
        x_nom = Dynamics.build_state(model, zeros(3), q_nom, v_nom, ω_nom)
        #why is the state zero? should it follow referecne trajectory?
        
        if costfun == :QuatLQR
            cost_nom = QuatLQRCost(Q*dt, R*dt, x_nom, w=0.0)
        elseif costfun == :ErrorQuad
            cost_nom = ErrorQuadratic(model, Diagonal(Q_diag[rm_quat])*dt, R*dt, x_nom)
        else
            cost_nom = LQRCost(Q*dt, R*dt, x_nom, u0)
        end

        traj = traj_ref[1:7:end] #about 50/7 = 7 points
        # how many waypoints do you leave behind
        # waypoints
        times = round.(Int, range(1, stop=101, length=length(traj)))
        println("length of traj is $(length(traj))")
        # times = [33, 66, 101]
        Qw_diag = Dynamics.fill_state(model, 1e3, 0.0,0.0,0.0)
                                        #    x   q:1*sq v:1 w:1
        Qf_diag = Dynamics.fill_state(model, 10., 100*sq, 10, 10)
        xf = Dynamics.build_state(model, traj[end], UnitQuaternion(I), zeros(3), zeros(3))

        costs = map(1:length(traj)) do i
            r = traj[i]
            xg = Dynamics.build_state(model, r, q_nom, v_nom, ω_nom)
            if times[i] == N
                Q = Diagonal(Qf_diag)
                w = 40.0
            else
                Q = Diagonal(1e-3*Qw_diag) * dt
                w = 0.1
            end
            if costfun == :QuatLQR
                QuatLQRCost(Q, R*dt, xg, w=w)
            elseif costfun == :ErrorQuad
                Qd = diag(Q)
                ErrorQuadratic(model, Diagonal(Qd[rm_quat]), R, xg)
            else
                LQRCost(Q, R, xg, u0)
            end
        end
        # println(times)
        costs_all = map(1:N) do k
            # println("k is $k")
            i = findfirst(x->(x ≥ k), times)
            # println(i)
            if k ∈ times
                costs[i]
                # println("Using goodcost")
            else
                cost_nom
                
            end
        end

        # println(typeof(costs_all))
        # println(typeof(costs_all[1]))
        # println(typeof(costs[1]))
        # println(typeof(cost_nom))
        
        obj = Objective(costs_all)

        # Initialization
        U_hover = [copy(u0) for k = 1:N-1] # initial hovering control trajectory

        # Constraints
        conSet = ConstraintList(n,m,N)
        if normcon
            if use_rot == :slack
                add_constraint!(conSet, QuatSlackConstraint(), 1:N-1)
            else
                add_constraint!(conSet, QuatNormConstraint(), 1:N-1)
                u0 = [u0; (@SVector [1.])]
            end
        end
        bnd = BoundConstraint(n,m, u_min=0.0, u_max=12.0)
        add_constraint!(conSet, bnd, 1:N-1)

        # Problem
        prob = Problem(model, obj, x0, tf, xf=xf, constraints=conSet)
        initial_controls!(prob, U_hover)
        rollout!(prob)

        return prob, opts
end
