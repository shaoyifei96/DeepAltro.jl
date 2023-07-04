function Quadrotor_kr(Rot=UnitQuaternion{Float64}; traj_ref, vel_ref, FM_ref, obstacles, t_vec, ini_rot, ini_w,
        costfun=:QuadraticCost, normcon=false)
        costfun == :QuatLQR ? sq = 0 : sq = 1
        model = RobotZoo.Quadrotor{Rot}()
        n,m = RD.dims(model)
        kf, km = model.kf, model.km
        L = model.motor_dist
        fM = @SMatrix [
            kf   kf   kf   kf;
            0    L*kf 0   -L*kf;
        -L*kf 0    L*kf 0;
            km  -km   km  -km;]
        MM = inv(fM)

        opts = SolverOptions(
            penalty_scaling=100.,
            penalty_initial=1.0,
            projected_newton=false,
            iterations_outer = 60
        )

        # discretization
        N = length(traj_ref) # number of knot points
        u0 = @SVector fill(0.5*9.81/4, m) #TODO: Change to actual vehicle mass

        tf = convert(Float64, t_vec[end]) # Assigned from SplineTrajectory segment 0 total_time
        println("tf: ", tf)
        # what is a reasonable longest final time!?
        dt = tf/(N-1) # total time

        # Initial condition
        x0_pos = traj_ref[1] #this is planning in local coordinates
        # x0 = RobotDynamics.build_state(model, x0_pos, ini_rot, vel_ref[1], ini_w)
        x0 = RobotDynamics.build_state(model, x0_pos,UnitQuaternion(I) , vel_ref[1], zeros(3))

        # cost
        # costfun == :QuatLQR ? sq = 0 : sq = 1
        # rm_quat = @SVector [1,2,3,4,5,6,8,9,10,11,12,13]                                    #       x         q:1e-5*sq      v 1e-3    w
        # Q = Diagonal(Q_diag)
        R = Diagonal(@SVector fill(0.1,m)) 
        q_nom = UnitQuaternion(I)
        ω_nom = zeros(3)
        # x_nom = Dynamics.build_state(model, zeros(3), q_nom, v_nom, ω_nom)
        #why is the state zero? should it follow referecne trajectory?
        
        # if costfun == :QuatLQR
        #     cost_nom = QuatLQRCost(Q*dt, R*dt, x_nom, w=0.0)
        # elseif costfun == :ErrorQuad
        #     cost_nom = ErrorQuadratic(model, Diagonal(Q_diag[rm_quat])*dt, R*dt, x_nom)
        # else
        #     cost_nom = LQRCost(Q*dt, R*dt, x_nom, u0)
        # end
        U_ref = [MM * v for v in FM_ref] # inital reference input control
        U_hover = [copy(u0) for k = 1:N-1]

        traj = traj_ref #about 50/7 = 7 points
        # how many waypoints do you leave behind
        # waypoints
        times = 1:N #round.(Int, range(1, stop=101, length=length(traj)))
        println("length of traj is $(length(traj))")
        # times = [33, 66, 101]
        #works well to reduce rotation a little 
        Qw_diag = Dynamics.fill_state(model, 0.03, 0,0, 0.01) #no waypoint cost since only the final point matters
        # Qw_diag = Dynamics.fill_state(model, 10, 0.0,0,0) #to check correctness of things
                                        #    x   q:1*sq v:1 w:1
        Qf_diag = Dynamics.fill_state(model, 1., 0.0, 1, 0.0)
        xf = Dynamics.build_state(model, traj[end], UnitQuaternion(I), vel_ref[end], zeros(3))

        costs = map(1:length(traj)) do i
            xg = Dynamics.build_state(model, traj[i], q_nom, vel_ref[i], ω_nom)
            if times[i] == N
                Q = Diagonal(Qf_diag)
                # w = 4.0
            else
                Q = Diagonal(Qw_diag) * dt
                # w = 0.1
            end
            # if costfun == :QuatLQR
            #     QuatLQRCost(Q, R*dt, xg, w=w)
            # elseif costfun == :ErrorQuad
            #     Qd = diag(Q)
            #     ErrorQuadratic(model, Diagonal(Qd[rm_quat]), R, xg)
            # else

            # LQRCost(Q, R, xg, U_hover[i]) #cost wrt ref input or hover input?
            if costfun == :QuatLQR
                TO.QuatLQRCost(Q, R, xg, u0)
            else
                LQRCost(Q, R, xg, u0)
            end
            
            # end
        end

        # println(typeof(costs_all))
        # println(typeof(costs_all[1]))
        # println(typeof(costs[1]))
        # println(typeof(cost_nom))
        
        obj = Objective(costs)

        # Initialization
    
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



        #need to build constraints: depending on which position is in polytope, move to next one 
        #if 1 in poly, 2 in poly, 3 not in poly, then assign poly constraint to 1 and 2
        polytope_counter = 1
        traj_counter = 1
        
        #make a double for loop to check if each point is in each polytope
        debug_flag = false
        if debug_flag 
            result_mat = zeros(Int8, length(traj), length(obstacles))
            for i = 1:length(traj)
                for j = 1:length(obstacles)
                    result_mat[i,j] = all(obstacles[j][1]*traj[i]-obstacles[j][2] .<= 0.0)
                end
            end
            display(result_mat)
            println("End Matrix")
            # exit()
        end
        # while traj_counter <= length(traj) && polytope_counter <= length(obstacles)
        #     poly = obstacles[polytope_counter]
        #     if  all(poly[1]*traj[traj_counter]-poly[2] .<= 0.0)
        #         println("""Adding constraint, point $traj_counter is in polytope $polytope_counter""")
        #         add_constraint!(conSet, LinearConstraint(n,m,-poly[1], -poly[2],Inequality(),1:3), traj_counter)
        #         traj_counter += 1
        #     else
        #         polytope_counter += 1
        #     end
        # end
        # if traj_counter < length(traj)
        #     println("Warning: not all traj initialized points are in polytopes")
        # end
        conSettemp = ConstraintList(n,m,N)
        start_end_mat = zeros(Int8,0,2)
        polytope_counter = 1
        polytope_reverse_vec = zeros(Bool, length(obstacles))
        exit_flag = false
        while !exit_flag && polytope_counter <= length(obstacles)
            obs_start = 0
            obs_end   = 0
            poly_res_start = ones(Bool, length(obstacles[polytope_counter][2])) #false if point not in polytope
            poly_res_end   = ones(Bool, length(obstacles[polytope_counter][2]))
            poly = obstacles[polytope_counter]
            if polytope_reverse_vec[polytope_counter]
                poly= (-poly[1], -poly[2])
            end

            for traj_counter = 1:length(traj)
                poly_res = poly[1]*traj[traj_counter]-poly[2] .<= 0.0
                if all(poly_res) && obs_start == 0
                    obs_start = traj_counter
                    if traj_counter > 1
                        poly_res_start = poly[1]*traj[traj_counter-1]-poly[2] .<= 0.0
                    end
                    continue
                    # this line shows when traj first enters the polytope 
                    # what constraints are satisfied for the previous point
                elseif all(poly_res) && obs_start != 0
                    if traj_counter == length(traj)
                        obs_end = length(traj)
                        exit_flag = true
                    end
                     #inside polytope
                elseif !all(poly_res) && obs_start != 0
                    poly_res_end = poly[1]*traj[traj_counter]-poly[2] .<= 0.0
                    obs_end = traj_counter-1
                    break
                elseif !all(poly_res) && obs_start == 0
                    continue #have not entered polytope yet
                end
                #if inside polytope all the way to end
                
            end
            println("polytope $polytope_counter")
            println(poly_res_start, poly_res_end)
            poly_res_all = poly_res_start .& poly_res_end
            println(poly_res_all)
            if obs_start == 0
                println("Warning: polytope $polytope_counter reverse = $(polytope_reverse_vec[polytope_counter]) not in trajectory")
                if polytope_reverse_vec[polytope_counter]
                    println("Error: traj not in polytope")
                    polytope_counter += 1
                    if polytope_counter > length(obstacles)
                        break
                    end
                end
                polytope_reverse_vec[polytope_counter] = true
                continue
            end
            start_end_mat = [start_end_mat; [obs_start, obs_end]']
            newcon = LinearConstraint(n,m,poly[1][poly_res_all,:], poly[2][poly_res_all],Inequality(),1:3) 
            add_constraint!(conSettemp, newcon,1)
            polytope_counter += 1
        end
        display(start_end_mat)
        println(polytope_reverse_vec)
        for i = 1:length(conSettemp)
            if i == 1
                start_idx = start_end_mat[i,1]
            else 
                if start_end_mat[i-1,2] > start_end_mat[i,1]
                    start_idx = (start_end_mat[i-1,2] + start_end_mat[i,1]) ÷ 2
                else
                    start_idx = start_end_mat[i,1]
                end
            end

            if i == length(conSettemp)
                last_idx = start_end_mat[i,2]
            else
                if start_end_mat[i,2] > start_end_mat[i+1,1]
                    last_idx = (start_end_mat[i,2] + start_end_mat[i+1,1]) ÷ 2
                else
                    last_idx = start_end_mat[i,2]
                end
            end
            add_constraint!(conSet, conSettemp[i], start_idx:last_idx)
            println("adding constraint $i from $start_idx to $last_idx")
        end


        ## this function does not take into account of two sided walls
        # polytope_counter = length(obstacles)
        # traj_counter = length(traj)
        # last_idx = length(traj)
        # while traj_counter >= 0 && polytope_counter >= 1
        #     poly = obstacles[polytope_counter]
        #     if traj_counter == 0
        #         poly_res = zeros(Bool, length(poly[2]))
        #     else
        #         poly_res = poly[1]*traj[traj_counter]-poly[2] .<= 0.0
        #     end
        #     if  all(poly_res)
        #         traj_counter -= 1
        #     else
        #         if traj_counter >= last_idx
        #             polytope_counter -= 1
        #         else
        #             newcon = LinearConstraint(n,m,-poly[1][.!poly_res,:], -poly[2][.!poly_res],Inequality(),1:3) 
        #             add_constraint!(conSet, newcon, (traj_counter+1):last_idx)
        #             println("""Adding constraint $(traj_counter+1) to $last_idx is in polytope $polytope_counter""")
        #             polytope_counter -= 1
        #             last_idx = traj_counter
        #         end
        #     end
        # end
        # if traj_counter < last_idx
        #     newcon = LinearConstraint(n,m,-poly[1][.!poly_res,:], -poly[2][.!poly_res],Inequality(),1:3) 
        #     add_constraint!(conSet, newcon, 1:last_idx)
        #     println("""Adding constraint 1 to $last_idx is in polytope $polytope_counter""")
        # end
        
# exit()
        # Problem
        prob = Problem(model, obj, x0, tf, xf=xf, constraints=conSet)
        initial_controls!(prob, U_hover)
        rollout!(prob)

        return prob, opts
end
