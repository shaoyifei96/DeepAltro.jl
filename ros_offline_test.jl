#!!!need to run ros_offline_prep first if first time running!!!
if ! isinteractive()
    include("./ros_offline_prep.jl")
end

using .kr_planning_msgs.msg: TrajectoryDiscretized, SplineTrajectory
using .visualization_msgs.msg
using .decomp_ros_msgs.msg
using .geometry_msgs.msg
using .std_msgs.msg: ColorRGBA, Float64Msg, Int32Msg

using Printf
using StaticArrays
using Statistics

using Test
using Altro
using BenchmarkTools
using TrajectoryOptimization
using RobotDynamics
using RobotZoo
using StaticArrays, LinearAlgebra
using JLD2
# using FileIO
using ForwardDiff
using FiniteDiff
using Random
using SparseArrays
using Revise
const RD = RobotDynamics
const TO = TrajectoryOptimization
const Dynamics = RobotDynamics

using Altro: ALTROSolver
using Rotations
using Plots

using Profile
using CSV
using DataFrames
online_mode = false

include("./yifei/quadrotor_kr.jl")

wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)
function problem_solving(msg, obstacles)
    println("path received")
    if isempty(obstacles) # in master branch we have no obstacles as polytopes
        println("no obstacles")
    else
        # println(fieldnames(traj))
        wpts =Vector{SVector{3, Float64}}([])
        vel =Vector{SVector{3, Float64}}([])
        FM =Vector{SVector{4, Float64}}([])
        #acc is also avaiable but not used
        for point in msg.pos
            push!(wpts, @SVector[point.x, point.y, point.z])
        end #This correctly convert messages to SVector

        for point in msg.vel
            push!(vel, @SVector[point.x, point.y, point.z])
        end 
        for i in 1:length(msg.moment)
            push!(FM, @SVector[msg.thrust[i], msg.moment[i].x, msg.moment[i].y, msg.moment[i].z])
        end 
        #debug: TODO: remove
        Rot=UnitQuaternion{Float64}
        ini_rot = UnitQuaternion{Float64}(msg.inital_attitude.w, msg.inital_attitude.x, msg.inital_attitude.y, msg.inital_attitude.z)
        ini_w   = @SVector [msg.initial_omega.x, msg.initial_omega.y, msg.initial_omega.z]
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
        U_hover = [MM * v for v in FM]
        # end of debug

        solver = ALTROSolver(Quadrotor_kr(traj_ref = wpts, vel_ref = vel, FM_ref = FM, 
        obstacles = obstacles, t_vec = msg.t, ini_rot = ini_rot,ini_w = ini_w)..., verbose=0)
        Z0 = deepcopy(get_trajectory(solver))
        TO.initial_trajectory!(solver,Z0)
        solve!(solver)
        X = states(solver)     # alternatively states(prob)
        U = controls(solver)   # alternatively controls(prob)
        ilqr = Altro.get_ilqr(solver)
        T = gettimes(solver)
        K = ilqr.K  # feedback gain matrices
        d = ilqr.d  # feedforward gains. Should be small.

        stats = Altro.stats(solver) 
        # println(fieldnames(typeof()))    
        # println()
        x_plot = [v[1] for v in X]
        y_plot = [v[2] for v in X]
        z_plot = [v[3] for v in X]
        vx_plot = [v[8] for v in X]
        vy_plot = [v[9] for v in X]
        vz_plot = [v[10] for v in X]

        x_ref_plot = [v[1] for v in wpts]
        y_ref_plot = [v[2] for v in wpts]
        z_ref_plot = [v[3] for v in wpts]

        vx_ref_plot = [v[1] for v in vel]
        vy_ref_plot = [v[2] for v in vel]
        vz_ref_plot = [v[3] for v in vel]

        u0 = 0.5*9.81/4
        u1 = [u[1].-u0 for u in U]
        u2 = [u[2].-u0 for u in U]
        u3 = [u[3].-u0 for u in U]
        u4 = [u[4].-u0 for u in U]

        u1_ref = [u[1].-u0 for u in U_hover]
        u2_ref = [u[2].-u0 for u in U_hover]
        u3_ref = [u[3].-u0 for u in U_hover]
        u4_ref = [u[4].-u0 for u in U_hover]
        plotting_bool = true
        if plotting_bool
            p = plot(T,x_plot, label="refined path x", lc=:red,ls=:solid)
            # p2 = plot(U[:,1], label="refined ctrl")
            plot!(T,y_plot, label="refined path y",lc=:green,ls=:solid)
            plot!(T,z_plot, label="refined path z",lc=:blue,ls=:solid)
            plot!(msg.t,x_ref_plot,label="original path x",lc=:red,ls=:dash)
            plot!(msg.t,y_ref_plot,label="original path y",lc=:green,ls=:dash)
            plot!(msg.t,z_ref_plot,label="original path z",lc=:blue,ls=:dash)
            # legend!()

            p2 = plot(T,vx_plot, label="refined vx",lc=:red,ls=:solid)
            plot!(T,vy_plot, label="refined vy",lc=:green,ls=:solid)
            plot!(T,vz_plot, label="refined vz",lc=:blue,ls=:solid)
            plot!(msg.t,vx_ref_plot,label="ref vx",lc=:red,ls=:dash)
            plot!(msg.t,vy_ref_plot,label="ref vy",lc=:green,ls=:dash)
            plot!(msg.t,vz_ref_plot,label="ref vz",lc=:blue,ls=:dash)
            # legend!()

            p3 = plot(T[1:end-1],u1, label="u1",lc=:red,ls=:solid)
            plot!(T[1:end-1],u2, label="u2",lc=:green,ls=:solid)
            plot!(T[1:end-1],u3, label="u3",lc=:blue,ls=:solid)
            plot!(T[1:end-1],u4, label="u4",lc=:black,ls=:solid)
            plot!(msg.t,u1_ref,label="u1 ref",lc=:red,ls=:dash)
            plot!(msg.t,u2_ref,label="u2 ref",lc=:green,ls=:dash)
            plot!(msg.t,u3_ref,label="u3 ref",lc=:blue,ls=:dash)
            plot!(msg.t,u4_ref,label="u4 ref",lc=:black,ls=:dash)


            p4 = plot(x_plot, y_plot, z_plot, label="refined path",lc=:green,ls=:solid, aspect_ratio = :equal)
            plot!(x_ref_plot, y_ref_plot, z_ref_plot, label="original path",lc=:black,ls=:dash)
            
            savefig(p, "refined_path.png")
            savefig(p2, "refined_vel.png")
            savefig(p3, "control_inputs.png")
            savefig(p4, "refined_path_3d.png")

        end
    end
    dt = msg.t[end] / 99.0
    # return controls, total ilqr iterations,   TODO: total time may be not correct when doing time optimization
    return X, stats.cost[end], (norm(diff(u1))+norm(diff(u2))+norm(diff(u3))+norm(diff(u4)))/dt, norm(u1)+norm(u2)+norm(u3)+norm(u4) , stats.iterations, Integer(stats.status), T[end], stats.tsolve
end

problem_id_v = zeros(Int64, 0) 
solve_time_v = zeros(Float64, 0)
u_smooth_v = zeros(Float64, 0)
u_norm_v = zeros(Float64, 0)
traj_time_v = zeros(Float64, 0)
solver_iter_v = zeros(Int64, 0)
solver_status_v = zeros(Int64, 0)
traj_cost_v = zeros(Float64, 0)
# traj jerk cost = zeros(Float64, 0)
if online_mode
    pub = Publisher{Marker}("iLQR_path", queue_size=1)
    init_node("path_refiner")
end
for i = 30:129
    println("problem $i ")
    data_dict = load("/home/yifei/Documents/optimal_ctrl/Altro_og/Altro.jl/data/problems_starting_from_yuwei/problem$i.jld2")# msg obstacles total_iter
    msg = data_dict["msg"]
    obstacles = data_dict["obstacles"]
    total_iter = data_dict["total_iter"]
    
    X, traj_cost, u_smooth, u_norm, solver_iter, status, traj_time, solve_time = problem_solving(msg, obstacles)
    push!(u_smooth_v, u_smooth)
    push!(problem_id_v, i)
    push!(solve_time_v, solve_time)
    push!(u_norm_v, u_norm)
    push!(traj_time_v, traj_time)
    push!(solver_iter_v, solver_iter)
    push!(solver_status_v, status)
    push!(traj_cost_v,traj_cost)
    send_msg = Marker()
    send_msg.type = 4
    send_msg.header = msg.header
    send_msg.scale.x = 0.2
    
    send_msg.points = Point[]
    for x in X
        push!(send_msg.points, Point(x[1], x[2], x[3]))
    end
    send_msg.color = ColorRGBA(1.0, 0.0, 0.0, 1.0)
    if online_mode
        publish(pub, send_msg) 
        # msg = obj.msg 
        # obstacles = obj.obstacles
        #also avaiable: total_iter
        spin()
        wait_for_key("press any key to continue")
    end
    
end
df = DataFrame(
problem_id = problem_id_v,    
solve_time = solve_time_v, 
u_norm_over_hover = u_norm_v,
traj_time = traj_time_v,
solver_iter = solver_iter_v,
solver_status = solver_status_v,
u_smooth = u_smooth_v,
traj_cost = traj_cost_v
               )
CSV.write("ros_offline_stats.csv", df)

