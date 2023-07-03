#!/usr/bin/env julia
# run(`source /home/yifei/ws/devel/setup.bash`)
import Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()
using RobotOS
@rosimport kr_planning_msgs.msg: TrajectoryDiscretized, SplineTrajectory
@rosimport visualization_msgs.msg: Marker
@rosimport decomp_ros_msgs.msg: Polyhedron, PolyhedronArray
@rosimport geometry_msgs.msg: Point
@rosimport std_msgs.msg: ColorRGBA, Float64, Int32
rostypegen()
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
using ForwardDiff
using FiniteDiff
using Random
using SparseArrays
const RD = RobotDynamics
const TO = TrajectoryOptimization
const Dynamics = RobotDynamics

using Altro: ALTROSolver
using Rotations
using Plots



include("./yifei/quadrotor_kr.jl")

obstacles = []# Ref{PolyhedronArray}() # Store it into a `Ref` which it's mutable but has a constant type
# const traj_received = Ref{SplineTrajectory}()

macro debugtask(ex)
    quote
      try
        $(esc(ex))
      catch e
        bt = stacktrace(catch_backtrace())
        io = IOBuffer()
        showerror(io, e, bt)
        errstr = String(take!(io))
        RobotOS.logfatal("Error: $errstr")
        exit()
      end
    end
  end

# logging: loginfo("$(RobotOS.get_caller_id()) I heard $(msg.data)")
function path_callback(msg::TrajectoryDiscretized, pub_obj::Publisher{Marker})
    #ToDo: compile to obj save startup time https://stackoverflow.com/questions/73599900/julia-seems-to-be-very-slow
    @debugtask begin
        println("path received")
        global obstacles
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

            solver = ALTROSolver(Quadrotor_kr(traj_ref = wpts, vel_ref = vel, FM_ref = FM, obstacles = obstacles, t_vec = msg.t)..., verbose=0)
            Z0 = deepcopy(get_trajectory(solver))
            TO.initial_trajectory!(solver,Z0)
            solve!(solver)
            X = states(solver)     # alternatively states(prob)
            U = controls(solver)   # alternatively controls(prob)
            ilqr = Altro.get_ilqr(solver)
            T = gettimes(solver)
            K = ilqr.K  # feedback gain matrices
            d = ilqr.d  # feedforward gains. Should be small.

            plotting_bool = true
            if plotting_bool
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

                u1 = [u[1] for u in U]
                u2 = [u[2] for u in U]
                u3 = [u[3] for u in U]
                u4 = [u[4] for u in U]

                u1_ref = [u[1] for u in U_hover]
                u2_ref = [u[2] for u in U_hover]
                u3_ref = [u[3] for u in U_hover]
                u4_ref = [u[4] for u in U_hover]


                p = plot(T,x_plot, label="refined path x")
                # p2 = plot(U[:,1], label="refined ctrl")
                plot!(T,y_plot, label="refined path y")
                plot!(T,z_plot, label="refined path z")
                plot!(msg.t,x_ref_plot,label="original path x")
                plot!(msg.t,y_ref_plot,label="original path y")
                plot!(msg.t,z_ref_plot,label="original path z")
                # legend!()

                p2 = plot(T,vx_plot, label="refined vx")
                plot!(T,vy_plot, label="refined vy")
                plot!(T,vz_plot, label="refined vz")
                plot!(msg.t,vx_ref_plot,label="ref vx")
                plot!(msg.t,vy_ref_plot,label="ref vy")
                plot!(msg.t,vz_ref_plot,label="ref vz")
                # legend!()

                p3 = plot(T[1:end-1],u1, label="u1")
                plot!(T[1:end-1],u2, label="u2")
                plot!(T[1:end-1],u3, label="u3")
                plot!(T[1:end-1],u4, label="u4")
                plot!(msg.t,u1_ref,label="u1 ref")
                plot!(msg.t,u2_ref,label="u2 ref")
                plot!(msg.t,u3_ref,label="u3 ref")
                plot!(msg.t,u4_ref,label="u4 ref")

                savefig(p, "refined_path.png")
                savefig(p2, "refined_vel.png")
                savefig(p3, "control_inputs.png")
            end

            send_msg = Marker()
            send_msg.type = 4
            send_msg.header = msg.header
            send_msg.scale.x = 0.2
            
            send_msg.points = Point[]
            for x in X
                push!(send_msg.points, Point(x[1], x[2], x[3]))
            end
            send_msg.color = ColorRGBA(1.0, 0.0, 0.0, 1.0)
            publish(pub_obj, send_msg)
        end
    end
    # println(x)
    # pt_msg = Point(msg.x, msg.y, 0.0)
    # 
end

function ros_to_polyhedron(msg::Polyhedron)
    @debugtask begin
        poly = Matrix{Float64}(undef,6,0)
        numPoints = length(msg.points)

        for i in 1:numPoints
            pt =  [msg.points[i].x,
                        msg.points[i].y,
                        msg.points[i].z]
            n =  [msg.normals[i].x,
                        msg.normals[i].y,
                        msg.normals[i].z]
            ptn = vcat(pt, n)
            poly = hcat(poly, ptn)
            # print(size(poly))
        end
        A,b = pt_normal_to_halfspace(poly,numPoints)
        return A, b
    end
end

function obs_callback(msg::PolyhedronArray)
    @debugtask begin
        polys = []

        for polyhedron in msg.polyhedrons
            push!(polys, ros_to_polyhedron(polyhedron))
        end
        global obstacles
        obstacles = polys
        # println(polys)
    end
end

function pt_normal_to_halfspace(poly, size_n)
    @debugtask begin
        A = zeros(size_n, 3)
        b = zeros(size_n)
        p0 = mean(poly[1:3,:], dims=2)
        for ii in 1:size_n
            p = poly[1:3,ii]
            n = poly[4:6,ii]
            c = dot(n, p)
            # if dot(n, p0) - c > 0 # this section seems to make mistakes.
            #     n = -n
            #     c = -c
            # end
            n = vec(n)
            A[ii,:] = n
            b[ii] = c
        end
        return A,b 
    end
end

# function spline_callback(msg::SplineTrajectory)
#     @debugtask begin
#         println("traj updated")
#         traj_received[] = msg
#         # println(traj_received.data[1].t_total)
#     end
# end


# function foo_callback(msg)
#   @debugtask begin
#      error("whoops")
#   end
# end
# function callback(msg::Pose2D, pub_obj::Publisher{Point})
#     pt_msg = Point(msg.x, msg.y, 0.0)
#     publish(pub_obj, pt_msg)
# end

# function loop(pub_obj)
#     loop_rate = Rate(5.0)%
#     while ! is_shutdown()
#         # npt = Point(rand(), rand(), 0.0)
#         # publish(pub_obj, npt)
#         rossleep(loop_rate)
#     end
# end

function main()
    init_node("path_refiner")
    pub = Publisher{Marker}("iLQR_path", queue_size=1)
    # sub = Subscriber{Pose2D}("pose", callback, (pub,), queue_size=10)
    
    sub1 = Subscriber{TrajectoryDiscretized}("/spline_traj_samples", path_callback, (pub,), queue_size=1)
    sub2 = Subscriber{PolyhedronArray}("/quadrotor/local_plan_server/trajectory_planner/sikangpolyhedron", obs_callback; queue_size=1)
    # sub3 = Subscriber{SplineTrajectory}("/quadrotor/local_plan_server/trajectory", spline_callback; queue_size=1)
    # sub3 = Subscriber{Marker} # PlanTwoPointActionGoal  ## we want this to have initial velocity and attitude but it doesnt have it
    spin()
end

if ! isinteractive()
    main()
end