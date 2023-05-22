#!/usr/bin/env julia
# run(`source /home/yifei/ws/devel/setup.bash`)
import Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()
using RobotOS
@rosimport planning_ros_msgs.msg: PlanTwoPointActionGoal
@rosimport visualization_msgs.msg: Marker
@rosimport decomp_ros_msgs.msg: PolyhedronArray
@rosimport geometry_msgs.msg: Point
@rosimport std_msgs.msg: ColorRGBA
rostypegen()
using .planning_ros_msgs.msg
using .visualization_msgs.msg
using .decomp_ros_msgs.msg
using .geometry_msgs.msg
using .std_msgs.msg

using Printf
using StaticArrays

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

include("./yifei/quadrotor_kr.jl")

const obstacles = Ref{PolyhedronArray}() # Store it into a `Ref` which it's mutable but has a constant type

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
function path_callback(msg::Marker, pub_obj::Publisher{Marker})
    #ToDo: compile to obj save startup time https://stackoverflow.com/questions/73599900/julia-seems-to-be-very-slow
    @debugtask begin
        println("path received")
        if isempty(obstacles)
            println("obstacles not yet assigned, do nothing")
        else 
            wpts =Vector{SVector{3, Float64}}([])
            for point in msg.points
                push!(wpts, @SVector[point.x, point.y, point.z])
            end #This correctly convert messages to SVector
            solver = ALTROSolver(Quadrotor_kr(traj = wpts, obstacles = obstacles)..., verbose=0)
            Z0 = deepcopy(get_trajectory(solver))
            TO.initial_trajectory!(solver,Z0)
            solve!(solver)
            X = states(solver)     # alternatively states(prob)
            U = controls(solver)   # alternatively controls(prob)
            ilqr = Altro.get_ilqr(solver)
            K = ilqr.K  # feedback gain matrices
            d = ilqr.d  # feedforward gains. Should be small.

            send_msg = deepcopy(msg)
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

function obs_callback(msg::PolyhedronArray)
    @debugtask begin
        println("obstacles updated")
        obstacles = msg
    end
end



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
    pub = Publisher{Marker}("iLQR_path", queue_size=10)
    # sub = Subscriber{Pose2D}("pose", callback, (pub,), queue_size=10)
    
    sub1 = Subscriber{Marker}("/quadrotor/local_plan_server/trajectory_planner/kino_astar_list", path_callback, (pub,), queue_size=1)
    sub2 = Subscriber{PolyhedronArray}("/quadrotor/local_plan_server/trajectory_planner/sikangpolyhedron", obs_callback; queue_size=10)
    # sub3 = Subscriber{Marker} # PlanTwoPointActionGoal  ## we want this to have initial velocity and attitude but it doesnt have it
    spin()
end

if ! isinteractive()
    main()
end