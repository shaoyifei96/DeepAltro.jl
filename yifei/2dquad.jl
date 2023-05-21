using Altro
using TrajectoryOptimization
using RobotDynamics
using StaticArrays, LinearAlgebra
const TO = TrajectoryOptimization
const RD = RobotDynamics

using ForwardDiff  # needed for @autodiff
using FiniteDiff   # needed for @autodiff
using DelimitedFiles
using Random
import YAML
using Distributions
config = YAML.load_file("/home/yifei/Documents/optimal_ctrl/ESE5460-Deep-iLQR/configs/Quad2D_datagen.yml")
# println(config)
# println(config["controller"]["horizon"])
# TODO: make all parameters to be in terms of file


RD.@autodiff struct quad2d <: RD.ContinuousDynamics end
RD.state_dim(::quad2d) = 6
RD.control_dim(::quad2d) = 2
g = 9.81 #m/s^2
m = 1.0 #kg 
ℓ = 0.3 #meters
J = 0.2*m*ℓ*ℓ
function RD.dynamics(::quad2d, x,u)
    θ = x[3]
    
    ẍ = (1/m)*(u[1] + u[2])*sin(θ)
    ÿ = (1/m)*(u[1] + u[2])*cos(θ) - g
    θ̈ = (1/J)*(ℓ/2)*(u[2] - u[1])
    # println(x)
    x =@SVector [x[4]; x[5]; x[6] ; ẍ; ÿ; θ̈]
    # println(x)
end

# Specify the default method to be used when calling the dynamics
#   options are `RD.StaticReturn()` or `RD.InPlace()`
RD.default_signature(::quad2d) = RD.StaticReturn()

# Specify the default method for evaluating the dynamics Jacobian
#   options are `RD.ForwardAD()`, `RD.FiniteDifference()`, or `RD.UserDefined()`
RD.default_diffmethod(::quad2d) = RD.ForwardAD()

model = quad2d()
dmodel = RD.DiscretizedDynamics{RD.RK4}(model)
n,m = size(model)    # get state and control dimension
N = 101              # number of time steps (knot points). Should be odd.
tf = 3.0             # total time (sec)
dt = tf / (N-1)      # time step (sec)

umin = [0.2*m*g; 0.2*m*g]
umax = [0.6*m*g; 0.6*m*g]

for data_idx = 1:config["datagen"]["num_samples"]
    # init_state = np.array([np.om.uniform(-3,3), np.random.uniform(-1,2), np.random.uniform(-math.pi/2.,math.pi/2.), np.random.uniform(-2,2), np.random.uniform(-2,2), np.random.uniform(-math.pi/4.,math.pi/4.)])
    x0 = SA_F64[rand(Uniform(-3,3)),rand(Uniform(-1,2)),rand(Uniform(-π/2, π/2)),rand(Uniform(-2,2)),rand(Uniform(-2,2)),rand(Uniform(-π/4, π/4))]   # start at the origin
    println(x0)
    xf = SA_F64[0,0,0,0,0,0]  # goal state

    Q  = Diagonal(SA[1.0,1.0,1.0,1.0,1.0,1.0])
    R  = Diagonal(SA[1e-3, 1e-3])
    Qf = 100Q
    obj = LQRObjective(Q,R,Qf,xf,N)

    cons = ConstraintList(n,m,N)

    # Goal constraint (first no constraints)
    # goal = GoalConstraint(xf)
    # add_constraint!(cons, goal, N)



    # bnd = BoundConstraint(n,m, u_min=umin, u_max=umax)
    # add_constraint!(cons, bnd, 1:N-1)  # add to all but the last time step
    # obsx = 0.7
    # obsy = 2
    # obsr = 0.15
    # obs = CircleConstraint(n, SA_F64[obsx], SA_F64[obsy], SA[obsr],1,2)#dim1 x, dim2 z
    # add_constraint!(cons, obs, 1:N-1)

    prob = Problem(model, obj, x0, tf, xf=xf, constraints=cons)
    initial_controls!(prob, [@SVector [0.5*m*g, 0.5*m*g] for k = 1:N-1])
    rollout!(prob)  

    opts = SolverOptions()
    opts.cost_tolerance = 1e-5
    opts.projected_newton = false
    # opts.line_search_lower_bound = 1e-8
    # opts.line_search_upper_bound = 0.1
    opts.verbose = 0
    println(opts.line_search_decrease_factor)
    # Create a solver, adding in additional options
    solver = ALTROSolver(prob, opts, show_summary=true; problem_idx = data_idx)
    # set_options!(solver, verbose = 5)
    solve!(solver)
    status(solver)

    X = states(solver)     # alternatively states(prob)
    U = controls(solver)   # alternatively controls(prob)
    ilqr = Altro.get_ilqr(solver)
    K = ilqr.K  # feedback gain matrices
    d = ilqr.d  # feedforward gains. Should be small.

end
# exit(0)

using MeshCat
using RobotZoo: Quadrotor, PlanarQuadrotor
using CoordinateTransformations, Rotations, Colors, StaticArrays, RobotDynamics

function set_mesh!(vis, model::L;
        scaling=1.0, color=colorant"black"
    ) where {L <: Union{Quadrotor, PlanarQuadrotor}} 
    # urdf_folder = joinpath(@__DIR__, "..", "data", "meshes")
    urdf_folder = @__DIR__
    # if scaling != 1.0
    #     quad_scaling = 0.085 * scaling
    obj = joinpath(urdf_folder, "quadrotor_scaled.obj")
    if scaling != 1.0
        error("Scaling not implemented after switching to MeshCat 0.12")
    end
    robot_obj = MeshFileGeometry(obj)
    mat = MeshPhongMaterial(color=color)
    setobject!(vis["robot"]["geom"], robot_obj, mat)
    if hasfield(L, :ned)
        model.ned && settransform!(vis["robot"]["geom"], LinearMap(RotX(pi)))
    end
end

function visualize!(vis, model::PlanarQuadrotor, x::StaticVector)
    py,pz = x[1], x[2]
    θ = x[3]
    settransform!(vis["robot"], compose(Translation(0,py,pz), LinearMap(RotX(-θ))))
end

function visualize!(vis, model, tf::Real, X)
    fps = Int(round((length(X)-1)/tf))
    anim = MeshCat.Animation(fps)
    n = RobotDynamics.state_dim(model)
    for (k,x) in enumerate(X)
        atframe(anim, k) do
            x = X[k]
            visualize!(vis, model, SVector{n}(x)) 
        end
    end
    setanimation!(vis, anim)
end

using GeometryBasics
using Plots
p = plot(X[:,1])
display(p)
# return


# vis = Visualizer()
# model = PlanarQuadrotor()
# set_mesh!(vis, model)
# render(vis)

# setobject!(vis, Sphere(Point(0,obsx,obsy), obsr))

# # setobject!(vis, Sphere(Point(xf[1],0,xf[2]), 0.1))
# # settransform!(vis, Translation(xf[1:3]))
# visualize!(vis, model, tf, X)



# print(X)
# X1 = [SVector{6}(x) for x in eachcol(xhist1)];
# X2 = [SVector{6}(x) for x in eachcol(xhist2)];