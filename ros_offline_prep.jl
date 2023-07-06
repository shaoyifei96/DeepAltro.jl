import Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()
using RobotOS
@rosimport kr_planning_msgs.msg: TrajectoryDiscretized, SplineTrajectory
@rosimport visualization_msgs.msg: Marker
@rosimport decomp_ros_msgs.msg: Polyhedron, PolyhedronArray
@rosimport geometry_msgs.msg: Point
@rosimport std_msgs.msg: ColorRGBA, Float64, Int32
rostypegen()