df = DataFrame(Problem=String[],Algorithm=String[],x=Vector[],lbd=Float64[],termination=Integer[],lbds=[],global_solves=Integer[],y_disc_size=Float64[], Differentiable=Vector[],y_disc=Matrix[])
 
push!(df,(case,"BF",sol_BF.x_list[end],sol_BF.lbd_list[end],sol_BF.status, sol_BF.lbd_list, sol_BF.global_solves, length(sol_BF.y_disc[:,1]), [], sol_BF.y_disc))