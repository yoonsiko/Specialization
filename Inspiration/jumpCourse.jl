#import Pkg
#Pkg.add("JuMP")
#Pkg.add("HiGHS") 

# Once JuMP is installed, to use JuMP in your programs write:
using JuMP

# We also need to include a Julia package which provides an appropriate solver. 
# We want to use HiGHS.Optimizer here which is provided by the HiGHS.jl package.
using HiGHS
# using Ipopt 

# JuMP builds problems incrementally in a Model object. Create a model by passing an optimizer to the Model function:
model = Model(HiGHS.Optimizer)

# Variables are modeled using @variable:
@variable(model, x >= 0)

# They can have lower and upper bounds
@variable(model, 0 <= y <= 3)

# The objective is set using @objective:
@objective(model, Min, 12x + 20y)

# Constraints are modeled using @constraint. Here, c1 and c2 are the names of our constraint.
# c1 : 6x+8y≥100.0
@constraint(model, c1, 6x + 8y >= 100)

# c2 : 7x+12y≥120.0
@constraint(model, c2, 7x + 12y >= 120)

# Call print to display the model:
print(model)

# To solve the optimization problem, call the optimize! function.
optimize!(model)

# The ! after optimize is part of the name. It's nothing special. 
# Julia has a convention that functions which mutate their arguments should end in !. A common example is push!.

# Now let's see what information we can query about the solution.
# termination_status tells us why the solver stopped:
@show termination_status(model)

# In this case, the solver found an optimal solution.
# Check primal_status to see if the solver found a primal feasible point:
@show primal_status(model)

# and dual_status to see if the solver found a dual feasible point:
@show dual_status(model)

# Now we know that our solver found an optimal solution, and that it has a primal and a dual solution to query.
# Query the objective value using objective_value:
@show objective_value(model)

# the primal solution using value:
@show value(x)
@show value(y)

# and the dual solution using shadow_price:
@show shadow_price(c1)
@show shadow_price(c2)


println("\nModel Basics")
# Create a model by passing an optimizer:
model2 = Model(HiGHS.Optimizer)

# Alternatively, call set_optimizer at any point before calling optimize!:
model3 = Model()
set_optimizer(model3, HiGHS.Optimizer)

# For some solvers, you can also use direct_model, which offers a more efficient connection to the underlying solver:
model4 = direct_model(HiGHS.Optimizer())
# Warning, Some solvers do not support direct_model!

println("\nSolver Options")
# Pass options to solvers with optimizer_with_attributes:
model5 = Model(optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false))

# You can also pass options with set_optimizer_attribute

model6 = Model(HiGHS.Optimizer)
set_optimizer_attribute(model6, "output_flag", false)

println("\nSolution basics")
# We saw above how to use termination_status and primal_status to understand the solution returned by the solver.

# However, only query solution attributes like value and objective_value if there is an available solution.
# Here's a recommended way to check:

function solve_infeasible()
    model = Model(HiGHS.Optimizer)
    @variable(model, 0 <= x <= 1)
    @variable(model, 0 <= y <= 1)
    @constraint(model, x + y >= 3)
    @objective(model, Max, x + 2y)
    optimize!(model)
    if termination_status(model) != OPTIMAL
        @warn("The model was not solved correctly.")
        return nothing
    end
    return value(x), value(y)
end

solve_infeasible()


println("\nVariable bounds")
# Let's create a new empty model to explain some of the variable syntax:

model7 = Model()
@variable(model7, free_x)
@variable(model7, keyword_x, lower_bound = 1, upper_bound = 2)
has_upper_bound(keyword_x)
upper_bound(keyword_x)
# lower_bound(free_x) <- this wont work


println("\nContainers")
# We have already seen how to add a single variable to a model using the @variable macro. 
# Now let's look at ways to add multiple variables to a model.

# JuMP provides data structures for adding collections of variables to a model.
# These data structures are referred to as containers and are of three types: Arrays, DenseAxisArrays, and SparseAxisArrays.
println("\nArrays")
# JuMP arrays are created when you have integer indices that start at 1:
@variable(model7, a[1:2, 1:2])

# Create an n-dimensional variable x∈Rn with bounds l≤x≤u (l,u∈Rn) as follows:
n = 10
l = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
u = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

@variable(model, l[i] <= x[i = 1:n] <= u[i])


# We can also create variable bounds that depend upon the indices:
@variable(model, y[i = 1:2, j = 1:2] >= 2i + j)


println("\nDenseAxisArrays")
# DenseAxisArrays are used when the indices are not one-based integer ranges. 
# The syntax is similar except with an arbitrary vector as an index as opposed to a one-based range:
@variable(model, z[i = 2:3, j = 1:2:3] >= 0)

# Indices do not have to be integers. They can be any Julia type:
@variable(model, w[1:5, ["red", "blue"]] <= 1)

println("\nSparseAxisArrays")
# SparseAxisArrays are created when the indices do not form a Cartesian product.
# For example, this applies when indices have a dependence upon previous indices (called triangular indexing):
@variable(model, u[i = 1:2, j = i:3])

# We can also conditionally create variables by adding a comparison 
# check that depends upon the named indices and is separated from the indices by a semi-colon ;:
@variable(model, v[i = 1:9; mod(i, 3) == 0])

println("\nIntegrality")
# JuMP can create binary and integer variables.
# Binary variables are constrained to the set {0,1}, and integer variables are constrained to the set Z.

# Create an integer variable by passing Int:
@variable(model, integer_x, Int)

# or setting the integer keyword to true:
@variable(model, integer_z, integer = true)

# Create a binary variable by passing Bin:
@variable(model, binary_x, Bin)

# or setting the binary keyword to true:
@variable(model, binary_z, binary = true)

print("\nConstraint basics")
# We'll need a need a new model to explain some of the constraint basics:
model = Model()
@variable(model, x)
@variable(model, y)
@variable(model, z[1:10]);

print("\nContainers")
# Just as we had containers for variables, JuMP also provides Arrays, DenseAxisArrays,
# and SparseAxisArrays for storing collections of constraints. Examples for each container type are given below.

# Arrays
@constraint(model, [i = 1:3], i * x <= i + 1)

# DenseAxisArrays
@constraint(model, [i = 1:2, j = 2:3], i * x <= j + 1)

# SparseAxisArrays
@constraint(model, [i = 1:2, j = 1:2; i != j], i * x <= j + 1)

print("\nConstraints in a loop")
# We can add constraints using regular Julia loops:
for i in 1:3
    @constraint(model, 6x + 4y >= 5i)
end

# or use for each loops inside the @constraint macro:
@constraint(model, [i in 1:3], 6x + 4y >= 5i)

# We can also create constraints such as ∑10i=1zi≤1:
@constraint(model, sum(z[i] for i in 1:10) <= 1)


println("\nObjective functions")
# Set an objective function with @objective:
model = Model(HiGHS.Optimizer)
@variable(model, x >= 0)
@variable(model, y >= 0)
@objective(model, Min, 2x + y)

# Create a maximization objective using Max:
@objective(model, Max, 2x + y)

# Calling @objective multiple times will over-write the previous objective. 
# This can be useful when you want to solve the same problem with different objectives.

println("\nVectorized syntax")
# We can also add constraints and an objective to JuMP using vectorized linear algebra.
# We'll illustrate this by solving an LP in standard form i.e.
vector_model = Model(HiGHS.Optimizer)

A = [
    1 1 9 5
    3 5 0 8
    2 0 6 13
]

b = [7; 3; 5]

c = [1; 3; 5; 2]

@variable(vector_model, x[1:4] >= 0)
@constraint(vector_model, A * x .== b)
@objective(vector_model, Min, c' * x)

optimize!(vector_model)

objective_value(vector_model)
