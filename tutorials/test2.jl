# Keyword: modular
# no. variables = no. constraints
# For testing -> equation solver not optimizer (<- later)
model = Model(Ipopt.Optimizer)
model = initialise_model(par)
model = reactor(model, par)
model = sep(model, par)
for i = eachindex(par.n_comp)
@NLconstraint(model, sep_input[i] - reactor_output[i] == 0)
end

Base.@kwdef mutable struct pb_SIP
    l_bdn::Vector
    u_bdn::Vector
    max_iter::Int = 50
    print::Int = 1 # 0 = no, 1 = iterates, 2 = trace also
    y_disc::Array{Float64}
    x_opt::Vector{Float64}
    greedy::Bool = true
    useBaron::Bool = false
end

pb = pb_SIP(y_low_ns, y_up_ns, iterations_limit, verbosity, y_disc, x_opt, greedy, useBaron) # variables defined earlier

Base.@kwdef mutable struct outer_struct
    SIPpart::pb_SIP
end

os = outer_struct(pb)

"""
    memoize(foo::Function, n_outputs::Int)

Take a function `foo` and return a vector of length `n_outputs`, where each
element is a function that returns the `i`'th output of `foo`.

To avoid duplication of work, cache the most-recent evaluations of `foo`.
Because `foo_i` is auto-differentiated with ForwardDiff, our cache needs to
work when `x` is a `Float64` and a `ForwardDiff.Dual`.
"""function memoize(foo::Function, n_outputs::Int)
    last_x, last_f = nothing, nothing
    last_dx, last_dfdx = nothing, nothingfunction foo_i(i, x::T...) where {T<:Real}
        if T == Float64if x != last_x
                last_x, last_f = x, foo(x...)
            endreturn last_f[i]::T
        elseif x != last_dx
                last_dx, last_dfdx = x, foo(x...)
            endreturn last_dfdx[i]::T
        endendreturn [(x...) -> foo_i(i, x...) for i in 1:n_outputs]
end

memoized_foo = memoize(foo, 2)
enthalpy = @NLexpression(model, memoize(heavy_enthalpy,10))

function foo(x, y)
    global function_calls += 1
    common_term = x^2 + y^2
    term_1 = sqrt(1 + common_term)
    term_2 = common_term
    return term_1, term_2
end

model = Model(Ipopt.Optimizer)
set_silent(model)
@variable(model, x[1:2] >= 0, start = 0.1)
register(model, :foo_1, 2, memoized_foo[1]; autodiff = true)
register(model, :foo_2, 2, memoized_foo[2]; autodiff = true)
@NLobjective(model, Max, foo_1(x[1], x[2]))
@NLconstraint(model, foo_2(x[1], x[2]) <= 2)
function_calls = 0
optimize!(model)
[memoized_foo[1](973), memoized_foo[2](973)]