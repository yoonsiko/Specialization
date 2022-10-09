msg = "Hello World"
println(msg)

println("Macro, a special type of function")
@show 1 + 1
@show 1 - 2
@show 2 * 2
@show 4 / 5
@show 3^2

@show typeof(1)
@show typeof(1.0)

println("\nComplex numbers")
x = 2 + 1im
@show real(x)
@show imag(x)
@show typeof(x)
@show x * (1 - 2im)

println("\nBackslash pi then tab to write symbol of pi")
@show typeof(π)

println("\nBoolean")
@show 0.1 * 3 == 0.3
@show sin(2π / 3) == √3 / 2

println("\nBoolean with approx")
@show 0.1 * 3 ≈ 0.3
@show sin(2π / 3) ≈ √3 / 2
@show isapprox(sin(2π / 3), √3 / 2; atol = 1e-8)

@show 1e-300 ≈ 0.0
@show isapprox(1e-9, 0.0; atol = 1e-8)

println("\nIf you aren't careful, floating point arithmetic can throw up all manner of issues")
1 + 1e-16 == 1

println("\nVector and Matrix")

b = [5, 6]
A = [1.0 2.0; 3.0 4.0]
@show typeof(b)
@show typeof(A)
@show x = A \ b
@show A * x
@show A * x ≈ b
# @show b * b  This wont work
@show b' * b
@show b * b'

println("\nOther common types")
@show typeof("This is Julia")
@show typeof("π is about 3.1415")
x = 123
println("The value of x is: $(x)")


println("\nSymbols")
println("The value of x is: $(eval(:x))")
typeof(:x)
String(:abc)
Symbol("abc")


println("\nTuples")
t = ("hello", 1.2, :foo)
@show typeof(t)
@show t[2]
a, b, c = t
@show b
t = (word = "hello", num = 1.2, sym = :foo)
@show t.word

println("\nDictionaries")
d1 = Dict(1 => "A", 2 => "B", 4 => "D")
d1[2]
Dict("A" => 1, "B" => 2.5, "D" => 2 - 3im)
Dict("A" => 1.0 + 0.0im, "B" => 2.5 + 0.0im, "D" => 2.0 - 3.0im)
d2 = Dict("A" => 1, "B" => 2, "D" => Dict(:foo => 3, :bar => 4))
d2["B"]
d2["D"][:foo]

println("\nStructs")
struct MyStruct
    x = 12
    y = "a"
    z = 1
end

a = MyStruct(1, "a", Dict(2 => 3))
a.x

# Structs are immutable by default, but we can construct a mutable struct:
mutable struct MyStructMutable
    x::Int = 12
    y::String = "a"
    z::Dict{Int,Int} = Dict(2 => 3)
end

a = MyStructMutable(1, "a", Dict(2 => 3))
a.x

println("\nLoops")
# Julia has native support for for-each style loops with the syntax for <value> in <collection> end:
for i in 1:5
    println(i)
end

# Ranges are constructed as start:stop and start:step:stop
for i in 1.2:1.1:5.6
    println(i)
end

# This for-each loop also works with dictionaries:
for (key, value) in Dict("A" => 1, "B" => 2.5, "D" => 2 - 3im)
    println("$(key): $(value)")
end


println("\nControl flow")
# Julia control flow is similar to Matlab, 
# using the keywords if-elseif-else-end, and the logical operators || and && for or and and respectively:

for i in 0:5:15
    if i < 5
        println("$(i) is less than 5")
    elseif i < 10
        println("$(i) is less than 10")
    else
        if i == 10
            println("the value is 10")
        else
            println("$(i) is bigger than 10")
        end
    end
end


println("\nComprehensions")
# Similar to languages like Haskell and Python, 
# Julia supports the use of simple loops in the construction of arrays and dictionaries, called comprehensions.

# A list of increasing integers:
[i for i in 1:5]

# Matrices can be built by including multiple indices:
[i * j for i in 1:5, j in 5:10]

# Conditional statements can be used to filter out some values:
[i for i in 1:10 if i % 2 == 1]

# A similar syntax can be used for building dictionaries:
Dict("$(i)" => i for i in 1:10 if i % 2 == 1)


println("\nFunctions")

# A simple function is defined as follows:
function print_hello()
    return println("hello")
end
print_hello()

# Arguments can be added to a function:
function print_it(x)
    return println(x)
end
print_it("hello")
print_it(1.234)
print_it(:my_id)

# Optional keyword arguments are also possible:
function print_it(x; prefix = "value:")
    return println("$(prefix) $(x)")
end
print_it(1.234)
print_it(1.234; prefix = "val:")

# The keyword return is used to specify the return values of a function:
function mult(x; y = 2.0)
    return x * y
end

mult(4.0)
mult(4.0; y = 5.0)

println("\nAnonymous functions")
# The syntax input -> output creates an anonymous function. These are most useful when passed to other functions
# For example: 
f = x -> x^2
f(2)

map(x -> x^2, 1:4)

println("\nType parameters")
# We can constrain the inputs to a function using type parameters,
# which are :: followed by the type of the input we want. 

# For example:
function foo(x::Int)
    return x^2
end

function foo(x::Float64)
    return exp(x)
end

function foo(x::Number)
    return x + 1
end

@show foo(2)
@show foo(2.0)
@show foo(1 + 1im)
 
# But what happens if we call foo with something we haven't defined it for?
# foo([1, 2, 3]) This gives an error

println("\nBroadcasting")
# In the example above, we didn't define what to do if f was passed a Vector.
# Luckily, Julia provides a convenient syntax for mapping f element-wise over arrays!
# Just add a . between the name of the function and the opening
# (. This works for any function, including functions with multiple arguments. 
#For example:
f.([1, 2, 3])

println("\nMutable vs immutable objects")
# Some types in Julia are mutable, which means you can change the values inside them.
# A good example is an array. You can modify the contents of an array without having to make a new array.

# In contrast, types like Float64 are immutable. You can't modify the contents of a Float64.
# This is something to be aware of when passing types into functions. For example:

function mutability_example(mutable_type::Vector{Int}, immutable_type::Int)
    mutable_type[1] += 1
    immutable_type += 1
    return
end

mutable_type = [1, 2, 3]
immutable_type = 1

mutability_example(mutable_type, immutable_type)

println("mutable_type: $(mutable_type)")
println("immutable_type: $(immutable_type)")

# Because Vector{Int} is a mutable type, modifying the variable inside the 
# function changed the value outside of the function. In contrast,
# the change to immutable_type didn't modify the value outside the function.

#You can check mutability with the isimmutable function:
isimmutable([1, 2, 3]) # false
isimmutable(1) #true


println("\nPackages")
using Random  # The equivalent of Python's `from Random import *`
import Random  # The equivalent of Python's `import Random`

Random.seed!(33)

[rand() for i in 1:10]

using Pkg
#Pkg.add("JuMP")

#Pkg.add("https://github.com/user-name/MyPackage.jl.git")

