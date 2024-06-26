using Printf, ModelingToolkit, DomainSets, DifferentialEquations, MethodOfLines, Plots, LaTeXStrings

# Function to convert a number to scientific notation with 10^power format
function scientific_notation(value)
    if value == 0
        return "0"
    else
        exponent = floor(log10(abs(value)))
        coefficient = value / (10.0^exponent)
        return @sprintf("\$ %.2f \\times 10^{%d} \$", coefficient, Int(exponent))
    end
end

# Constants
τ_cr = 0.05 # Critical Shields Number
g = 9.81 # Gravitational Acceleration m/s^2
κ = 0.41 # Von Karman constant
q_w = 30 # Discharge per width
r_0 = 0.5 # Immobile Clast Radius
D = 0.25 # Mobile Layer Grain Size
D_min = 0.25
λ = 30 # Roughness patch spacing
k = (2*pi)/λ # Wavenumber
ρ_s = 2650 # Sediment density
ρ = 1000 # Fluid density
θ = 0.04996 # Angle of slope in radians (7% 0.06989)
γ = 1

# Dimensionless Parameters
ζ = (q_w^(2.0/3.0))/((g^(1.0/3.0))*r_0)
α = (ρ*(g*q_w)^(2.0/3.0) * κ^2.0) / ((ρ_s - ρ)*g*D)
β = (λ*κ^2.0)/(r_0)
ψ = λ / r_0

@parameters x, t
@variables u(..) [bounds=(0,Inf)]
@variables η(..) 

# Linear Operators
Dx = Differential(x)
Dt = Differential(t)

# Periodic Roughness
r = (1-(D_min/r_0))*cos(2.0*pi*x) + (1+(D_min/r_0))

# Term inside the Logarithm
log_term = (log(ζ / (exp(1)*u(x, t) * r)))

# Fluid Momentum
eq1 = ζ * u(x,t) * Dx(u(x,t)) ~ -ζ * expand_derivatives(Dx(1/u(x,t))) - Dx(η(x,t)) - β*u(x,t)^3.0 * (1.0 - ((2/ζ) * u(x,t) * r) + ((1/ζ^(2.0)) * u(x,t)^(2.0) * r^(2.0))) * ((log_term + ((1/ζ) * u(x,t) * r))^(-2.0)) + ψ*tan(θ)

# Mass Conservation
log_dx = Dx((α * u(x,t)^2.0) * (1 - ((2/ζ) * u(x,t) * r) + ((1/ζ^(2.0)) * u(x,t)^(2.0) * r^(2.0))) * (log_term + ((1/ζ) * u(x,t) * r))^(-2.0) - τ_cr)
eq2B = -expand_derivatives(log_dx^γ) * (1/α)
τ = ((α*u(x,t)^2.0) * (1 - ((2/ζ) * u(x,t) * r) + ((1/ζ^(2.0)) * u(x,t)^(2.0) * r^(2.0))) * (log_term + (u(x,t)*r / ζ))^(-2.0))
f(x,t) = ((τ >= τ_cr) * eq2B)
eq2 = Dt(η(x,t)) ~ f(x,t)

# Combine into a system
eq = [eq1, eq2]

# Domain
x_start = t_start = 0.0
x_end = 3
t_end = 5.0

domains = [x ∈ IntervalDomain(x_start,x_end), t ∈ IntervalDomain(t_start,t_end)]

# Periodic Boundary Conditions
η0(x, t) = 0.0
u0(x, t) = 2.0

bcs = [η(x, 0.0) ~ η0(x, 0.0),
        η(x_start, t) ~ η(x_end, t),
        u(x_start,t) ~ u(x_end,t),
        u(x,0.0) ~ u0(x,0.0)]

@named pdesys = PDESystem(eq, bcs, domains, [x, t], [u(x, t), η(x, t)])

# Finite Difference Method
discretization = MOLFiniteDifference([x => 100], t)

# Define the differential equations problem
prob = discretize(pdesys, discretization);

# Solve PDE system
sol1 = solve(prob, Rosenbrock23(), saveat = 0.1, maxiters=1e5)

disct = sol1[t]
discx = sol1[x]
discu = sol1[u(x,t)]
discη = sol1[η(x,t)]

anim = @animate for i in 1:length(disct)
    η = discη[:, i]
    format_time = scientific_notation(disct[i])
    p1 = plot(discx, η, ylim=(-2,4), xlim=(0,3), title=latexstring("\$t^*\$ = $(format_time)"), 
              color="blue", xlabel=latexstring("Dimensionless Distance, \$x^*\$"), ylabel=latexstring("Dimensionless Elevation, \$η^{*}(x^*,t^*)\$"),
        label=latexstring("\$q_w\$ = $(q_w)"), legend=:topleft)
end
gif(anim, "Step_pools.gif", fps=50)
