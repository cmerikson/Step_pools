using DifferentialEquations, DomainSets, ModelingToolkit, MethodOfLines, Printf, Plots, LaTeXStrings, DataFrames, CSV

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

# Function to creat animation from simulation results
function create_animations(df::DataFrame, peaks::Tuple{Float64, Float64})
    FilteredData = filter(row -> row.Peak1 == peaks[1] && row.Peak2 == peaks[2], df)
    anim = @animate for time in unique(FilteredData.t)
        Data = filter(row -> row.t == time, FilteredData)
        #format_time = scientific_notation(Data.t)
        
        p1 = plot(Data.x, Data.η, ylim=(-1,1), color="#9d162e", title=latexstring("\$t^*\$ = $(time)"), xlabel=latexstring("Dimensionless Distance, \$x^*\$"), ylabel=latexstring("Dimensionless Elevation, \$η^{*}(x^*,t^*)\$"), label="Elevation")
        hline!([0.0], linestyle=:dash, color="gray", label="")
        p2 = plot(Data.x, Data.u, ylim=(0.5,0.8), color="#003c73", title=latexstring("\$t^*\$ = $(time)"), xlabel=latexstring("Dimensionless Distance, \$x^*\$"), ylabel=latexstring("Dimensionless Velocity, \$u^{*}(x^*,t^*)\$"), label="Velocity")
        plot!(twinx(), Data.x, Data.roughness, ylim=(2,7), yticks=[], color="black", label="Roughness")
        
        plot(p1, p2, layout=(2, 1), size=(800, 600))
    end
    return anim
end

# Constants
τ_cr = 0.05 # Critical Shields Number
g = 9.81 # Gravitational Acceleration m/s^2
κ = 0.41 # Von Karman constant
r_0 = 0.5 # Immobile Clast Radius
D = 0.25 # Mobile Layer Grain Size
D_min = 0.25
λ = 30 # Roughness patch spacing
k = (2*pi)/λ # Wavenumber
ρ_s = 2650 # Sediment density
ρ = 1000 # Fluid density
θ = 0.04996 # Angle of slope in radians (7% 0.06989)
γ = 1

function compute_dimensionless(q_w)
    ζ = (q_w^(2.0/3.0))/((g^(1.0/3.0))*r_0)
    α = (ρ*(g*q_w)^(2.0/3.0) * κ^2.0) / ((ρ_s - ρ)*g*D)
    β = (λ*κ^2.0)/(r_0)
    ψ = λ / r_0
    return ζ, α, β, ψ
end

# Numerical Solver
function Gaussian_DualPeak_Model(Peak1::Float64, Peak2::Float64, Discharge::Float64, timestep::Float64, simtime::Float64)
    q_w = Discharge # Discharge per width

    ζ, α, β, ψ = compute_dimensionless(q_w)
    
    @parameters x, t
    @variables u(..) [bounds=(0,Inf)]
    @variables η(..) 
    
    # Linear Operators
    Dx = Differential(x)
    Dt = Differential(t)
    
    # Gaussian Peaked Roughness
    peak1 = (1.0 - D/r_0)*exp((-(x - Peak1)^2.0)/0.1) + (1.0 + D/r_0)
    peak2 = (1.0 - D/r_0)*exp((-(x - Peak2)^2.0)/0.1) + (1.0 + D/r_0)

    r = peak1 + peak2

    log_term = (log(ζ / (exp(1)*u(x, t) * r)))

    eq1 = ζ * u(x,t) * Dx(u(x,t)) ~ -ζ * expand_derivatives(Dx(1/u(x,t))) - Dx(η(x,t)) - β*u(x,t)^3.0 * (1.0 - ((2/ζ) * u(x,t) * r) + ((1/ζ^(2.0)) * u(x,t)^(2.0) * r^(2.0))) * ((log_term + ((1/ζ) * u(x,t) * r))^(-2.0)) + ψ*tan(θ)
    
    log_dx = Dx((α * u(x,t)^2.0) * (1 - ((2/ζ) * u(x,t) * r) + ((1/ζ^(2.0)) * u(x,t)^(2.0) * r^(2.0))) * (log_term + ((1/ζ) * u(x,t) * r))^(-2.0) - τ_cr)
    eq2B = -expand_derivatives(log_dx^γ) * (1/α)
    τ = ((α*u(x,t)^2.0) * (1 - ((2/ζ) * u(x,t) * r) + ((1/ζ^(2.0)) * u(x,t)^(2.0) * r^(2.0))) * (log_term + (u(x,t)*r / ζ))^(-2.0))
    f(x,t) = ((τ >= τ_cr) * eq2B)
    eq2 = Dt(η(x,t)) ~ f(x,t)
    
    eq = [eq1, eq2]

    # Domain
    x_start = t_start = 0.0
    x_end = 3.0
    t_end = simtime
    
    domains = [x ∈ IntervalDomain(x_start,x_end), t ∈ IntervalDomain(t_start,t_end)]
    
    # Periodic Boundary Conditions
    η0(x, t) = 0.0
    u0(x, t) = 2.0
    
    bcs = [η(x, 0.0) ~ η0(x, 0.0),
            η(x_start, t) ~ η(x_end, t),
            u(x_start,t) ~ u(x_end,t),
            u(x,0.0) ~ u0(x,0.0)]
    
    @named pdesys = PDESystem(eq, bcs, domains, [x, t], [u(x, t), η(x, t)])

    discretization = MOLFiniteDifference([x => 100], t)

    prob = discretize(pdesys, discretization);

    sol1 = solve(prob, Rosenbrock23(), saveat = timestep, maxiters=1e5)

    disct = sol1[t]
    discx = sol1[x]
    discu = sol1[u(x,t)]
    discη = sol1[η(x,t)]

    # Ensure the matrices have the same dimensions
    rows, cols = size(discu)
    @assert size(discη) == (rows, cols)
    
    # Create vectors for each matrix
    discu_vec = vec(discu')
    discη_vec = vec(discη')
    
    # Create additional columns for time slices and spatial indices
    time_slices = repeat(range(0.0, stop=50.0, length=cols), outer=rows)
    spatial_indices = repeat(range(0.0, stop=3.0, length=rows), inner=cols)
    
    discr = (1.0 - D/r_0).*exp.((-(discx .- Peak1).^2.0)./0.1) .+ (1.0 + D/r_0) .+ ((1.0 - D/r_0).*exp.((-(discx .- Peak2).^2.0)./0.1) .+ (1.0 + D/r_0) .+ (1.0 - D/r_0))
    
    roughness = repeat(discr, inner=cols)
    
    df = DataFrame(
        t = time_slices,
        x = spatial_indices,
        u = discu_vec,
        η = discη_vec,
        roughness = roughness,
    )
    return df
end

function run_simulations(peaks::Vector{Tuple{Float64, Float64}}, Discharge::Float64, timestep::Float64, simtime::Float64)
    results = []
    for (Peak1, Peak2) in peaks
        df = Gaussian_DualPeak_Model(Peak1, Peak2, Discharge, timestep, simtime)
        df[:, :Peak1] .= Peak1
        df[:, :Peak2] .= Peak2
        push!(results, df)
    end
    return vcat(results...)
end

peaks = [(0.0,3.0), (0.5, 2.5), (1.0, 2.0), (1.25, 1.75)]
timestep = 0.1
simtime = 50.0
Discharge = 50.0

simulation_results = run_simulations(peaks, Discharge, timestep, simtime)

CSV.write("Gaussian_Results.csv", simulation_results)

for i in peaks
    anim = create_animations(simulation_results, i)
    gif(anim, "simulation_results_$(i)peak.gif", fps=30)
end

