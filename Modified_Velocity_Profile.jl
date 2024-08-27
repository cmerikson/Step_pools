using DifferentialEquations, DomainSets, ModelingToolkit, MethodOfLines, Printf, Plots, LaTeXStrings, DataFrames, CSV, Symbolics, NonlinearSolve, StatsBase, NaNMath

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

# Function to create animation from simulation results
function create_animations(df::DataFrame, peaks::Tuple{Float64, Float64, Union{Nothing, Float64}})
    println("Creating animation for peaks: $peaks")
    FilteredData = filter(row -> row.Peak1 == peaks[1] && row.Peak2 == peaks[2] &&
                          (peaks[3] === nothing ? row.Peak3 === NaN : row.Peak3 == peaks[3]), df)
    println("Filtered data has $(nrow(FilteredData)) rows")

    if nrow(FilteredData) == 0
        println("No data available for peaks: $peaks")
        return nothing
    end

    anim = @animate for time in unique(FilteredData.t)
        Data = filter(row -> row.t == time, FilteredData)
        
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
ρ_s = 2650 # Sediment density
ρ = 1000 # Fluid density
γ = 1

# Reference angle of slope in radians (10% 0.09967, 8% 0.07983, 7% 0.06989 5% 0.04996)

function compute_dimensionless(q_w,λ,r_0,D)
    ζ = (q_w^(2.0/3.0))/((g^(1.0/3.0))*r_0)
    α = (ρ*(g*q_w)^(2.0/3.0) * κ^2.0) / ((ρ_s - ρ)*g*D)
    β = (λ*κ^2.0)/(r_0)
    ψ = λ / r_0
    return ζ, α, β, ψ
end

# Numerical Solver
function Gaussian_DualPeak_Model(Peak1::Float64, Peak2::Float64, Discharge::Float64, Roughness_Scale::Float64, Length_Scale::Float64, D50::Float64, Slope::Float64, timestep::Float64, simtime::Float64; Peak3::Union{Nothing, Float64}=nothing, u_init::Float64=1.0, Bankfull::Bool=false, Reference::String=nothing)
    q_w = Discharge # Discharge per width
    θ = Slope
    r_0 = Roughness_Scale
    D = D50
    λ = Length_Scale

    ζ, α, β, ψ = compute_dimensionless(q_w,λ,r_0,D)
    
    @parameters x, t
    @variables u(..) [bounds=(0,Inf)]
    @variables η(..) 
    
    # Linear Operators
    Dx = Differential(x)
    Dt = Differential(t)
    
    # Gaussian Peaked Roughness
    peak1 = (1.0 - D/r_0)*exp((-(x - Peak1)^2.0)/0.05) + (1.0 + D/r_0)
    peak2 = (1.0 - D/r_0)*exp((-(x - Peak2)^2.0)/0.05) + (1.0 + D/r_0)
    
    if Peak3 !== nothing
        peak3 = (1.0 - D/r_0)*exp((-(x - Peak3)^2.0)/0.05) + (1.0 + D/r_0)
        r = gau = peak1 + peak2 + peak3

        #r = ifelse(η(x,t) > gau, η(x,t), gau)
    else
        r = gau = peak1 + peak2
        
        #r= ifelse(η(x,t) > gau, η(x,t), gau)
    end

    log_termA = (NaNMath.log((1.0/2.0) + ((1.0*ζ) / (2.0*u(x, t) * r))))
    log_termB = (NaNMath.log(1.0 + ((1.0*ζ) / (u(x, t) * r))))

    eq1 = ζ * u(x,t) * Dx(u(x,t)) ~ -ζ * expand_derivatives(Dx(1/u(x,t))) - Dx(η(x,t)) - β*u(x,t)^3.0 * (1 - ((2/ζ) * u(x,t) * r) + ((1/ζ^(2.0)) * u(x,t)^(2.0) * r^(2.0))) * ((((1/(1.0*ζ))*log_termA*(u(x, t) * r)) + log_termB + ((1.0/ζ)*(u(x, t) * r) * log(2.0*exp(1))) - 1.0)^(-2.0)) + ψ*tan(θ)
    
    log_dx = Dx((α * u(x,t)^2.0) * (1 - ((2/ζ) * u(x,t) * r) + ((1/ζ^(2.0)) * u(x,t)^(2.0) * r^(2.0))) * (((1/(1.0*ζ))*log_termA*(u(x, t) * r)) + log_termB + ((1.0/ζ)*(u(x, t) * r) * log(2.0*exp(1))) - 1.0)^(-2.0) - τ_cr)
    eq2B = -expand_derivatives(log_dx^γ) * (1/α)
    τ = (α*u(x,t)^2.0) * ((1 - ((2/ζ) * u(x,t) * r) + ((1/ζ^(2.0)) * u(x,t)^(2.0) * r^(2.0))) * (((1/(1.0*ζ))*log_termA*(u(x, t) * r)) + log_termB + ((1.0/ζ)*(u(x, t) * r) * log(2.0*exp(1))) - 1.0)^(-2.0))
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
    u0(x, t) = u_init
    
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
    time_slices = repeat(range(0.0, stop=simtime, length=cols), outer=rows)
    spatial_indices = repeat(range(0.0, stop=3.0, length=rows), inner=cols)
    
    discr = (1.0 - D/r_0).*exp.((-(discx .- Peak1).^2.0)./0.1) .+ (1.0 + D/r_0) .+ ((1.0 - D/r_0).*exp.((-(discx .- Peak2).^2.0)./0.1) .+ (1.0 + D/r_0))
    if Peak3 !== nothing
        discr .+= ((1.0 - D/r_0).*exp.((-(discx .- Peak3).^2.0)./0.1) .+ (1.0 + D/r_0))
    end
    
    roughness = repeat(discr, inner=cols)

    τ_values = [Symbolics.value(substitute(τ, Dict(x => i, t => j, u(x, t) => k, η(x,t) => m))) for (i, j, k, m) in zip(spatial_indices, time_slices, discu_vec, discη_vec)]

    function calculate_bankfull(τ_cr::Float64)

        f(u,r) = Symbolics.value(α) .* u .* u .* (1 .- (2 .* u .* Symbolics.value(r) ./ Symbolics.value(ζ)) .+ (u .* u .* Symbolics.value(r).^2 ./ Symbolics.value(ζ).^2)) .* ((NaNMath.log.(Symbolics.value(ζ) ./ (exp(1) .* u .* Symbolics.value(r)))) .+ (u .* Symbolics.value(r) ./ Symbolics.value(ζ))).^(-2.0) .- τ_cr
        
        u0 = repeat([1.0],length(discr))

        problem = NonlinearProblem(f, u0, discr)
        sol = solve(problem)
        mean_u = geomean(sol.u)
        return mean_u
    end

    if Bankfull
        bankfull_u = calculate_bankfull(τ_cr)
        df = DataFrame(
        t = time_slices,
        x = spatial_indices,
        u = discu_vec,
        η = discη_vec,
        roughness = roughness,
        tau = τ_values,
        q = q_w,
        slope = θ,
        r0 = r_0,
        bankfull = repeat([bankfull_u], length(time_slices)),
        Reference = Reference
        )
    else
        df = DataFrame(
        t = time_slices,
        x = spatial_indices,
        u = discu_vec,
        η = discη_vec,
        roughness = roughness,
        tau = τ_values,
        q = q_w,
        slope = θ,
        r0 = r_0,
        Reference = Reference
        )
    end
    
    df[:, :Peak1].=Peak1
    df[:, :Peak2].=Peak2
    df[:, :Peak3] .= Peak3 === nothing ? NaN : Peak3
 
    return df
end

function run_simulations(peaks::Vector{Tuple{Float64, Float64, Union{Nothing, Float64}}}, Discharge::Float64, Roughness_Scale::Float64, Length_Scale::Float64, D50::Float64, Slope::Float64, timestep::Float64, simtime::Float64; u_init::Float64=1.0, Bankfull::Bool=false, Reference::String=nothing)
    results = [Gaussian_DualPeak_Model(Peak1, Peak2, Discharge, Roughness_Scale, Length_Scale, D50, Slope, timestep, simtime, Peak3=Peak3, Bankfull=Bankfull, Reference=Reference) for (Peak1, Peak2, Peak3) in peaks]
    return vcat(results...)
end

function create_animations(simulation_results, output_directory)
    for i in peaks
        anim = create_animations(simulation_results, i)
        gif(anim, string(output_directory,"/simulation_results_$(i[1])_$(i[2])_$(i[3] === nothing ? "none" : string(i[3]))_peak.gif"), fps=30)
    end
end
