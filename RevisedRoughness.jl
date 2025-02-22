using DifferentialEquations, DomainSets, ModelingToolkit, MethodOfLines, Printf, Plots, LaTeXStrings, DataFrames, CSV, Symbolics, NonlinearSolve, StatsBase, NaNMath, ColorSchemes, ProgressLogging

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

# Create Plot
function PlotEvolution(disct, discx, discη, roughness, simulation_time, time_step)
    colors = [get(ColorSchemes.oslo, i / 4) for i in 0:3]
    slices = simulation_time / time_step

    # Initial Bed surface
    plot(discx, discη[:,1], xlim=(0,20), ylim=(-1,1), title=latexstring("\$η^*\$ Evolution"), 
    color=colors[4], xlabel=latexstring("Dimensionless Distance, \$x^*\$"), ylabel=latexstring("Dimensionless Elevation, \$η^{*}(x^*,t^*)\$"),
    label=latexstring("\$t^{*} = 0\$"), linestyle=:dash, legend=:topleft)

    plot!(discx,discη[:,Int(round(slices*0.25))],label=latexstring("\$t^{*} = \$ $(scientific_notation(disct[Int(round(slices*0.25))]))"), color=colors[3], linewidth=4)
    plot!(discx,discη[:,Int(round(slices*0.5))],label=latexstring("\$t^{*} = \$ $(scientific_notation(disct[Int(round(slices*0.5))]))"), color=colors[2], linewidth=4)
    #plot!(discx,discη[:,Int(round(slices*0.75))],label=latexstring("\$t^{*} = \$ $(scientific_notation(disct[Int(round(slices*0.75))]))"), color="#003c73")
    plot!(discx,discη[:,end],label=latexstring("\$t^{*} = \$ $(scientific_notation(disct[end]))"), color=colors[1], linewidth=4)
    p = plot!(twinx(), discx, roughness, color="grey", label="Roughness", ylim=(1,400), yticks=false)
    return p
end

function PlotFrame(df::DataFrame, simulation_time; ybounds=[2,11], variable="elevation")
    colors = [get(ColorSchemes.oslo, i / 4) for i in 0:3]

    a = filter(row -> row.t == minimum(df.t),df)
    b = filter(row -> row.t == round(maximum(df.t) * 0.25), df)
    c = filter(row -> row.t == round(maximum(df.t) * 0.5), df)
    d = filter(row -> row.t == maximum(df.t),df)
    if variable == "elevation"
        # Initial Bed surface
        plot(a.x, a.η, xlim=(0,20), ylim=(-1,1), title=latexstring("\$η^*\$ Evolution", ), 
        color=colors[4], xlabel=latexstring("Dimensionless Distance, \$x^*\$"), ylabel=latexstring("Dimensionless Elevation, \$η^{*}(x^*,t^*)\$"),
        label=latexstring("\$t^{*} = 0\$"), linestyle=:dash, legend=:topleft)
        plot!(twiny(), xticks=false)

        plot!(b.x,b.η,label=latexstring("\$t^{*} = \$ $(scientific_notation(b[1,"t"]))"), color=colors[3], linewidth=4)
        plot!(c.x,c.η,label=latexstring("\$t^{*} = \$ $(scientific_notation(c[1,"t"]))"), color=colors[2], linewidth=4)
        plot!(d.x,d.η,label=latexstring("\$t^{*} = \$ $(scientific_notation(d[1,"t"]))"), color=colors[1], linewidth=4)
        p = plot!(twinx(), a.x, a.roughness, color="grey", label="Roughness", ylim=(ybounds[1],ybounds[2]), yticks=false)
        return p
    elseif variable == "velocity"
         # Initial Bed surface
         plot(a.x, a.u, xlim=(0,20), ylim=(0,2), title=latexstring("\$u^*\$ Evolution", ), 
         color=colors[4], xlabel=latexstring("Dimensionless Distance, \$x^*\$"), ylabel=latexstring("Dimensionless Velocity, \$u^{*}(x^*,t^*)\$"),
         label=latexstring("\$t^{*} = 0\$"), linestyle=:dash, legend=:topleft)
         plot!(twiny(), xticks=false)
 
         plot!(b.x,b.u,label=latexstring("\$t^{*} = \$ $(scientific_notation(b[1,"t"]))"), color=colors[3], linewidth=4)
         plot!(c.x,c.u,label=latexstring("\$t^{*} = \$ $(scientific_notation(c[1,"t"]))"), color=colors[2], linewidth=4)
         plot!(d.x,d.u,label=latexstring("\$t^{*} = \$ $(scientific_notation(d[1,"t"]))"), color=colors[1], linewidth=4)
         p = plot!(twinx(), a.x, a.roughness, color="grey", label="Roughness", ylim=(ybounds[1],ybounds[2]), yticks=false)
         return p
    else
        println("Error: argument 'variable' must be 'elevation' or 'velocity'.")
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
g = 9.81 # Gravitational Acceleration m/s^2
κ = 0.41 # Von Karman constant
ρ_s = 2650 # Sediment density
ρ = 1000 # Fluid density
γ = 1.5

# Reference angle of slope in radians (10% 0.09967, 8% 0.07983, 7% 0.06989 5% 0.04996)

"""
compute_dimensionless(q_w,λ,r_0,D)

Compute the dimenionless numbers of the system.

# Returns
- The dimenionless numbers in the order: ζ, α, β, ψ
"""
function compute_dimensionless(q_w,λ,r_0,D)
    ζ = (q_w^(2.0/3.0))/((g^(1.0/3.0))*r_0)
    α = (ρ*(g*q_w)^(2.0/3.0) * κ^2.0) / ((ρ_s - ρ)*g*D)
    β = (λ*κ^2.0)/(r_0)
    ψ = λ / r_0
    return ζ, α, β, ψ
end

function TuneTurbulence(D::Float64; Factor::Float64=3.0)
    return (30/(Factor*D))*D
end

function ComputeStress(u, r, α, ζ, ν)
    Numerator = α .* u .* ( (ζ^2.0 ./ (u.^2.0)) .- ((2.0 .* r .* ζ) ./ u) .+ r.^2.0)
    Denominator = (((ζ./u) .* (NaNMath.log.(1 .+ ((ν*ζ)./(u.*r))) .- 1)) .+ ((r./ν) .* (ν - ((1+ν) * NaNMath.log(1+ν)) .+ NaNMath.log.(1 .+ ((ν*ζ) ./ (u .* r)))))).^2.0
    return Numerator ./ Denominator
end

struct Parameters
    ζ::Float64
    β::Float64
    ψ::Float64
    ν::Float64
    r0::Float64
    θ::Float64
end

function equilibrium_u0!(f, u, p::Parameters)
    u0 = u[1]
    F = u0^3.0 * (( (p.ζ^2.0 / u0^2.0) - ((2.0 * p.r0 * p.ζ)/u0) + p.r0^2.0) 
    / ( ((p.ζ/u0)*(log(1.0 + ((p.ν*p.ζ)/(u0 * p.r0))) - 1.0)) 
    + ((p.r0/p.ν) * (p.ν - ((1+p.ν) * log(1.0 + p.ν)) 
    + log(1.0 + ((p.ν*p.ζ)/(u0*p.r0))))) )^(2.0))

    f[1] = p.β * F - p.ψ * tan(p.θ)
end

"""
EqulibriumVelocity(f, u_guess, p)

Compute the equilibrium velocity given a function `f`, an initial guess `u_guess`, 
and a vector of parameters `p`.

# Arguments
- `f`: The function to be used in the computation.
- `u_guess`: Initial guess for the equilibrium.
- `p`: A vector of parameters in the following order:
      `[ζ, β, ψ, ν, r0, θ]`

# Returns
- The computed equilibrium velocity.
"""
function EqulibriumVelocity(f, u_guess, p::Parameters)
    # p order: [ζ, β, ψ, ν, r0, θ]
    prob = NonlinearProblem(equilibrium_u0!, [u_guess], p)
    sol = solve(prob, NewtonRaphson())
    return sol[1]
end

# Numerical Solver
function Gaussian_DualPeak_Model(Peak1::Float64, Peak2::Float64, Discharge::Float64, Roughness_Scale::Float64, Length_Scale::Float64,
     D50::Float64, Slope::Float64, timestep::Float64, simtime::Float64, ν::Float64; 
     Peak3::Union{Nothing, Float64}=nothing, u_init::Float64=1.0, τ_cr::Float64=0.05, σ::Float64=0.05, A::Float64=2.0,
     Reference::String=nothing, figure::Bool=false, flat::Bool=false, k::Float64=3.0)

    q_w = Discharge # Discharge per width
    θ = Slope
    r_0 = Roughness_Scale
    D = D50
    λ = Length_Scale
    τ_cr = τ_cr # Critical Shields Number

    c = u_init + sqrt(q_w / u_init)

    ζ, α, β, ψ = compute_dimensionless(q_w,λ,r_0,D)
    
    @parameters x, t
    @variables u(..) [bounds=(0,Inf)]
    @variables η(..) 
    
    # Linear Operators
    Dx = Differential(x)
    Dt = Differential(t)
    
    # Gaussian Peaked Roughness
    if !flat
        peak1 = (A/r_0)*exp((-(x - Peak1)^2.0)/σ)
        peak2 = (A/r_0)*exp((-(x - Peak2)^2.0)/σ)
    
        if Peak3 !== nothing
            peak3 = (A/r_0)*exp((-(x - Peak3)^2.0)/σ)
            r = gau = (peak1 + peak2 + peak3 + (k*D/r_0))
        else
            r = gau = (peak1 + peak2 + (k*D/r_0))
        end
    else
        r = D/r_0
    end

    eq1 = ζ * u(x,t) * Dx(u(x,t)) ~ -ζ * expand_derivatives(Dx(1/u(x,t))) - Dx(η(x,t)) - ((β*u(x,t)^3.0 * ( (ζ^2 / (u(x,t))^2) - ((2.0 * r * ζ)/u(x,t)) + (r^(2.0)))) / ( ((ζ/u(x,t))*(NaNMath.log(1 + ((ν*ζ)/(u(x,t) * r))) - 1.0)) + ((r/ν) * (ν - ((1+ν) * NaNMath.log(1 + ν)) + NaNMath.log(1 + ((ν*ζ)/(u(x,t)*r))))) )^(2.0)) + ψ*tan(θ)
    
    eq2B = -expand_derivatives(Dx((α * u(x,t)^2.0) * (( (ζ^2 / (u(x,t))^2.0) - ((2.0 * r * ζ)/u(x,t)) + r^2.0) / ( ((ζ/u(x,t))*(NaNMath.log(1 + ((ν*ζ)/(u(x,t) * r))) - 1.0)) + ((r/ν) * (ν - ((1+ν) * NaNMath.log(1 + ν)) + NaNMath.log(1 + ((ν*ζ)/(u(x,t)*r))))) )^(2.0)) - τ_cr)) * (1/α)

    τ = (α * u(x,t)^2.0) * (( (ζ^2 / (u(x,t))^2) - ((2.0 * r * ζ)/u(x,t)) + r^2.0) / ( ((ζ/u(x,t))*(NaNMath.log(1 + ((ν*ζ)/(u(x,t) * r))) - 1.0)) + ((r/ν) * (ν - ((1+ν) * NaNMath.log(1 + ν)) + NaNMath.log(1 + ((ν*ζ)/(u(x,t)*r))))) )^(2.0))

    f(x,t) = ((τ >= τ_cr) * eq2B)^γ
    eq2 = Dt(η(x,t)) ~ f(x,t)
    
    eq = [eq1, eq2]

    # Domain
    x_start = t_start = 0.0
    x_end = 20.0
    t_end = simtime
    
    domains = [x ∈ IntervalDomain(x_start,x_end), t ∈ IntervalDomain(t_start,t_end)]

    η0(x, t) = 0.0
    u0(x, t) = u_init

    if flat
        # Periodic Boundary Conditions
        bcs = [η(x, 0.0) ~ η0(x, 0.0),
                η(x_start, t) ~ η(x_end, t),
                u(x_start,t) ~ u(x_end,t),
                u(x,0.0) ~ u0(x,0.0)]
    else
        #=
        bcs = [η(x, 0.0) ~ η0(x, 0.0),
                η(x_start, t) ~ η(x_end, t),
                u(x_start,t) ~ u(x_end,t),
                u(x,0.0) ~ u0(x,0.0)]
        =#
        bcs = [
            # Initial conditions for all x at t=0
            η(x, 0.0) ~ η0(x, 0.0),
            u(x, 0.0) ~ u0(x, 0.0),
        
            # Dirichlet at x_start
            u(x_start, t) ~ u_init,
            #η(x_start, t) ~ 0.0,
        
            # Sommerfeld
            Dx(u(x_end,t)) ~ (-1/c) * Dx(u(x_end,t)),
            #Dx(u(x_start,t)) ~ (-1/c) * Dx(u(x_start,t)),

            # Neumann
            #Dx(η(x_end,t)) ~ 0.0
            ]
    end
    
    @named pdesys = PDESystem(eq, bcs, domains, [x, t], [u(x, t), η(x, t)])

    discretization = MOLFiniteDifference([x => 250], t)

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
    spatial_indices = repeat(range(0.0, stop=20.0, length=rows), inner=cols)
    
    discr = (((A*D)/r_0) .* exp.((-(discx .- Peak1).^2.0)./σ)) .+ (((A*D)/r_0) .* exp.((-(discx .- Peak2).^2.0)./σ)) .+ (k*D/r_0)
    if Peak3 !== nothing
        discr .+= (((A*D)/r_0) .* exp.((-(discx .- Peak3).^2.0)./σ))
    end
    
    roughness = repeat(discr, inner=cols)

    #τ_values = [Symbolics.value(substitute(τ, Dict(x => i, t => j, u(x, t) => k, η(x,t) => m))) for (i, j, k, m) in zip(spatial_indices, time_slices, discu_vec, discη_vec)]

    df = DataFrame(
        t = time_slices,
        x = spatial_indices,
        u = discu_vec,
        η = discη_vec,
        roughness = roughness,
        #tau = τ_values,
        q = q_w,
        slope = θ,
        r0 = r_0,
        ν = ν,
        Reference = Reference
        )

    τ_values = ComputeStress(df.u, df.roughness, α, ζ, ν)
    df.tau = vec(τ_values')
    
    df[:, :Peak1].=Peak1
    df[:, :Peak2].=Peak2
    df[:, :Peak3] .= Peak3 === nothing ? NaN : Peak3

    if figure
        p = PlotEvolution(disct, discx, discη, r, simtime, timestep)
        return df, p
    else
        return df
    end
end

function process_simulations(data::DataFrame)
    # Step 1: Create 'Run' as a unique group identifier based on the combination of columns
    data.Run = map(x -> string(x), groupby(data, [:r0, :slope, :q, :Peak1, :Peak2, :Peak3]).groups)
    data.Run = categorical(data.Run)
    data = filter(row -> 5.0 < row.x < 15.0, data)

    # Step 2: Group data by 'Run' and 't' and calculate 'mean_u'
    grouped = groupby(data, [:Run, :t])
    mean_u_df = transform(grouped, [:u, :q] => ((u, q) -> mean(u .* (9.81 .* q) .^ (1/3))) => :mean_u)

    # Step 3: Ungroup to proceed with further transformations
    ungrouped = DataFrame(mean_u_df)  # Convert GroupedDataFrame to a normal DataFrame

    # Step 4: Calculate 'sigma' as the standard deviation of 'η'
    sigma_df = transform(ungrouped, :η => (η -> std(η)) => :sigma)

    # Step 5: Calculate 'h' as q / mean_u
    h_df = transform(sigma_df, [:q, :mean_u] => ((q, mean_u) -> q ./ mean_u) => :h)

    # Step 6: Calculate 'SigmaRatio' as h / sigma and 'DRatio' as h / r0
    ratio_df = transform(h_df, [:h, :sigma] => ((h, sigma) -> h ./ sigma) => :SigmaRatio,
                         [:h, :r0] => ((h, r0) -> h ./ r0) => :DRatio)

    # Step 7: Filter rows where t > 0
    filtered_df = filter(:t => t -> t > 0, ratio_df)

    # Step 8: Calculate 'Friction' as mean_u / sqrt(h * 9.81 * slope)
    friction_df = transform(filtered_df, [:mean_u, :h, :slope] => ((mean_u, h, slope) -> mean_u ./ sqrt.(h .* 9.81 .* slope)) => :Friction)

    return friction_df
end

function ExclusionZone(Data::DataFrame, H::Float64, slope::Float64)
    data = process_simulations(Data)
    final_state = filter(r -> r.t == maximum(r.t), data)
    f = (Statistics.mean(final_state.Friction))^(-2)
    U_0 = mean(final_state.mean_u)
    Γ = U_0^2 / (9.81 * sin(slope))
    ϵ = (H/Γ)^(1/3)
    min_length = (ϵ / f)*H
    return min_length
end

function run_simulations(peaks::Vector{Tuple{Float64, Float64, Union{Nothing, Float64}}}, Discharge::Float64, Roughness_Scale::Float64, Length_Scale::Float64, D50::Float64, Slope::Float64, timestep::Float64, simtime::Float64, ν; u_init::Float64=1.0, τ_cr::Float64=0.05, Reference::String=nothing)
    results = [Gaussian_DualPeak_Model(Peak1, Peak2, Discharge, Roughness_Scale, Length_Scale, D50, Slope, timestep, simtime, ν, u_init=u_init, τ_cr=τ_cr, Peak3=Peak3, Reference=Reference) for (Peak1, Peak2, Peak3) in peaks]
    return vcat(results...)
end

function create_animations_list(simulation_results, peaks, output_directory, fps)
    for i in peaks
        anim = create_animations(simulation_results, i)
        gif(anim, string(output_directory,"/simulation_results_$(i[1])_$(i[2])_$(i[3] === nothing ? "none" : string(i[3]))_peak.gif"), fps=fps)
    end
end
