# If not already installed, first hit "]" within a Julia REPL. Then type:
# add Catalyst DifferentialEquations Plots Latexify

using Catalyst
using DifferentialEquations
using Latexify
using OrdinaryDiffEq
using ModelingToolkit
using DataFrames
using CSV
using DiffEqCallbacks
using Distributions
using CairoMakie
using FileIO
using JLD2

include("rxns.jl")
include("drawp.jl")
include("readnhanes.jl")
include("mh.jl")
include("mh_sim.jl")
#include("mh_sim_limited.jl")
include("select_vps.jl")
indexof(sym,syms) = findfirst(isequal(sym),syms)

JLD2_PP_FILE = "./log/midD.jld2"

highTGcsv = "./log/parameters-highTG.csv"
midTGcsv = "./log/parameters.csv"
midDcsv = "./log/parameters-median-dist.csv"
veryhitgcsv = "./log/veryhitg.csv"

csvfile = veryhitgcsv

df = CSV.read(csvfile,DataFrame)

rn = get_rxns()

rn_odesys = structural_simplify(convert(ODESystem,rn;combinatoric_ratelaws=false))
psym = parameters(rn_odesys)

@variables begin
    t
    a100(t) = df.TV[df.Parameter .== "a100"][1]
    tg_a100(t) = df.TV[df.Parameter .== "tg_a100"][1]
    ch_a100(t) = df.TV[df.Parameter .== "ch_a100"][1]
    a1(t) = df.TV[df.Parameter .== "a1"][1]
    tg_a1(t) = df.TV[df.Parameter .== "tg_a1"][1]
    ch_a1(t) = df.TV[df.Parameter .== "ch_a1"][1]
    a48(t) = df.TV[df.Parameter .== "a48"][1]
    tg_a48(t) = df.TV[df.Parameter .== "tg_a48"][1]
    ch_a48(t) = df.TV[df.Parameter .== "ch_a48"][1]
    fa_lcy(t) = df.TV[df.Parameter .== "fa_lcy"][1]
    fa_ler(t) = df.TV[df.Parameter .== "fa_ler"][1]
    tg_lcy(t) = df.TV[df.Parameter .== "tg_lcy"][1]
    tg_ler(t) = df.TV[df.Parameter .== "tg_ler"][1]
    ch_l(t) = df.TV[df.Parameter .== "ch_l"][1]
    ldlr(t) = df.TV[df.Parameter .== "ldlr"][1]
    pk9(t) = df.TV[df.Parameter .== "pk9"][1]
    ctpi_gut(t) = df.TV[df.Parameter .== "ctpi_gut"][1]
    ctpi_cent(t) = df.TV[df.Parameter .== "ctpi_cent"][1]
    ctpi_Q1(t) =  df.TV[df.Parameter .== "ctpi_Q1"][1]
end

@parameters begin
    ks_a100 = df.TV[df.Parameter .== "ks_a100"][1]
    ks_a1 = df.TV[df.Parameter .== "ks_a1"][1]
    kcl_a1 = df.TV[df.Parameter .== "kcl_a1"][1]
    krct = df.TV[df.Parameter .== "krct"][1]
    ksrb1 = df.TV[df.Parameter .== "ksrb1"][1]
    klpl = df.TV[df.Parameter .== "klpl"][1]
    kldlr_a100 = df.TV[df.Parameter .== "kldlr_a100"][1]
    kldlr_a48 = df.TV[df.Parameter .== "kldlr_a48"][1]
    ks_a48 = df.TV[df.Parameter .== "ks_a48"][1]
    kctp = df.TV[df.Parameter .== "kctp"][1]
    ks_fa = df.TV[df.Parameter .== "ks_fa"][1]
    kdnl = df.TV[df.Parameter .== "kdnl"][1]
    kd_fa = df.TV[df.Parameter .== "kd_fa"][1]
    kest_lcy = df.TV[df.Parameter .== "kest_lcy"][1]
    kest_ler = df.TV[df.Parameter .== "kest_ler"][1]
    klip = df.TV[df.Parameter .== "klip"][1]
    ker = df.TV[df.Parameter .== "ker"][1]
    ks_ch = df.TV[df.Parameter .== "ks_ch"][1]
    kd_ch = df.TV[df.Parameter .== "kd_ch"][1]
    ks_ldlr = df.TV[df.Parameter .== "ks_ldlr"][1]
    kd_ldlr = df.TV[df.Parameter .== "kd_ldlr"][1]
    ks_pk9 = df.TV[df.Parameter .== "ks_pk9"][1]
    kcl_pk9 = df.TV[df.Parameter .== "kcl_pk9"][1]
    vd_ler = df.TV[df.Parameter .== "vd_ler"][1]
    vd_lcy = df.TV[df.Parameter .== "vd_lcy"][1]
    vd_p = df.TV[df.Parameter .== "vd_p"][1]
    f_lpl_h = df.TV[df.Parameter .== "f_lpl_h"][1]
    alpha_a48 = df.TV[df.Parameter .== "alpha_a48"][1]
    beta_a48 = df.TV[df.Parameter .== "beta_a48"][1]
    alpha_a100 = df.TV[df.Parameter .== "alpha_a100"][1]
    beta_a100 = df.TV[df.Parameter .== "beta_a100"][1]
    f_ins = df.TV[df.Parameter .== "f_ins"][1]
    K_tgler = df.TV[df.Parameter .== "K_tgler"][1]
    K_chl = df.TV[df.Parameter .== "K_chl"][1]
    cetp_scale = df.TV[df.Parameter .== "cetp_scale"][1]
    pk9_scale = df.TV[df.Parameter .== "pk9_scale"][1]
    ka_ctpi = df.TV[df.Parameter .== "ka_ctpi"][1]
    kel_ctpi = df.TV[df.Parameter .== "kel_ctpi"][1]
    vd_ctpi_cent = df.TV[df.Parameter .== "vd_ctpi_cent"][1]
    vd_ctpi_Q1 = df.TV[df.Parameter .== "vd_ctpi_Q1"][1]
    Q_ctpi = df.TV[df.Parameter .== "Q_ctpi"][1]
end


p0 = [
    ks_a100 => df.TV[df.Parameter .== "ks_a100"][1],
    ks_a1 => df.TV[df.Parameter .== "ks_a1"][1],
    kcl_a1 => df.TV[df.Parameter .== "kcl_a1"][1],
    krct => df.TV[df.Parameter .== "krct"][1],
    ksrb1 => df.TV[df.Parameter .== "ksrb1"][1],
    klpl => df.TV[df.Parameter .== "klpl"][1],
    kldlr_a100 => df.TV[df.Parameter .== "kldlr_a100"][1],
    kldlr_a48 => df.TV[df.Parameter .== "kldlr_a48"][1],
    ks_a48 => df.TV[df.Parameter .== "ks_a48"][1],
    kctp => df.TV[df.Parameter .== "kctp"][1],
    ks_fa => df.TV[df.Parameter .== "ks_fa"][1],
    kdnl => df.TV[df.Parameter .== "kdnl"][1],
    kd_fa => df.TV[df.Parameter .== "kd_fa"][1],
    kest_lcy => df.TV[df.Parameter .== "kest_lcy"][1],
    kest_ler => df.TV[df.Parameter .== "kest_ler"][1],
    klip => df.TV[df.Parameter .== "klip"][1],
    ker => df.TV[df.Parameter .== "ker"][1],
    ks_ch => df.TV[df.Parameter .== "ks_ch"][1],
    kd_ch => df.TV[df.Parameter .== "kd_ch"][1],
    ks_ldlr => df.TV[df.Parameter .== "ks_ldlr"][1],
    kd_ldlr => df.TV[df.Parameter .== "kd_ldlr"][1],
    ks_pk9 => df.TV[df.Parameter .== "ks_pk9"][1],
    kcl_pk9 => df.TV[df.Parameter .== "kcl_pk9"][1],
    vd_ler => df.TV[df.Parameter .== "vd_ler"][1],
    vd_lcy => df.TV[df.Parameter .== "vd_lcy"][1],
    vd_p => df.TV[df.Parameter .== "vd_p"][1],
    f_lpl_h => df.TV[df.Parameter .== "f_lpl_h"][1],
    alpha_a48 => df.TV[df.Parameter .== "alpha_a48"][1],
    beta_a48 => df.TV[df.Parameter .== "beta_a48"][1],
    alpha_a100 => df.TV[df.Parameter .== "alpha_a100"][1],
    beta_a100 => df.TV[df.Parameter .== "beta_a100"][1],
    f_ins => df.TV[df.Parameter .== "f_ins"][1],
    K_tgler => df.TV[df.Parameter .== "K_tgler"][1],
    K_chl => df.TV[df.Parameter .== "K_chl"][1],
    cetp_scale => df.TV[df.Parameter .== "cetp_scale"][1],
    pk9_scale => df.TV[df.Parameter .== "pk9_scale"][1],
    ka_ctpi => df.TV[df.Parameter .== "ka_ctpi"][1],
    kel_ctpi => df.TV[df.Parameter .== "kel_ctpi"][1],
    vd_ctpi_cent => df.TV[df.Parameter .== "vd_ctpi_cent"][1],
    vd_ctpi_Q1 => df.TV[df.Parameter .== "vd_ctpi_Q1"][1],
    Q_ctpi => df.TV[df.Parameter .== "Q_ctpi"][1]
    ]

u0 = [
    a100 => df.TV[df.Parameter .== "a100"][1],
    tg_a100 => df.TV[df.Parameter .== "tg_a100"][1],
    ch_a100 => df.TV[df.Parameter .== "ch_a100"][1],
    a1 => df.TV[df.Parameter .== "a1"][1],
    tg_a1 => df.TV[df.Parameter .== "tg_a1"][1],
    ch_a1 => df.TV[df.Parameter .== "ch_a1"][1],
    a48 => df.TV[df.Parameter .== "a48"][1],
    tg_a48 => df.TV[df.Parameter .== "tg_a48"][1],
    ch_a48 => df.TV[df.Parameter .== "ch_a48"][1],
    fa_lcy => df.TV[df.Parameter .== "fa_lcy"][1],
    fa_ler => df.TV[df.Parameter .== "fa_ler"][1],
    tg_lcy => df.TV[df.Parameter .== "tg_lcy"][1],
    tg_ler => df.TV[df.Parameter .== "tg_ler"][1],
    ch_l => df.TV[df.Parameter .== "ch_l"][1],
    ldlr => df.TV[df.Parameter .== "ldlr"][1],
    pk9 => df.TV[df.Parameter .== "pk9"][1],
    ctpi_gut => df.TV[df.Parameter .== "ctpi_gut"][1],
    ctpi_cent => df.TV[df.Parameter .== "ctpi_cent"][1],
    ctpi_Q1 =>  df.TV[df.Parameter .== "ctpi_Q1"][1]
]

# Solve to baseline:
tspan = (0.0,10000.0) # [days]
cb = TerminateSteadyState(1e-8,1e-6) # callback to stop at steadystate, more efficient for M-H search
op = ODEProblem(rn_odesys, u0, tspan, p0, callback=cb)
sol = solve(op, Rodas5())

pfit_name = [ks_a100, ks_a1, kcl_a1,
krct, ksrb1, klpl, kldlr_a100,
kldlr_a48, ks_a48, kdnl, kd_fa,
kest_lcy, kest_ler, klip,
ker, ks_ch, kd_ch, ks_ldlr, kd_ldlr,
ks_pk9, kcl_pk9, f_lpl_h, alpha_a48,
beta_a48, alpha_a100, beta_a100, f_ins]

fiti = Vector{Int64}(undef,size(pfit_name)[1])
for (idx, val) in enumerate(pfit_name)
    fiti[idx] = indexof(val,psym)
end

#=
function mh(n_pp::Integer, n_param_fit::Integer, log_flag::Bool,
        p::Vector{Float64}, plb::Vector{Float64},pub::Vector{Float64}, 
        Sim_and_Score::Function, OdeF::Function)
=#


# Options for S.S. MH search:
n_param_fit = size(pfit_name)[1]
log_flag = false
vd_lcy_idx = indexof(vd_lcy, psym)
vd_ler_idx = indexof(vd_ler, psym)
p = copy(op.p)
plb = 0.2 .* p
pub = 4.0 .* p
OdeF = op
#Sim_and_Score = 
d_obs = return_nhanes()
xlb = sol[:,end] .* 0.1
xub = sol[:,end] .* 10.0
ctpidx = indexof(cetp_scale,psym)
pk9idx = indexof(pk9_scale,psym)
FScore(a,opn) = mh_sim(a,opn,xlb,xub,d_obs,ctpidx, pk9idx, vd_lcy_idx, vd_ler_idx)

EofX = p[fiti]
StdofX = (pub[fiti] .- plb[fiti])./8.0 
VarofX = StdofX.^2

mu = log.(EofX.^2 ./ sqrt.(VarofX .+ EofX.^2))
sigma = log.(1. .+ VarofX./EofX.^2) .* I(n_param_fit)

#sigmaX = varX.*I(n_param_fit)
d_param = MvNormal(mu,sigma)

#= 
function drawp(p0::Vector{Float64}, dparam::Distribution,
    fiti::Vector{Int64}, log_flag::Bool)::Vector{Float64}
    p = Vector{Float64}(undef,size(p0)[1])
    =#
DrawF(pn) = drawp(pn, d_param, fiti, true)

# Plausible Patient search (commented out by default). This step is slow.
# Each call to 'mh' launches a Metropolis-Hastings seach for new PPs.
# Multiple starts were needed to well-populate the 4D distribution.
#=
p_pp_midD = mh(160, p, FScore, DrawF, OdeF)
#FileIO.save("./log/midD-PP.jld2","p_pp_midD",p_pp_midD)

p_pp_hitg = mh(160, p, FScore, DrawF, OdeF)
FileIO.save("./log/highTG-PP.jld2","p_pp_hitg",p_pp_hitg)

p_pp_midtg = mh(16000, p, FScore, DrawF, OdeF)
FileIO.save("./log/midTG-PP.jld2","p_pp_midtg",p_pp_midtg)

p_pp_vhitg = mh(8000, p, FScore, DrawF, OdeF)
FileIO.save("./log/vhiTG-PP.jld2","p_pp_vhitg",p_pp_vhitg)
=#

# Load pre-run plausible patient (PP) files. Each file represents a different
# prior seeding a Metropolis-Hastings simluation.
p_pp_hitg = FileIO.load("./log/highTG-PP.jld2","p_pp_hitg")
p_pp_midtg = FileIO.load("./log/midTG-PP.jld2","p_pp_midtg")
p_pp_midD = FileIO.load("./log/midD-PP.jld2","p_pp_midD")
p_pp_vhitg = FileIO.load("./log/vhiTG-PP.jld2","p_pp_vhitg")

# Combine the various PP runs into one larger PP population:
p_pp = [p_pp_midtg p_pp_hitg p_pp_midD p_pp_vhitg] 

# Pre-define vectors of observables:
pn = Vector{Float64}(undef,size(p)[1]) # Number of parameters
obs_pp = Matrix{Float64}(undef,4,size(p_pp)[2]) # Matrix of observables
ss_pp = Matrix{Float64}(undef,size(u0)[1],size(p_pp)[2]) # Matrix of S.S. values for PPs
lf_prct = Vector{Float64}(undef,size(p_pp)[2]) # LF% for each PP

# Re-simulate every PP, save their "observable" states for Virtual Patient (VP) selection:
for ii = 1:size(p_pp)[2]
    pn .= copy(p) # Playing it safe not to overwrite the original PP population
    opn = op
    opn.p .= p_pp[:,ii]
    soln  = solve(opn, Rodas5())
    obs_pp[1:4,ii] = GetObs(soln,p_pp[:,ii],vd_lcy_idx,vd_ler_idx) #exp.(GetObs(soln))

	# Save the steady state of each PP (starting point for therapy simulations):
    ss_pp[:,ii] .= soln[end]

    lf_end = soln[tg_lcy*p_pp[vd_lcy_idx,ii] + tg_ler*p_pp[vd_ler_idx,ii]][end] # mmols-TG, total fat in hepatocyte

    lf_prct[ii] = lf_end .* 885.0/1000.0 / 1500.0 * 100.0 # g-TG/g-LW --> % LF # TODO: This calculation might be unnneeded now

end

# Select Virtual Patients from the Plausible Population using the ARM Method (Allen et al. CPT:PSP. 2016):
selidx, ppincld = select_vps(d_obs,obs_pp')

# CETPi simulations, we simulate a 600mg QD dose using a fixed inhibition of CETP for 4 weeks:

tc_pcb_cetpi = Vector{Float64}(undef,size(p_pp)[2]) 
hdl_pcb_cetpi = Vector{Float64}(undef,size(p_pp)[2])
tg_pcb_cetpi = Vector{Float64}(undef,size(p_pp)[2])
lf_pcb_cetpi = Vector{Float64}(undef,size(p_pp)[2])  
for ii = 1:size(p_pp)[2]
    pn .= copy(p)
    opn = op
    opn.p .= p_pp[:,ii]
    opn.p[ctpidx] = 0.45
    op_cepti = remake(opn,u0=ss_pp[:,ii],tspan=(0.0, 28.0))
    sol_cetpi = solve(op_cepti, Rodas5())

    tc_pcb_cetpi[ii] = 100. * (sol_cetpi[ch_a1+ch_a100+ch_a48][end]/sol_cetpi[ch_a1+ch_a100+ch_a48][1] - 1.0)
    hdl_pcb_cetpi[ii] = 100. * (sol_cetpi[ch_a1][end]/sol_cetpi[ch_a1][1] - 1.0)
    tg_pcb_cetpi[ii] = 100. * (sol_cetpi[tg_a1+tg_a100+tg_a48][end]/sol_cetpi[tg_a1+tg_a100+tg_a48][1] - 1.0)

    lf_base = sol_cetpi[tg_lcy*p_pp[vd_lcy_idx,ii] + tg_ler*p_pp[vd_ler_idx,ii]][1]
    lf_final = sol_cetpi[tg_lcy*p_pp[vd_lcy_idx,ii] + tg_ler*p_pp[vd_ler_idx,ii]][end]

    lf_pcb_cetpi[ii] = 100.0 * (lf_final/lf_base - 1.0)

end

#PCSK9 Antibody, we simulate a fixed inhibition of PCSK9 for 12 weeks:

tc_pcb_pk9ab = Vector{Float64}(undef,size(p_pp)[2]) 
hdl_pcb_pk9ab = Vector{Float64}(undef,size(p_pp)[2])
tg_pcb_pk9ab = Vector{Float64}(undef,size(p_pp)[2])
lf_pcb_pk9ab = Vector{Float64}(undef,size(p_pp)[2]) 
pk9_pcb_pk9ab = Vector{Float64}(undef,size(p_pp)[2]) 

# Loop through each Plausible Patient and add on the PCSK9 antibody:
for ii = 1:size(p_pp)[2]
    pn = p_pp[:,ii]
    pn[pk9idx] = 5 # 60% reduction in PCSK9 at week 12 (fasting)
    op_pk9ab = remake(op,p = pn, u0=ss_pp[:,ii],tspan=(0.0, 84.0)) 
    
    
    # ODE solve with the antibody:
    sol_pk9ab = solve(op_pk9ab, Rodas5())
    
	# Express key outputs as % change from baseline (PCSK9, PTCh, HDL-Ch, PTG):
    pk9_pcb_pk9ab[ii] = 100.0 * (sol_pk9ab[pk9][end]/sol_pk9ab[pk9][1] - 1.0)
    tc_pcb_pk9ab[ii] = 100. * (sol_pk9ab[ch_a1+ch_a100+ch_a48][end]/sol_pk9ab[ch_a1+ch_a100+ch_a48][1] - 1.0)
    hdl_pcb_pk9ab[ii] = 100. * (sol_pk9ab[ch_a1][end]/sol_pk9ab[ch_a1][1] - 1.0)
    tg_pcb_pk9ab[ii] = 100. * (sol_pk9ab[tg_a1+tg_a100+tg_a48][end]/sol_pk9ab[tg_a1+tg_a100+tg_a48][1] - 1.0)

    lf_base = sol_pk9ab[tg_lcy*p_pp[vd_lcy_idx,ii] + tg_ler*p_pp[vd_ler_idx,ii]][1]
    lf_final = sol_pk9ab[tg_lcy*p_pp[vd_lcy_idx,ii] + tg_ler*p_pp[vd_ler_idx,ii]][end]
	
	# Save the % change in liver fat:
    lf_pcb_pk9ab[ii] = 100.0 * (lf_final/lf_base - 1.0)
end


# Sensitivity Analaysis
# We perform a simple linearized SA by stepping up and down key parameters from baseline by a fixed
# percentage and simulating back to SS on each Virtual Patient. We then scale the outputs for visualization as a heatmap
# on a single scale: 

# Select a subset of parameters for the SA, the list below should be all of the biological parameters of the model
# the would naturally vary between individuals, some are "druggable", some are unlikely targets:
pSA = [ks_a100,ks_a1,kcl_a1,krct,ksrb1,klpl,kldlr_a100,
kldlr_a48,ks_a48,kctp,ks_fa,kdnl,kd_fa,kest_lcy,kest_ler,klip,
ker,ks_ch,kd_ch,ks_ldlr,kd_ldlr,ks_pk9,kcl_pk9,f_lpl_h,alpha_a48,
beta_a48,alpha_a100,beta_a100,f_ins]


# Labels for a heatmap for the above, this is only for plotting purposes and having a readable y-axis:
sdSA = ["ApoB100 Synth.","ApoA1 Synth.","ApoA1 CL","ABCA1 Activity","SRB1 Activity","LPL Activity","LDLR-Ab100",
"LDLR-Ab48","ApoB48 Synth.","CETP","LFA Uptake","DNL","LFA Oxid.","EstCyt","EstER",
"Lipol","ER Uptake of TG","Ch Synth.","Ch CL","LDLR Synth","LDLR CL","PCSK9 Synth.","PCSK9 CL",
"fLPL-Hepat","TG Diet","Ch Diet","TG VLDL","LMTP Activity","Ins Sens"]

# Find the index of each of the parameters, this is necessary because ModelingToolkit doesn't keep a consistent
# ordering of parametesrs:
SA_np = size(pSA)[1]
pSA_idx = Vector{Int64}(undef,SA_np)
for (ii,val) in enumerate(pSA)
    pSA_idx[ii] = indexof(val,psym)
end

# Fold change for each parameter from baseline, 1/step_up is the left-hand side:
step_up = 1.1

# Pre-declare some needed outputs:
p_pp_sel = p_pp[:,selidx]
pp_sel = ss_pp[:,selidx]
sa_med = Matrix{Float64}(undef,SA_np,4)

# Loop through each parameter above for the SA:
for ii = 1:SA_np
    sa_pp = Matrix{Float64}(undef,size(p_pp_sel)[2],4)
    
    # Loop through each Virtual Patient for each parameter in the SA:
    for jj = 1:size(p_pp_sel)[2]
    
    	# Save the original parameters and step up and down:
        pn2 = p_pp_sel[:,jj]
        pn2[pSA_idx[ii]] = step_up*p_pp_sel[pSA_idx[ii],jj]
        pn1 = p_pp_sel[:,jj]
        pn1[pSA_idx[ii]] = 1.0/step_up*p_pp_sel[pSA_idx[ii],jj]
		
		# Simulate the step up:
        opn = remake(op,p=pn2,u0=pp_sel[:,jj])
        soln2  = solve(opn, Rodas5())
        simobs2 = exp.(GetObs(soln2,pn2,vd_lcy_idx,vd_ler_idx))

		# Simulate the step down:
        opn = remake(op,p=pn1,u0=pp_sel[:,jj])
        soln1  = solve(opn, Rodas5())
        simobs1 = exp.(GetObs(soln1,pn1,vd_lcy_idx,vd_ler_idx))
		
		# Simulate the base case:
   	    pn0 = p_pp_sel[:,jj]
        remake(opn,p = pn0, u0 = pp_sel[:,jj])
        soln0 = solve(opn,Rodas5())
        simobs0 = exp.(GetObs(soln0,pn0,vd_lcy_idx,vd_ler_idx))
		
		# Express the SA outputs (the key observables) as d(ln(obs))/d(lnp):
        sa_pp[jj,:] = @. ((simobs2 - simobs1)/simobs0)/(step_up - 1.0/step_up)

    end
    # Take the median value of the population for reporting purposes, other metrics may
    # be equally (or more) valid here. We selected this one for the poster since it
    # seemed reflective of the population:
    sa_med[ii,1:4] = median(sa_pp, dims = 1)
end

# Scale the SA for reporting and plotting:
sa_min = minimum(sa_med,dims=1)
sa_max = maximum(sa_med,dims=1)

hm_scal = Matrix{Float64}(undef,SA_np,4) # Pre-allocate the matrix, rows = number of parameters, cols = number of observables

# Re-scale the medians to be between 0-->1:
hm_scal = (sa_med .- sa_min)./(sa_max .- sa_min)
std_hm_scal = std(hm_scal, dims = 1)
mean_hm_scal = mean(hm_scal, dims = 1)

# Find the "standout" values for a heatmap plot, this is just for display
# purposes of showing the top few rows:
hm_scal_scal = (hm_scal .- mean_hm_scal)./std_hm_scal
scal_r = Vector{Bool}(undef,size(hm_scal_scal)[1])
for ii = 1:size(hm_scal_scal)[1]
	# We set a cutoff here for plotting, 1.2 was selected just to show ~the top ten on the poster
	# but does not have any special significance beyond there:
    if ((maximum(hm_scal_scal[ii,:]) >= 1.2) || (minimum(hm_scal_scal[ii,:]) <= -1.2))
        maximum(hm_scal_scal[ii,:])
        minimum(hm_scal_scal[ii,:])
        scal_r[ii] = true
    else
        scal_r[ii] = false
    end
end

# Sub select the group of scaled SAs that met our threshold above, this should now be
# ~10 rows:
hm_scal_plot = hm_scal[scal_r,:]
hm_scal_plot

# Plotting:
include("figure1to3.jl")