# Step 1: Load the model as @reaction_network from Catalyst.
using Catalyst
include("../src/rxns.jl")
rn = get_rxns()

# Step 2: Convert the model to an ODESystem and simplify it.
rn_odesys = structural_simplify(convert(ODESystem,rn;combinatoric_ratelaws=false))
psym = parameters(rn_odesys)

# Step 3: Initialize the inital conditions and parameters of the model from a CSV file.
using CSV, DataFrames
veryhitgcsv = "./log/veryhitg.csv"
csvfile = veryhitgcsv
df = CSV.read(csvfile,DataFrame)
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
# Doublicate for ODEProblem.
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


# Step 4: Solve ODE.
using DiffEqCallbacks, OrdinaryDiffEq
tspan = (0.0,10000.0) # [days]
cb = TerminateSteadyState(1e-8,1e-6) # callback to stop at steadystate, more efficient for M-H search
op = ODEProblem(rn_odesys, u0, tspan, p0, callback=cb)
sol = solve(op, Rodas5())

# Step 5: Visualize solution.
using Plots
p = plot(sol, legend = :outerright)
