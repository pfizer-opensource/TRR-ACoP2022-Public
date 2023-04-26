#=
Helper function for the M-H search.

This function needs to know how to take an input parameter vector and ODE function
    and return a score the M-H algorithm to make its next move on.
=#

function mh_sim(pin::Vector,OdeF::ODEProblem,
    xlb::Vector,xub::Vector,d_obs::Distribution,
    ctpidx::Int64, pk9idx::Int64, vd_lcy_idx::Int64, vd_ler_idx::Int64)

    OdeFn = OdeF
    OdeFn = remake(OdeFn,p=pin)#copy(p)

    sol = solve(OdeFn,Rodas5()) # Simulate the model, the timing etc. should be pre-definied in OdeF
    sim_obs = GetObs(sol,pin,vd_lcy_idx,vd_ler_idx) # Get the observables from the solution construct

    OdeCTP = OdeFn

    pctp = copy(pin)
    pctp[ctpidx] = 0.3

    OdeCTP = remake(OdeCTP,tspan=(0.0,28.0),u0=sol[end],p=pctp)
    solCTP = solve(OdeCTP,Rodas5())

    tc_pcb = 100.0*(solCTP[ch_a1+ch_a100+ch_a48][end]/solCTP[ch_a1+ch_a100+ch_a48][1] - 1.0)
    hdl_pcb = 100.0*(solCTP[ch_a1][end]/solCTP[ch_a1][1] - 1.0)
    tg_pcb = 100.0*(solCTP[tg_a1+tg_a100+tg_a48][end]/solCTP[tg_a1+tg_a100+tg_a48][1] - 1.0)
    ctp_err = (tc_pcb - 0.48)^2 + (hdl_pcb - 22.)^2 + (tg_pcb + 4.0)^2 # DeGrooth 2002, 600mg, 4wks

    ppk9 = copy(pin)
    ppk9[pk9idx] = 3.5 # Shooting for -60 --> -83% inhibition
    OdePK9 = remake(OdeFn, tspan = (0.0,84.0), u0=sol[end],p=ppk9)
    solPK9 = solve(OdePK9,Rodas5())
    tc_pcb = 100.0*(solPK9[ch_a1+ch_a100+ch_a48][end]/solPK9[ch_a1+ch_a100+ch_a48][1] - 1.0)
    hdl_pcb = 100.0*(solPK9[ch_a1][end]/solPK9[ch_a1][1] - 1.0)
    tg_pcb = 100.0*(solPK9[tg_a1+tg_a100+tg_a48][end]/solPK9[tg_a1+tg_a100+tg_a48][1] - 1.0)
    pk9_err = (tc_pcb + 40.)^2 + (hdl_pcb - 7.)^2 + (tg_pcb - 17.0)^2 # Zhang et al. BMC Medicine. 2015.

    # Score the solution, if any X is out-of-bounds, default to eps()
    if all(sol[:,end] .>= xlb) && all(sol[:,end] .<= xub)
        q = pdf(d_obs,sim_obs)*exp(-0.01*ctp_err)*exp(-0.01*pk9_err) # slope CETPi/PK9 Antibody exp() was by feel
    else
        q = eps()
    end

    return q
end

function GetObs(sol::ODESolution, p::Vector{Float64}, vd_lcy_idx::Int64, vd_ler_idx::Int64)::Vector{Float64}
    # Helper function that extracts the relevant part of the solution structure, simple in our case.
    #sim_obs = [Float64.(getLTGprct(sol[:,end],p)[1]),sol[5,end]]
    lnHDL = log.(sol[ch_a1][end])
    lnTC = log.(sol[ch_a1 + ch_a100 + ch_a48][end])
    lnTG = log.(sol[tg_a1 + tg_a100 + tg_a48][end])

    LF = sol[tg_lcy*p[vd_lcy_idx] + tg_ler*p[vd_ler_idx]][end].* 885.0/1000.0 / 1500.0 * 100.0 # [%]

    sim_obs = [lnHDL;lnTC;lnTG;log(LF)]
    
    return sim_obs
end
