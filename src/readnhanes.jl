using ReadStatTables
using DataFrames
using CairoMakie
using Distributions
using StatsBase
using CSV

function lognstat(m,v)
    mu = log(m^2/(sqrt(v+m^2)))
    sigma = sqrt(log(v/m^2+1))

    return (mu,sigma)
 end
 
function joint_tg_dist(frac_fast_tg::Float64=0.65)::Distribution
    # Creates a MvLogNormal distribution from the data below based on the mean and covariance of the data.
   
    # Data from captured from:
    # Anna Kotronen, Jukka Westerbacka, Robert Bergholm, Kirsi H. Pietiläinen, Hannele Yki-Järvinen, 
    # Liver Fat in the Metabolic Syndrome, 
    # The Journal of Clinical Endocrinology & Metabolism, 
    # Volume 92, Issue 9, 1 September 2007, Pages 3490–3497, 
    # https://doi.org/10.1210/jc.2007-0482
    df = CSV.read("./log/kotronen-liver-fat.csv",DataFrame)

    #x[:,2] = x[:,2]./1#frac_fast_tg
    #xl = log.(x)
    #mu = vec(mean(xl,dims=1))
    #cov = StatsBase.cov(xl)
    d_joint = fit(MvNormal,log.([df.LTG_PRCT';df.PTG_mM'])) # Multivariate log-normal distribution

    return d_joint    
end


function return_nhanes(dims=[1;2;3;4])

dfHDL = DataFrame(readstat("public/P_HDL.XPT"))
dfTG =  DataFrame(readstat("public/P_TRIGLY.XPT"))
dfCh = DataFrame(readstat("public/P_TCHOL.XPT"))

df = dropmissing(innerjoin(dfHDL,dfTG,dfCh, on=:SEQN))


dltg = joint_tg_dist()
dN = fit(MvNormal,[log.(df.LBDHDDSI)';log.(df.LBDTCSI)';log.(df.LBDTRSI)'])

mu_new = [dN.μ;dltg.μ[1]]
sigma_new = Matrix{Float64}(undef,size(mu_new)[1],size(mu_new)[1])
rho_new = Matrix{Float64}(undef,size(mu_new)[1],size(mu_new)[1])
sigma_new[1:3,1:3] .= dN.Σ
rho_new[1:3,1:3] .= cov2cor(dN.Σ)
rho_new[1,4] = 0.04 # From C3711005, not shown otherwise
rho_new[4,1] = rho_new[1,4]
rho_new[2,4] = 0.14 # From C3711005, not shown otherwise
rho_new[4,2] = rho_new[2,4]
rho_new[3,4] = cov2cor(dltg.Σ)[1,2]
rho_new[4,3] = cov2cor(dltg.Σ)[1,2]
shdl = sqrt(sigma_new[1,1])
stc = sqrt(sigma_new[2,2])
stg = sqrt(sigma_new[3,3])
sltg = sqrt(dltg.Σ[1,1])


sigma_new[1,4] = sltg*shdl*rho_new[1,4]
sigma_new[4,1] = sltg*shdl*rho_new[4,1]

sigma_new[2,4] = sltg*stc*rho_new[2,4]
sigma_new[4,2] = sltg*stc*rho_new[4,2]
sigma_new[3,4] = sltg*stg*rho_new[3,4]
sigma_new[4,3] = sltg*stg*rho_new[4,3]
sigma_new[4,4] = sltg^2

if dims == [1;2;3]
    d = fit(MvNormal,[log.(df.LBDHDDSI)';log.(df.LBDTCSI)';log.(df.LBDTRSI)'])
elseif dims == [1;2]
    d = fit(MvNormal,[log.(df.LBDHDDSI)';log.(df.LBDTCSI)'])
elseif dims == [1;3]
    d = fit(MvNormal,[log.(df.LBDHDDSI)';log.(df.LBDTRSI)'])
elseif dims == [2;3]
    d = fit(MvNormal,[log.(df.LBDTCSI)';log.(df.LBDTRSI)'])
elseif dims == [3;4]
    d = MvNormal(mu_new[3:4],sigma_new[3:4,3:4])
else
    d = MvNormal(mu_new,sigma_new)
end


return d

end

