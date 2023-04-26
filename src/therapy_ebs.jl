using Distributions

function return_therapy_ebs()

nv = 1000

# DeGrooth, Dalcetrapib 600mg, 4wks:
obs_labels = ["HDL-Ch";"TPCh";"TPTG"]
ther_labels = ["Dalcetrapib_600mg_QD";"Evolocumab_140mg_Q2W"]

#df = DataFrame

delta_600mg = [0.32;0.0;-0.1] # [mM]
delta_0mg = [0.04;0.0;0.0] # [mM]
std_delta_600mg = [0.22;0.6;0.5] # [mM]
std_delta_0mg = [0.15;0.5;0.4] # [mM]

base_600mg = [1.21;5.7;1.7]
base_0mg = [1.16;5.6;1.5]
std_base_600mg = [0.25;1.0;0.9]
std_base_0mg = [0.23;1.1;0.7]

dcpb_N = 50
rho = 0.0

variates = Matrix{Float64}(undef,nv,size(obs_labels)[1])

 for ii = 1:size(delta_600mg)[1]
    # Create the delta distributions:

    # Create the base distributions:
    (mu,sigma) = lognstat(base_600mg[ii],std_base_600mg[ii]^2)
    CovMat = [std_delta_600mg[ii]^2 rho*std_delta_600mg[ii]*sigma;
                rho*std_delta_600mg[ii]*sigma   sigma^2]
    muMat = [delta_600mg[ii];mu]
    d_600mg = MvNormal(muMat,CovMat)
    
    (mu,sigma) = lognstat(base_0mg[ii],std_base_0mg[ii]^2)
    CovMat = [std_delta_0mg[ii]^2 rho*std_delta_0mg[ii]*sigma;
                rho*std_delta_0mg[ii]*sigma   sigma^2]
    muMat = [delta_0mg[ii];mu]
    d_0mg = MvNormal(muMat,CovMat)

    r600mg = rand(d_600mg,nv)
    r0mg = rand(d_0mg,nv)
    variates[:,ii] = @. 100.0*(r600mg[1,:]/exp(r600mg[2,:]) - r0mg[1,:]/exp(r0mg[2,:]))
    
 end

dcpb_eb_mean = mean(variates,dims=1)
dcpb_eb_CI_width = quantile(TDist(dcpb_N -1),1-0.05/2)*std(variates,dims=1)/sqrt(dcpb_N)

dcpb_eb_ci = [dcpb_eb_mean .- dcpb_eb_CI_width;
                dcpb_eb_mean .+ dcpb_eb_CI_width]


# Evolucamab
ev_N = (993+586)/2

# 95% CIs:
ev_mean = [6.90 -40.48 -17.35]
ev_95CI = [5.37 -45.33 -23.50;
            8.43 -35.62 -11.20]

ev_eb_CI = [(ev_mean .- ev_95CI[1,:]');
            (ev_95CI[2,:]' .- ev_mean)]



eb = [dcpb_eb_ci ev_95CI]#[dcpb_eb_CI ev_eb_CI]

return eb

end