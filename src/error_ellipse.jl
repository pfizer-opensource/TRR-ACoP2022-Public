function errorellipse(d::Distribution)
    # Returns a 2D error ellipse at the 95% prediction interval
    # Future exensions:
    #   Expansion to 3D+ (projections)
    #   Assumes distribution is MvLogNormal

    nvar = 100000 # number of points to estimate ellipse using
    nplot = nvar # number of points for plotting
    r = rand(d,nvar) # draw variates
    rcov = StatsBase.cov(r')
    reigval = LinearAlgebra.eigvals(rcov)
    reigvec = LinearAlgebra.eigvecs(rcov)

    reigval_max = findmax(reigval)[1]
    reigvec_max = reigvec[:,findmax(reigval)[2]]
    reigval_min = findmin(reigval)[1]
    reigvec_min = reigvec[:,findmin(reigval)[2]]
    rangle = atan(reigvec_max[2],reigvec_max[1])

    if (rangle < 0)
        rangle = rangle + 2*pi
    end

    rmean = mean(r,dims=2)

    chisquare_val = 2.4477
    theta_grid = LinRange(0,2*pi,nplot)
    phi = rangle
    X0 = rmean[1]
    Y0 = rmean[2]
    a = chisquare_val*sqrt(reigval_max)
    b = chisquare_val*sqrt(reigval_min)

    R = [cos(phi) sin(phi);
        -sin(phi) cos(phi)]

    ellipsexy = Matrix{Float64}(undef,nplot,2)
    ellipsexy[:,1] = a*cos.(theta_grid)
    ellipsexy[:,2] = b*sin.(theta_grid)

    ellipsexy = ellipsexy * R
    ellipsexy[:,1] = ellipsexy[:,1] .+ X0
    ellipsexy[:,2] = ellipsexy[:,2] .+ Y0

    return ellipsexy

end


# https://discourse.julialang.org/t/plot-ellipse-in-makie/82814/4
function getellipsepoints(cx, cy, rx, ry, θ)
	t = range(0, 2*pi, length=100)
	ellipse_x_r = @. rx * cos(t)
	ellipse_y_r = @. ry * sin(t)
	R = [cos(θ) sin(θ); -sin(θ) cos(θ)]
	r_ellipse = [ellipse_x_r ellipse_y_r] * R
	x = @. cx + r_ellipse[:,1]
	y = @. cy + r_ellipse[:,2]
	(x,y)
end

# https://discourse.julialang.org/t/plot-ellipse-in-makie/82814/4
function getellipsepoints(μ,Σ, confidence=0.95)
    
    # Modified, I would rather pass in a Distribution:
    μ = d.μ
    Σ = d.Σ


	quant = quantile(Chisq(2), confidence) |> sqrt
	cx = μ[1]
	cy =  μ[2]
	
	egvs = eigvals(Σ)
	if egvs[1] > egvs[2]
		idxmax = 1
		largestegv = egvs[1]
		smallesttegv = egvs[2]
	else
		idxmax = 2
		largestegv = egvs[2]
		smallesttegv = egvs[1]
	end

	rx = quant*sqrt(largestegv)
	ry = quant*sqrt(smallesttegv)
	
	eigvecmax = eigvecs(Σ)[:,idxmax]
	θ = atan(eigvecmax[2]/eigvecmax[1])
 	if θ < 0
		θ += 2*π
	end

	getellipsepoints(cx, cy, rx, ry, θ)
end