#=
Specialized implementation of a M-H algorithm for generating plausible patients.

Based on: Rieger et al. 2018. (Matlab implementation)

Note: we put *some of* the pieces in place for multi-threading below, but right now this algorithm only runs on a single thread.
Similar to StaticArrays and the ODE solver, there is likely a performance boost possible, but unneeded at the moment.

=#
using Base.Threads, ProximalOperators, LinearAlgebra

function mh(n_pp::Integer,
    p::Vector{Float64},
    Sim_and_Score::Function, 
    DrawP::Function,
    OdeF::ODEProblem)::Matrix{Float64}
    
    nth = 1 # nthreads()
    n_param = size(p,1) # total number of parameters

    BURNIN = 0 # number of PPs before we compute the final RHO matrix

    global p_pp = Matrix{Float64}(undef,n_param,1)
    p0 = copy(p) # Initial parameter vector


   for t in 1:nth
        mh_count = 0 # Loop counter
        pp_count = 0 # Number of PP counter
        q0 = 0.000 # Initial proposal score
        q1 = 0.000 # Next proposal score
        p1 = Vector{Float64}(undef,n_pp) # Next parameter vector1
        p2 = Vector{Float64}(undef,n_pp) # Next parameter vector2
        #SigmaMat = copy(SigmaMat0)
        #RhoMat = copy(RhoMat0)
        burn_count = 0

        if (mod(n_pp,nth) == 0)
            n_pp_th = n_pp/nth
        else
            if (t == nth)
                n_pp_th = n_pp - (nth-1)*floor(n_pp/nth)
            else
                n_pp_th = floor(n_pp/nth)
            end
        end
        n_pp_th = Int64(n_pp_th)

        p_pp_th = Matrix{Float64}(undef,n_param,n_pp_th)

        p_pp_BURNIN = Matrix{Float64}(undef,n_param,BURNIN)

        while (pp_count < n_pp_th)
            mh_count += 1
        
            # First time through we need to run two guesses:
            if mh_count == 1

                # Solve initial parameter proposal:
                p1 = DrawP(p0)
                q1 = Sim_and_Score(p1,OdeF)
            end

            if burn_count == BURNIN
                burn_count += 1 # prevent recalculation every other loop
            end

            p2 = DrawP(p0)
            q2 = Sim_and_Score(p2,OdeF)

            # Core of the M-H algorithm, probalistic acceptance of the next guess:
            # We add onto q0 in the denominator to avoid an odd corner case where we were accepting a lot of
            # low probability PPs, or in other words q1 must be "meaningfully" better than q0 or we are not interested.
            if (q2/(q1 + 10.0*eps()) > rand())
                # Accept the next solution, record a PP in various ways:
                if burn_count < BURNIN
                    burn_count += 1
                    p_pp_BURNIN[:,burn_count] = p2
                else
                    pp_count += 1 # iterate the number of PPs
                    p_pp_th[:,pp_count] = p2
                end
                q1 = copy(q2) # Reset for the next iteration
                p1 = copy(p2) # Reset for the next iteration
            end # if (q2/q1 > rand())
        end # while (pp_count < n_pp_th)

        global p_pp = hcat(p_pp,p_pp_th) # Save our PPs, this extra step seemed necessary if we ever wanted multi-threading.
    end
    p_pp = p_pp[:,2:end] # remove the first dummy row

    return p_pp
end