module Synthesis_wave
export synthesize_noise_field

include("./gfandsource.jl")
import .GFandSource
using Dierckx, Random, Distributions, FFTW

function synthesize_noise_field!(;receiverdata, omega, t_trace, source_activatetime, Receiver_ID, config, rho, velocity, attenuation, sourceloc, receiverloc, scatterloc, fMrand)


    Rx1 = [receiverloc.x[Receiver_ID], receiverloc.y[Receiver_ID]]

    #dispersion interpolation
    spl_v = Spline1D(velocity.period, velocity.v; k=2, bc="nearest")
    #attenuation interpolation
    spl_alpha = Spline1D(attenuation.period, attenuation.alpha; k=2, bc="nearest")

    gf1 = zeros(Complex{Float64}, sourceloc.numofsource, 2*(length(omega)-1))

    #synthesize main phase
    for j = 1:sourceloc.numofsource

        Sx1 =  [sourceloc.x[j], sourceloc.y[j]]

        for k = 1:length(omega)-1
            period_k = 2*pi ./ omega[k]
            gf1[j, k] = GFandSource.gf(Rx1, Sx1, omega[k], spl_v(period_k), rho, spl_alpha(period_k))
            #gf1[j, k] = GFandSource.gf(Rx1, Sx1, omega[k], 3000.0, rho, 0.0)
        end

        #remove components at omega = zero
        gf1[j, 1] = 0.0

    end

    #synthesize scatter waves

    for j = 1:sourceloc.numofsource

        Sx1 =  [sourceloc.x[j], sourceloc.y[j]]

        for k = 1:length(omega)-1
            
            
            gf_sc_1 = zeros(Complex{Float64}, scatterloc.numofscatter)
            period_k = 2*pi ./ omega[k]

            for l = 1:scatterloc.numofscatter

                Scatter_x =  [scatterloc.x[l], scatterloc.y[l]]
                gf_sc_1[l] = GFandSource.gf_scatter(Rx1, Sx1, omega[k], Scatter_x, spl_v(period_k), 0.5 * scatterloc.strength[l], rho, spl_alpha(period_k))
            end
            gf1[j,k] += sum(gf_sc_1)
        end
        #avoid NaN at omega = 0
        gf1[j,1] = 0.0
    end
    

    #convolved with Ricker
    u_temp = zeros(Float64, length(t_trace));
    gftemp = zeros(Complex{Float64}, 2*(length(omega)-1))

    #stack each activation
    for j = 1:sourceloc.numofsource
        for k = 1:config.Source_average_num_per_unithour

            #fMrand = rand(Normal(config.Source_peakfreq_mean, config.Source_peakfreq_variance), sourceloc.numofsource)
            if config.Sourcetype == "ricker"

                fM = fMrand[(j-1)*sourceloc.numofsource + k]

                for l = 1:length(omega)-1
                    gftemp[l] = gf1[j,l] * GFandSource.rickerWavelet(omega[l], 2*pi*fM)
                end

            else
                error("Source type should be 'ricker' for the moment.")
            end

            if  config.IsGaussianNoise

                """
                NOT IMPLEMENTED!!

                Random.seed!(noiseseed[noiseid])
                noisearray = 1 .- (rand(-100:100, NumofSource, length(omega)-1)./100) * noiselevel

                for i = 1:sourceloc.numofsource
                    for j = 1:length(omega)-1
                        gf1[i,j] = gf1[i,j] * noisearray[i,j]

                    end
                end
                """
            end

            """
            Exercise 0: Plot displacement in time domain
            """

            # take inverse fft
            #u1test = ifft(sum(gf1, dims=1))
            #u1test = ifft(sum(gf2, dims=1))
            
            u1test = ifft(gftemp)

            """
            trace1 =  scatter(;x=t, y=real(u1test[1, 1:length(omega)]), mode="lines", name="Source No.")
            layout = Layout(width=800,height=600, xaxis_range = [0, maximum(t)])
            p = plot(trace1, layout) 

            for l = 2:sourceloc.numofsource
                addtraces!(p, scatter(;x=t, y=real(u1test[l, 1:length(omega)]), mode="lines", name="Source No."))
            end
            """

            #find activate time id
            tid = findfirst(x -> x>=source_activatetime[j, k], t_trace);

            u_temp[tid:tid+length(omega)-1] += real(u1test[1:length(omega)]);
        end
    end

    #store in struct
    receiverdata.Receiver_ID = Receiver_ID
    receiverdata.loc_x = receiverloc.x[Receiver_ID]
    receiverdata.loc_y = receiverloc.y[Receiver_ID]
    receiverdata.t = t_trace
    receiverdata.nt = length(t_trace)
    receiverdata.u =  u_temp

end

end
