"""
High order cross correlation (C2, C3) exercise

March 07, 2019
Kurama Okubo
"""

using Plots, FFTW, LinearAlgebra, Random, Statistics, DSP, Printf

#------------------------------------------------#
lx = 1e4	# distance of sensors R1 = (-lx, 0) R2 = (+lc, 0) [m]
R = 3e4	# radius of source location (S = (Rcos(theta), Rsin(theta))) [Hz]
vel= 3e3 # wave speed between sensors [m/s]
rho = 3000.0	# density of medium [kg/m3]
Q = 1e9	# Attenuation term (quality factor) Q
dt = 0.01
T = 100
NumofSource = 100;

ifScatter = false #if scatter is added
#---parameters for scatters---# 
NumofScatter = 100;
scatternoiseseed = 1120;
scatterstrength = 5.0e6*ones(NumofScatter);

#-----------------------------#


ifRicker = true# If Ricker wavelet is convolved
ifcheckRicker = false #debug Ricker wavelet
ifsourcenoise = false #add noise to source
noiselevel = 5 #[%] perturbation in frequency domain
noiseseed=1160:10:1160
fM = 0.5 # [Hz] peak (center) frequency

Figdirroot = "/Users/kurama/Documents/kurama/research/Highorder_crosscorrelation/code/fig_withsourcenoise_modified_1/"
Filetype = "pdf"
#------------------------------------------------#

if ~isdir(Figdirroot)
    mkdir(Figdirroot)
end

#The far field Green function (Fichtner, 2015 (doi:10.1093/gji/ggv182) eq. 12)

"""
gf(x::Array{Float64,1},xi::Array{Float64,1},omega::Float64, vel::Float64, rho::Float64, Q::Float64)
returns the Green’s function solution for two-dimensional space, under the far field approximation
"""
function gf(x::Array{Float64,1},xi::Array{Float64,1},omega::Float64, vel::Float64, rho::Float64, Q::Float64)
    r = norm(x.-xi)
    g1 = (1/sqrt((8*pi*omega*(1/vel) * r))) * exp(-im * (omega* (1/vel) * r + pi/4)) * exp(-(omega* r)/(2*vel*Q)) 
    return g1
end

"""
gf_scatters(x::Array{Float64,1},xi::Array{Float64,1},omega::Float64, vel::Float64, rho::Float64, Q::Float64)
returns the Green’s function solution for two-dimensional space, under the far field approximation
"""
function gf_scatter(x::Array{Float64,1}, xi::Array{Float64,1}, omega::Float64, xscatter::Array{Float64,1}, vel::Float64, scatterstrength::Float64, rho::Float64, Q::Float64) 

    g0_sc_x = gf(x, xscatter, omega, vel, rho, Q)
    g0_xi_sc = gf(xscatter, xi, omega, vel, rho, Q)
    g1_scatter = g0_sc_x * (omega^2/vel^2) * scatterstrength * g0_xi_sc
    #g1_scatter = g0_sc_x * scatterstrength * g0_xi_sc
    #println((omega^2/vel^2))
    return g1_scatter
end



"""
rickerWavelet(omega)
returns zero-phase Ricker wavelet
"""
function rickerWavelet(omega::Float64, omega0::Float64)
    ricker = ((2*omega^2)/(sqrt(pi) * omega0^3)) * exp(-(omega^2 / omega0^2))
    return ricker
end

#location of receiver 1 and 2
Rx1 = [-lx, 0]
Rx2 = [lx, 0]

#location of receiver 3
Rx3 = [-2.0*lx, lx]

#location of source
theta = zeros(NumofSource, 1)
Sx = zeros(NumofSource, 2)
Scatter_x = zeros(NumofScatter, 2)
for i = 1:NumofSource
    theta[i] = 2*pi*(i-1)/NumofSource
    Sx[i,:] = R .* [cos(theta[i]), sin(theta[i])]  
end

#location of scatter
Random.seed!(scatternoiseseed)
noisearray = 2 .* (rand(-100:100, NumofScatter, 2)./100)
for i = 1:NumofScatter
    Scatter_x[i,:] = noisearray[i,:] .* [lx, lx]
    #Scatter_x[i,:] = [2.0*lx, -0.9*lx]
end


"""
Plot source and receiver and scatter
"""

plotdomain = maximum(abs.(Rx3))
p_all = plot(xlabel = "x [km]", 
    ylabel = "y [km]",
    title = "Geometry of source and receiver",
    xlim = (-2.0*plotdomain, 2.0*plotdomain)./1e3,
    ylim = (-2.0*plotdomain, 2.0*plotdomain)./1e3,
    )

#plot source
p_all = scatter!(Sx[:,1]./1e3, Sx[:,2]./1e3,
        markershape = :star6,
        markersize = 10,
        markercolor = :black,
        label = "source")

p_all = scatter!(Scatter_x[:,1]./1e3, Scatter_x[:,2]./1e3,
        markershape = :circle,
        markersize = 6,
        markercolor = :black,
        label = "scatter")

#plot receiver
p_all = scatter!(Rx1[1,:]./1e3, Rx1[2,:]./1e3,
        markershape = :dtriangle,
        markersize = 10,
        markercolor = :red,
        label = "receiver 1")
p_all = scatter!(Rx2[1,:]./1e3, Rx2[2,:]./1e3,
        markershape = :dtriangle,
        markersize = 10,
        markercolor = :blue,
        label = "receiver 2")
p_all = scatter!(Rx3[1,:]./1e3, Rx3[2,:]./1e3,
        markershape = :dtriangle,
        markersize = 10,
        markercolor = :green,
        label = "receiver 3")


plot(p_all, size = (800, 800), legend=:topright)


savefig(Figdirroot*"geom_source_and_receiver."*Filetype) 


"""
Synthesize signals
from v06. SNratio analysis
"""

NT = round(Int64, T/dt)
NumofFreq = NT
fs = 1/dt

#cc1 and cc2 array for SN analysis
cc1_12_SNanalysis = zeros(Float64, length(noiseseed), 2*NT-1)
cc2_12_SNanalysis = zeros(Float64, length(noiseseed), 2*NT-1)

global NT

#for noiseid = 1:length(noiseseed)
    noiseid = 1


    Figdir = Figdirroot*"noiseid_$noiseid/"

    if ~isdir(Figdir)
        mkdir(Figdir)
    end
    #omega
    omega = 2*pi.*(fs*(0:(NT/2))/NT);

    gf1 = zeros(Complex{Float64}, NumofSource, 2*(length(omega)-1))
    gf2 = zeros(Complex{Float64}, NumofSource, 2*(length(omega)-1))
    gf3 = zeros(Complex{Float64}, NumofSource, 2*(length(omega)-1))

    for i = 1:NumofSource
        for j = 1:length(omega)-1
            gf1[i,j] = gf(Rx1, Sx[i,:], omega[j], vel, rho, Q)
            gf2[i,j] = gf(Rx2, Sx[i,:], omega[j], vel, rho, Q)
            gf3[i,j] = gf(Rx3, Sx[i,:], omega[j], vel, rho, Q)
            
            #remove components at omega = zero

            gf1[i,1] = 0.0
            gf2[i,1] = 0.0
            gf3[i,1] = 0.0

        end
    end

    if ifScatter
        for i = 1:NumofSource
            for j = 1:length(omega)-1

                gf_sc_1 = zeros(Complex{Float64}, NumofScatter)
                gf_sc_2 = zeros(Complex{Float64}, NumofScatter)
                gf_sc_3 = zeros(Complex{Float64}, NumofScatter)

                for k = 1:NumofScatter
                    gf_sc_1[k] = gf_scatter(Rx1, Sx[i,:], omega[j], Scatter_x[k,:], vel, scatterstrength[k], rho, Q)    
                    gf_sc_2[k] = gf_scatter(Rx2, Sx[i,:], omega[j], Scatter_x[k,:], vel, scatterstrength[k], rho, Q)    
                    gf_sc_3[k] = gf_scatter(Rx3, Sx[i,:], omega[j], Scatter_x[k,:], vel, scatterstrength[k], rho, Q)    
                end

                #println("before add:$(gf1[i,j])\n")
                gf1[i,j] = gf1[i,j] + sum(gf_sc_1)
                gf2[i,j] = gf2[i,j] + sum(gf_sc_2)
                gf3[i,j] = gf3[i,j] + sum(gf_sc_3)
                #println("after add:$(gf1[i,j])\n")

            end

            #avoid NaN at omega = 0
            gf1[i,1] = 0.0
            gf2[i,1] = 0.0
            gf3[i,1] = 0.0
        end
    end

    #rptest = plot(real(gf1[1,:]), line=(:black, 1, :solid))
    #plot(rptest, size = (800, 800), legend=false)
    #sleep(100)

    if ifRicker
        for i = 1:NumofSource
            for j = 1:length(omega)-1
                gf1[i,j] = gf1[i,j] * rickerWavelet(omega[j], 2*pi*fM)
                gf2[i,j] = gf2[i,j] * rickerWavelet(omega[j], 2*pi*fM)
                gf3[i,j] = gf3[i,j] * rickerWavelet(omega[j], 2*pi*fM)

            end
        end
    end

    if ifsourcenoise

        Random.seed!(noiseseed[noiseid])
        noisearray = 1 .- (rand(-100:100, NumofSource, length(omega)-1)./100) * noiselevel

        for i = 1:NumofSource
            for j = 1:length(omega)-1
                gf1[i,j] = gf1[i,j] * noisearray[i,j]
                gf2[i,j] = gf2[i,j] * noisearray[i,j]
                gf3[i,j] = gf3[i,j] * noisearray[i,j]

            end
        end
    end

    #make it symmetric
    """
    for i = 1:1
        for j = 1:length(omega)-1
            gf1[i, Int((length(omega)-1) + j)] = real(gf1[i, Int(length(omega) - j)]) + im*-imag(gf1[i, Int(length(omega) - j)]);
            gf2[i, Int((length(omega)-1) + j)] = real(gf2[i, Int(length(omega) - j)]) + im*-imag(gf2[i, Int(length(omega) - j)]);
            gf3[i, Int((length(omega)-1) + j)] = real(gf3[i, Int(length(omega) - j)]) + im*-imag(gf3[i, Int(length(omega) - j)]);
        end
    end
    """

    if ifcheckRicker
        """
        check Ricker wavelet
        """
        rwave = zeros(Complex{Float64}, 2*(length(omega)-1))
        for j = 1:length(omega)
            rwave[j] = rickerWavelet(omega[j], 2*pi*fM);
        end
        fs_ricker = omega ./ (2*pi)
        rp1 = plot(fs_ricker, real(rwave)[1:length(omega)], line=(:black, 1, :solid),
            xlabel = "Frequency[Hz]", 
            ylabel = "Spectral density",
            title = "Ricker wavelet",
            xlim = (0.01*fM, 5*fM),
            ylim = (minimum(real(rwave)[1:length(omega)])+1e-3, 1.2*maximum(real(rwave)[1:length(omega)])),
            #xaxis=:log,
            #yaxis=:log
           )

        URicker = fftshift(ifft(rwave))
        t_ricker = dt*range(-(length(omega)-1), stop=(length(omega)-1))
        TD = sqrt(6)/(pi*fM)
        TR = TD/sqrt(3)
        rp2 = plot(t_ricker[1:end-1], real(URicker), line=(:black, 1, :solid),
            xlabel = "Time [s]", 
            ylabel = "u1 [m]",
            xlim = (-3*TD, 3*TD),
            )
        rp2 = plot!([TD, TD]./2, [-1, 1] .* 1e-3, line=(:red, 1, :dash))
        rp2 = plot!([TR, TR]./2, [-1, 1] .* 1e-3, line=(:blue, 1, :dash))
        plot(rp1, rp2, layout = (2,1), size = (800, 800), legend=false)
        savefig(Figdir*"RickerWacelet."*Filetype) 

    end

    """
    Exercise 0: Plot displacement in time domain
    """

    # take inverse fft
    u1test = ifft(sum(gf1, dims=1))
    t = dt .* range(0, stop=length(omega)-1)

    p1 = plot(t, real.(u1test[1:length(omega)]), line=(:black, 1, :solid),
        #marker = (:cross, 2, :green),
        xlabel = "Time [s]", 
        ylabel = "u1 [m]",
        title = "u1",
        #xlim = (0, T/2),
        #ylim = (1.2*minimum(real.(u1test[1:length(omega)])), maximum(real.(u1test[1:length(omega)])) * 1.2),
        size=(600,300)
        )

    plot(p1, layout = (1,1), size = (600, 300), legend=false)
    savefig(p1, Figdir*"u1."*Filetype) 

    # plot along azimuth

    signal_magnification = 5*1/maximum(real.(u1test))
    xticks = 0:45:360

    p_all = plot(xlabel = "azimuth [deg]", 
        ylabel = "Time [s]",
        title = "Displacement at R1",
        xlim = (-20, 380),
        ylim = (0, T/2),
        xticks = xticks
        )

    for i = 1:NumofSource
        u1test = ifft(gf1[i,:])
        u2test = ifft(gf2[i,:])
        azimuth =  rad2deg(theta[i])  
        p_all = plot!(azimuth.+ signal_magnification.*real.(u1test[1:length(omega)]), t, line=(:red, 1, :solid))
        p_all = plot!(azimuth.+ signal_magnification.*real.(u2test[1:length(omega)]), t, line=(:blue, 1, :solid))
    end

    #plot at 360
    u1last = ifft(gf1[1,:])
    u2last = ifft(gf2[1,:])
    azimuthlast = 360
    p_all = plot!(azimuthlast.+ signal_magnification.*real.(u1last[1:length(omega)]), t, line=(:red, 1, :solid))
    p_all = plot!(azimuthlast.+ signal_magnification.*real.(u2last[1:length(omega)]), t, line=(:blue, 1, :solid))


    plot(p_all, layout = (1,1), size = (600, 600), legend=false)
    savefig(p_all, Figdir*"u1allazimuth."*Filetype) 

    """
    Exercise 1: First order cross correlation C1:1->2
    """

    signal_magnification = 1e6
    maxplotT = 12
    xticks = 0:45:360

    cc1_12_pos = zeros(Complex{Float64}, NumofSource, NT)
    cc1_12_neg = zeros(Complex{Float64}, NumofSource, NT)

    sum_cc1_12_pos = zeros(Complex{Float64}, NT)
    sum_cc1_12_neg = zeros(Complex{Float64}, NT)

    ccid = -(NT-1):NT-1

    #modified summation
    for i = 1:NumofSource
        for j = 1:NT
            
            cc1_12_pos[i,j] = conj(gf1[i,j]) * gf2[i,j]
            cc1_12_neg[i,j] = gf1[i,j] * conj(gf2[i,j])

        end
    end

	sum_cc1_12_pos = conj(sum(gf1,dims=1)).*sum(gf2,dims=1)
	sum_cc1_12_neg = sum(gf1,dims=1).*conj(sum(gf2,dims=1))


    # plot along azimuth

    minplotT = -20
    maxplotT = 20
    xticks = 0:45:360

    p_all = plot(xlabel = "azimuth [deg]", 
        ylabel = "Lag time [s]",
        title = "C1 Coherency",
        xlim = (-20, 380),
        ylim = (minplotT, maxplotT),
        xticks = xticks
        )

    t_cc = dt .* range(-(NT-1), stop=NT-1)

    cc1_12_stack = zeros(Complex{Float64}, 2*NT-1)

    #estimate amplification coefficient
    cc1test = vcat(ifft(cc1_12_neg[1, 2:end])[end:-1:1], ifft(cc1_12_pos[1,:])) 
    signal_magnification = 5*1/maximum(real.(cc1test))

    for i = 1:NumofSource
    	cc1test = vcat(ifft(cc1_12_neg[i, 2:end])[end:-1:1], ifft(cc1_12_pos[i,:])) 
        azimuth =  rad2deg(theta[i])
        p_all = plot!(azimuth.+ signal_magnification.*real.(cc1test[1:end]), t_cc[1:end], line=(:black, 1, :solid))
        #global cc1_12_stack += cc1test
    end

    #cc1stack is WRONG. modified summation
    cc1_12_stack  = vcat(ifft(sum_cc1_12_neg[2:end])[end:-1:1], ifft(sum_cc1_12_pos[:])) 

    #plot at 360
    cc1test = vcat(ifft(cc1_12_neg[1, 2:end])[end:-1:1], ifft(cc1_12_pos[1,:]) ) 
    azimuthlast = 360
    p_all = plot!(azimuthlast.+ signal_magnification.*real.(cc1test[1:end]), t_cc[1:end], line=(:black, 1, :solid))

    #plot stack
    xplotrange = maximum(abs.(signal_magnification.*real.(cc1_12_stack)))

    p_stack = plot(signal_magnification.*real.(cc1_12_stack), t_cc, xlabel = "", 
        ylabel = "Lag time [s]",
        title = "Stacked C1:1-2",
        xlim = (-1.5*xplotrange, 1.5*xplotrange),
        ylim = (minplotT, maxplotT),
        xticks = []
        )

    layout = @layout [a{0.8w} b{0.2w}]
    plot(p_all, p_stack, layout = layout, yflip = true, size = (800, 600), legend=false)
    savefig(Figdir*"cc1_1-2."*Filetype) 
    
    global cc1_12_stack
    #add for SN analysis
    cc1_12_SNanalysis[noiseid, :] = real.(cc1_12_stack)

    """
    Exercise 2: Second order cross correlation C2:1->2
    """
    #CC1: 3->1, 3-2
    cc1_31_pos = zeros(Complex{Float64}, NumofSource, NT)
    cc1_31_neg = zeros(Complex{Float64}, NumofSource, NT)
    cc1_32_pos = zeros(Complex{Float64}, NumofSource, NT)
    cc1_32_neg = zeros(Complex{Float64}, NumofSource, NT)

    for i = 1:NumofSource
        for j = 1:NT
            
            cc1_31_pos[i,j] = conj(gf3[i,j]) * gf1[i,j]
            cc1_31_neg[i,j] = gf3[i,j] * conj(gf1[i,j])

            cc1_32_pos[i,j] = conj(gf3[i,j]) * gf2[i,j]
            cc1_32_neg[i,j] = gf3[i,j] * conj(gf2[i,j])
        end
    end

    #check cc1_31 and cc1_32
    # plot along azimuth
    minplotT = -20
    maxplotT = 20
    xticks = 0:45:360

    p_all = plot(xlabel = "azimuth [deg]", 
        ylabel = "Lag time [s]",
        title = "C1:3-1 and 3-2 Coherency",
        xlim = (-20, 380),
        ylim = (minplotT, maxplotT),
        xticks = xticks
        )

    t_cc = dt .* range(-(NT-1), stop=NT-1)

    cc1_31_stack = zeros(Complex{Float64}, 2*NT-1)
    cc1_32_stack = zeros(Complex{Float64}, 2*NT-1)

    #estimate amplification coefficient
    cc1_31test = vcat(ifft(cc1_31_neg[1, 2:end])[end:-1:1], ifft(cc1_31_pos[1,:])) 
    signal_magnification = 5*1/maximum(real.(cc1test))

    for i = 1:NumofSource
        cc1_31test = vcat(ifft(cc1_31_neg[i, 2:end])[end:-1:1], ifft(cc1_31_pos[i,:])) 
        cc1_32test = vcat(ifft(cc1_32_neg[i, 2:end])[end:-1:1], ifft(cc1_32_pos[i,:])) 
        azimuth =  rad2deg(theta[i])
        p_all = plot!(azimuth.+ signal_magnification.*real.(cc1_31test[1:end]), t_cc[1:end], line=(:red, 1, :solid))
        p_all = plot!(azimuth.+ signal_magnification.*real.(cc1_32test[1:end]), t_cc[1:end], line=(:blue, 1, :solid))
        global cc1_31_stack += cc1_31test
        global cc1_32_stack += cc1_32test
        #println("$i:   ", maximum(abs.(cc1_31_stack)))
    end

    global cc1_31_stack
    global cc1_32_stack

    #plot at 360
    cc1_31test = vcat(ifft(cc1_31_neg[1, 2:end])[end:-1:1], ifft(cc1_31_pos[1,:]) ) 
    cc1_32test = vcat(ifft(cc1_32_neg[1, 2:end])[end:-1:1], ifft(cc1_32_pos[1,:]) ) 
    azimuthlast = 360
    p_all = plot!(azimuthlast.+ signal_magnification.*real.(cc1_31test[1:end]), t_cc[1:end], line=(:red, 1, :solid))
    p_all = plot!(azimuthlast.+ signal_magnification.*real.(cc1_32test[1:end]), t_cc[1:end], line=(:blue, 1, :solid))

    xplotrange = maximum(abs.(signal_magnification.*real.(cc1_31_stack)))

    p_stack = plot(signal_magnification.*real.(cc1_31_stack), t_cc, line=(:red, 1, :solid),
    	xlabel = "", 
        ylabel = "Lag time [s]",
        xlim = (-1.5*xplotrange, 1.5*xplotrange),
        ylim = (minplotT, maxplotT),
        xticks = []
        )
    p_stack = plot!(signal_magnification.*real.(cc1_32_stack), t_cc, line=(:blue, 1, :solid))


    layout = @layout [a{0.8w} b{0.2w}]
    plot(p_all, p_stack, layout = layout, yflip = true, size = (800, 600), legend=false)
    savefig(Figdir*"cc1_3-1and3-2."*Filetype) 

    """
    Plot positive side C2POS and negative side C2NEG
    """

    #CC2POS and CC2NEG
    cc2POS_12_pos = zeros(Complex{Float64}, NumofSource, NT)
    cc2POS_12_neg = zeros(Complex{Float64}, NumofSource, NT)

    cc2NEG_12_pos = zeros(Complex{Float64}, NumofSource, NT)
    cc2NEG_12_neg = zeros(Complex{Float64}, NumofSource, NT)

    for i = 1:NumofSource
        for j = 1:NT
            
            cc2POS_12_pos[i,j] = conj(cc1_31_pos[i,j]) * cc1_32_pos[i,j]
            cc2POS_12_neg[i,j] = cc1_31_pos[i,j] * conj(cc1_32_pos[i,j])

    		cc2NEG_12_pos[i,j] = conj(cc1_31_neg[i,j]) * cc1_32_neg[i,j]
            cc2NEG_12_neg[i,j] = cc1_31_neg[i,j] * conj(cc1_32_neg[i,j])
        end
    end


    # plot along azimuth
    minplotT = -20
    maxplotT = 20
    xticks = 0:45:360

    p_all = plot(xlabel = "azimuth [deg]", 
        ylabel = "Lag time [s]",
        title = "C2 Coherency",
        xlim = (-20, 380),
        ylim = (minplotT, maxplotT),
        xticks = xticks
        )

    t_cc = dt .* range(-(NT-1), stop=NT-1)

    cc2POS_12_stack = zeros(Complex{Float64}, 2*NT-1)
    cc2NEG_12_stack = zeros(Complex{Float64}, 2*NT-1)

    #estimate amplification coefficient
    cc2_POStest = vcat(ifft(cc2POS_12_neg[1, 2:end])[end:-1:1], ifft(cc2POS_12_pos[1,:])) 
    signal_magnification = 5*1/maximum(real.(cc2_POStest))

    for i = 1:NumofSource
    	cc2_POStest = vcat(ifft(cc2POS_12_neg[i, 2:end])[end:-1:1], ifft(cc2POS_12_pos[i,:])) 
    	cc2_NEGtest = vcat(ifft(cc2NEG_12_neg[i, 2:end])[end:-1:1], ifft(cc2NEG_12_pos[i,:])) 
        azimuth =  rad2deg(theta[i])
        p_all = plot!(azimuth.+ signal_magnification.*real.(cc2_POStest[1:end]), t_cc[1:end], line=(:red, 1, :solid))
        p_all = plot!(azimuth.+ signal_magnification.*real.(cc2_NEGtest[1:end]), t_cc[1:end], line=(:blue, 1, :solid))
        global cc2POS_12_stack += cc2_POStest
        global cc2NEG_12_stack += cc2_NEGtest
    end

    global cc2POS_12_stack
    global cc2NEG_12_stack

    #plot at 360
    cc2_POStest = vcat(ifft(cc2POS_12_neg[1, 2:end])[end:-1:1], ifft(cc2POS_12_pos[1,:]) ) 
    cc2_NEGtest = vcat(ifft(cc2NEG_12_neg[1, 2:end])[end:-1:1], ifft(cc2NEG_12_pos[1,:]) ) 
    azimuthlast = 360
    p_all = plot!(azimuthlast.+ signal_magnification.*real.(cc2_POStest[1:end]), t_cc[1:end], line=(:red, 1, :solid))
    p_all = plot!(azimuthlast.+ signal_magnification.*real.(cc2_NEGtest[1:end]), t_cc[1:end], line=(:blue, 1, :solid))

    xplotrange = maximum(abs.(signal_magnification.*real.(cc2POS_12_stack)))

    #plot stack
    p_stack = plot(signal_magnification.*real.(cc2POS_12_stack), t_cc, line=(:red, 1, :solid),
    	xlabel = "", 
        ylabel = "Lag time [s]",
        xlim = (-1.5*xplotrange, 1.5*xplotrange),
        ylim = (minplotT, maxplotT),
        xticks = []
        )

    p_stack = plot!(signal_magnification.*real.(cc2NEG_12_stack), t_cc, line=(:blue, 1, :solid))

    p_stack_POSandNEG = plot(signal_magnification.*(real.(cc2POS_12_stack) .+ real.(cc2NEG_12_stack)), t_cc, line=(:black, 1, :solid),
    	xlabel = "", 
        ylabel = "Lag time [s]",
        title = "Sum up",
        xlim = (-1.5*xplotrange, 1.5*xplotrange),
        ylim = (minplotT, maxplotT),
        xticks = []
        )

    layout = @layout [a{0.8w} b{0.1w} c{0.1w}]
    plot(p_all, p_stack, p_stack_POSandNEG, layout = layout, yflip = true, size = (800, 600), legend=false)
    savefig(Figdir*"cc3_1-2."*Filetype)

    #add for SN analysis
    global cc2_12_SNanalysis[noiseid, :] = real.(cc2POS_12_stack) .+ real.(cc2NEG_12_stack)
#end

#Signal-to-Noise ratio analysis (Clarke et al. 2011 eq (1)-(3))
N_noise = length(noiseseed)

c1_SNvariation = zeros(Float64, 2*NT-1)
c1_SNsignal = zeros(Float64, 2*NT-1)
c2_SNvariation = zeros(Float64, 2*NT-1)
c2_SNsignal = zeros(Float64, 2*NT-1)

cc1_var_temp = zeros(Float64, 2*NT-1)
cc2_var_temp = zeros(Float64, 2*NT-1)

c1_SNvariation = std(cc1_12_SNanalysis, dims=1)
c2_SNvariation = std(cc2_12_SNanalysis, dims=1)

c1_SNsignal = abs.(hilbert(mean(cc1_12_SNanalysis, dims=1)))
c2_SNsignal = abs.(hilbert(mean(cc2_12_SNanalysis, dims=1)))

#SNR is averaged over t = [0, end]
cc1_SNR_mean = mean(c1_SNsignal./c1_SNvariation)
cc2_SNR_mean = mean(c2_SNsignal./c2_SNvariation)

cc1_SNR_max = maximum(c1_SNsignal./c1_SNvariation)
cc2_SNR_max = maximum(c2_SNsignal./c2_SNvariation)

#plot cc and SNR
minplotT = -25
maxplotT = 25

#normalized amplitude with each cc; focus on the SNR
xplotrange_cc1 = maximum(abs.(cc1_12_SNanalysis))
xplotrange_cc2 = maximum(abs.(cc2_12_SNanalysis))

t_cc = dt .* range(-(NT-1), stop=NT-1)

p_SNR_cc1 = plot(mean(cc1_12_SNanalysis, dims=1)'./xplotrange_cc1, t_cc, line=(:black, 1, :solid),
    #xlabel = "CC1:1-2 SNR ="*@sprintf("avg. %4.2f max. %4.2f", cc1_SNR_mean, cc1_SNR_max), 
    xlabel = "Normalized Cross-correlation function", 
    ylabel = "Lag time [s]",
    xlim = (-1.1, 1.1),
    ylim = (minplotT, maxplotT),
    size = (400, 800),
    label="C1",
    xticks = []
    )

p_SNR_cc1 = plot!(mean(cc2_12_SNanalysis, dims=1)'./xplotrange_cc2, t_cc, line=(:red, 1, :solid),
    #xlabel = "CC2:1-2 mean(SNR) ="*@sprintf("avg. %4.2f max. %4.2f", cc2_SNR_mean, cc2_SNR_max), 
    ylabel = "Lag time [s]",
    #xlim = (-1.5*xplotrange_cc2, 1.5*xplotrange_cc2),
    ylim = (minplotT, maxplotT),
    size = (400, 800),
    label="C2"
    #xticks = []
    )

#layout = @layout [a{1.0w}]
plot(p_SNR_cc1, yflip = true, size = (400, 800), legend=true)
savefig(Figdirroot*"cc_SNR."*Filetype)


