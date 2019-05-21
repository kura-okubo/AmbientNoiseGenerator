"""
Sensitivity kernel for C1, C2 and C3
March 07, 2019
Kurama Okubo
"""

using Plots, FFTW, LinearAlgebra
gr()

#------------------------------------------------#
L = 6e4 # domain size x
W = 3e4 # domain size y

r_dist = 3e4
#make grid
gridsize = 1e3 #[m]

vel = 3e3 # wave speed between sensors [m/s]
rho = 3000.0 # density of medium [kg/m3]
Q = 1e9	# Attenuation term (quality factor) Q
dt = 0.05
T = 40

ifScatter = false #if scatter is added
#---parameters for scatters---# 
NumofScatter = 100;
scatternoiseseed = 1120;
scatterstrength = 5.0e6*ones(NumofScatter);

ifRicker = true # If Ricker wavelet is convolved
ifcheckRicker = false #debug Ricker wavelet
fM = 0.5 # [Hz] peak (center) frequency

noiseid = 1

Figdirroot = "/Users/kurama/Documents/kurama/research/sensitivity_kernel/fig/"
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

#MODE C1

#location of receiver 1 and 2
Rx1 = [-(r_dist/2), 0]
Rx2 = [(r_dist/2), 0]

#uniform source distribution

plotdomain = maximum([L/2, W/2])
p_all = plot(xlabel = "x [km]", 
    ylabel = "y [km]",
    title = "Geometry of source and receiver",
    xlim = (-plotdomain, plotdomain)./1e3,
    ylim = (-plotdomain, plotdomain)./1e3,
    )


"""
Plot source and receiver and scatter
"""

#plot source
"""
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
"""

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

"""
p_all = scatter!(Rx3[1,:]./1e3, Rx3[2,:]./1e3,
        markershape = :dtriangle,
        markersize = 10,
        markercolor = :green,
        label = "receiver 3")
"""

plot(p_all, size = (800, 800), legend=:topright)


savefig(Figdirroot*"geom_source_and_receiver."*Filetype) 


#location of receiver 3
#Rx3 = [-2.0*lx, lx]

#location of source
"""
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

NT = round(Int64, T/dt)
NumofFreq = NT
fs = 1/dt

#cc1 array for SN analysis
cc1_12_SNanalysis = zeros(Float64, 1, 2*NT-1)


"""
Synthesize signals
from v06. SNratio analysis
"""

ngrid_x = round(Int, L/gridsize)
ngrid_y = round(Int, W/gridsize)
Total_ngrid = (ngrid_x+1)*(ngrid_y+1)

xgrid = collect(range(-L, stop=L, length=ngrid_x+1)) 
ygrid = collect(range(-W, stop=W, length=ngrid_y+1))

global NT

#gf3 = zeros(Complex{Float64}, NumofSource, 2*(length(omega)-1))

#location of source

Sx = zeros(Float64, Total_ngrid , 2)

icount = 0
for i = 1:ngrid_x+1
    for j = 1:ngrid_y+1
        
        global icount += 1
        Sx[icount,:] = [xgrid[i], ygrid[j]]

        #uniform source distribution

        #if (xgrid[i] == L/2 && ygrid[j] == 0) || (xgrid[i] == -L/2 && ygrid[j] == 0)
        #    global icount += 1
        #    Sx[icount,:] = [xgrid[i], ygrid[j]]
        #    println(Sx[icount, :])  
        #end
        #if (xgrid[i] == L/2 && ygrid[j] == 0) 
        #    global icount += 1
        #    Sx[icount,:] = [xgrid[i], ygrid[j]]
        #    println(Sx[icount, :])  
        #end

    end
end

NumofSource = icount

#omega
omega = 2*pi.*(fs*(0:(NT/2))/NT);

gf1 = zeros(Complex{Float64}, NumofSource, 2*(length(omega)-1))
gf2 = zeros(Complex{Float64}, NumofSource, 2*(length(omega)-1))


for i = 1:NumofSource
    for j = 1:length(omega)-1

        gf1[i,j] = gf(Rx1, Sx[i,:], omega[j], vel, rho, Q)
        gf2[i,j] = gf(Rx2, Sx[i,:], omega[j], vel, rho, Q)

        if isinf(gf1[i,j]) 
            gf1[i,j] = 0 + 0im
        end

        if isinf(gf2[i,j]) 
            gf2[i,j] = 0 + 0im
        end

        #gf3[i,j] = gf(Rx3, Sx[i,:], omega[j], vel, rho, Q)
        
        #remove components at omega = zero

        gf1[i,1] = 0.0
        gf2[i,1] = 0.0

        #gf3[i,1] = 0.0

    end
end

if ifScatter
    for i = 1:NumofSource
        for j = 1:length(omega)-1

            gf_sc_1 = zeros(Complex{Float64}, NumofScatter)
            gf_sc_2 = zeros(Complex{Float64}, NumofScatter)
            #gf_sc_3 = zeros(Complex{Float64}, NumofScatter)

            for k = 1:NumofScatter
                gf_sc_1[k] = gf_scatter(Rx1, Sx[i,:], omega[j], Scatter_x[k,:], vel, scatterstrength[k], rho, Q)    
                gf_sc_2[k] = gf_scatter(Rx2, Sx[i,:], omega[j], Scatter_x[k,:], vel, scatterstrength[k], rho, Q)    
                #gf_sc_3[k] = gf_scatter(Rx3, Sx[i,:], omega[j], Scatter_x[k,:], vel, scatterstrength[k], rho, Q)    
            end

            #println("before add:$(gf1[i,j])\n")
            gf1[i,j] = gf1[i,j] + sum(gf_sc_1)
            gf2[i,j] = gf2[i,j] + sum(gf_sc_2)
            #gf3[i,j] = gf3[i,j] + sum(gf_sc_3)
            #println("after add:$(gf1[i,j])\n")

        end

        #avoid NaN at omega = 0
        gf1[i,1] = 0.0
        gf2[i,1] = 0.0
        #gf3[i,1] = 0.0
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
            #gf3[i,j] = gf3[i,j] * rickerWavelet(omega[j], 2*pi*fM)

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
    savefig(Figdirroot*"RickerWacelet."*Filetype) 

end


"""
Exercise 0: Plot displacement in time domain
"""

# take inverse fft
u1test = ifft(sum(gf1, dims=1))
u2test = ifft(sum(gf2, dims=1))
#u1test = ifft(sum(gf2, dims=1))
t = dt .* range(0, stop=length(omega)-1)

p1 = plot(t, real.(u1test[1:length(omega)]), line=(:black, 1, :solid),
#marker = (:cross, 2, :green),
xlabel = "Time [s]", 
ylabel = "u1 [m]",
title = "u1",
label="u1",
#xlim = (0, T/2),
#ylim = (1.2*minimum(real.(u1test[1:length(omega)])), maximum(real.(u1test[1:length(omega)])) * 1.2),
size=(600,300)
)


p1 = plot!(t, real.(u2test[1:length(omega)]), line=(:red, 1, :dash),
#marker = (:cross, 2, :green),
label="u2"
#xlim = (0, T/2),
#ylim = (1.2*minimum(real.(u1test[1:length(omega)])), maximum(real.(u1test[1:length(omega)])) * 1.2),
)

plot(p1, layout = (1,1), size = (600, 300), legend=true)
savefig(p1, Figdirroot*"u1."*Filetype) 


"""
Exercise 1: First order cross correlation C1:1->2
"""

signal_magnification = 1e6
minplotT = -T/2
maxplotT = T/2

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

# plot along azimuth

t_cc = dt .* range(-(NT-1), stop=NT-1)

cc1_12_stack = zeros(Complex{Float64}, 2*NT-1)

#estimate amplification coefficient
cc1test = vcat(ifft(cc1_12_neg[1, :])[end:-1:2], ifft(cc1_12_pos[1,:])) 
signal_magnification = 5*1/maximum(real.(cc1test))


for i = 1:NumofSource
    cc1test = vcat(ifft(cc1_12_neg[i, 1:end])[end:-1:2], ifft(cc1_12_pos[i,:]))
    global cc1_12_stack += cc1test
end

#plot stack
yplotrange = maximum(abs.(signal_magnification.*real.(cc1_12_stack)))

p_stack = plot(t_cc, signal_magnification.*real.(cc1_12_stack), 
    xlabel = "Lag time [s]", 
    ylabel = "Coherency",
    title = "Stacked C1:1-2",
    xlim = (minplotT, maxplotT),
    ylim = (-1.5*yplotrange, 1.5*yplotrange)
    #yticks = []
    )

plot(p_stack, size = (800, 600), legend=false)

savefig(Figdirroot*"cc1_1-2."*Filetype) 

global cc1_12_stack
#add for SN analysis
cc1_12_SNanalysis[noiseid, :] = real.(cc1_12_stack)


