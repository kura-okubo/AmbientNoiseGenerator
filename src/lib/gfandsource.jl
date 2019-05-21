module GFandSource
export gf, gf_scatte, rickerWavelet, whiten
using FFTW, LinearAlgebra, DSP

"""
gf(x::Array{Float64,1},xi::Array{Float64,1},omega::Float64, vel::Float64, rho::Float64, alpha::Float64)
returns the Green’s function solution for two-dimensional space, under the far field approximation
"""
function gf(x::Array{Float64,1},xi::Array{Float64,1},omega::Float64, vel::Float64, rho::Float64, alpha::Float64)
    r = norm(x.-xi)
    #g1 = (1/sqrt((8*pi*omega*(1/vel) * r))) * exp(-im * (omega* (1/vel) * r + pi/4)) * exp(-(omega* r)/(2*vel*Q)) 
    g1 = (1/sqrt((8*pi*omega*(1/vel) * r))) * exp(-im * (omega* (1/vel) * r + pi/4)) * exp(-alpha * r) 
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


"""
    whiten(A, freqmin, freqmax, fs, pad=100)
Whiten spectrum of time series `A` between frequencies `freqmin` and `freqmax`.
Uses real fft to speed up computation.
Returns the whitened (single-sided) fft of the time series.
# Arguments
- `A::AbstractArray`: Time series.
- `fs::Real`: Sampling rate of time series `A`.
- `freqmin::Real`: Pass band low corner frequency.
- `freqmax::Real`: Pass band high corner frequency.
- `pad::Int`: Number of tapering points outside whitening band.
"""
function whiten(A::AbstractArray, freqmin::Float64, freqmax::Float64, fs::Float64;
                pad::Int=50)
    if ndims(A) == 1
        A = reshape(A,size(A)...,1) # if 1D array, reshape to (length(A),1)
    end

    # get size and convert to Float32
    M,N = size(A)
    if eltype(A) != Float32
        A = Float32.(A)
    end

    # get whitening frequencies
    freqvec = rfftfreq(M,fs)
    left = findfirst(x -> x >= freqmin, freqvec)
    right = findlast(x -> x <= freqmax, freqvec)
    low, high = left - pad, right + pad

    if low <= 1
        low = 1
        left = low + pad
    end

    if high > length(freqvec)
        high = length(freqvec)- 1
        right = high - pad
    end

    # take fft and whiten
    fftraw = rfft(A,1)
    # left zero cut-off
     for jj = 1:N
         for ii = 1:low
            fftraw[ii,jj] = 0. + 0.0im
        end
    end
    # left tapering
     for jj = 1:N
         for ii = low+1:left
            fftraw[ii,jj] = cos(pi / 2 + pi / 2* (ii-low-1) / pad).^2 * exp(
            im * angle(fftraw[ii,jj]))
        end
    end
    # pass band
     for jj = 1:N
         for ii = left:right
            fftraw[ii,jj] = exp(im * angle(fftraw[ii,jj]))
        end
    end
    # right tapering
     for jj = 1:N
         for ii = right+1:high
            fftraw[ii,jj] = cos(pi/2 * (ii-right) / pad).^2 * exp(
            im * angle(fftraw[ii,jj]))
        end
    end
    # right zero cut-off
     for jj = 1:N
         for ii = high+1:size(fftraw,1)
            fftraw[ii,jj] = 0. + 0.0im
        end
    end
    return fftraw
end


end
