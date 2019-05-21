module Readfile
export initlogo, sconfig, ssourceloc, sreceiverloc, svelocity, sscatterloc, sattenuation
export readsource, readconfig
using LinearAlgebra, Dates

#Define data type
abstract type Abstractconfig end
abstract type Abstractsourceloc end
abstract type Abstractreceiverloc end
abstract type Abstractvelocity end
abstract type Abstractscatter end
abstract type Abstractattenuation end

#Define original structure
mutable struct sconfig <: Abstractconfig

    #define contents
    in_source ::AbstractString
    in_receiver ::AbstractString
    in_velocity ::AbstractString 
    in_scatter ::AbstractString 
    in_attenuation ::AbstractString
    rho ::AbstractArray{Float64} 
    NumofTime_ID ::Int64
    Time_ID ::AbstractArray{Int64}
    init_day ::AbstractArray{Int64}
    end_day ::AbstractArray{Int64}
    sampling_frequency ::Float64
    minimnum_frequency ::Float64
    IsAttenuation ::Bool
    IsSpectralNormalization ::Bool
    IsGaussianNoise ::Bool
    Sourcetype ::AbstractString
    Synthesize_unit_duration :: Float64
    T_for_gf :: Float64
    Source_randomseed ::AbstractArray{Int64}
    Source_amp_mean ::Float64
    Source_amp_variance ::Float64
    Source_peakfreq_mean ::Float64
    Source_peakfreq_variance ::Float64
    Stacking_Frequency ::Float64
    Source_average_num_per_unithour ::Int64
    C3_truncate_alpha ::Float64
    C3_max_time_lag :: Float64
    IsSaveRawTrace ::Bool
    IsSaveCC ::Bool

    sconfig() = new()

end

mutable struct ssourceloc <: Abstractsourceloc

    #define contents
    
    x ::AbstractArray{Float64}
    y ::AbstractArray{Float64}
    numofsource :: Int64
    ssourceloc() = new()

end

mutable struct sreceiverloc <: Abstractreceiverloc

    #define contents
    x ::AbstractArray{Float64}
    y ::AbstractArray{Float64}
    numofreceiver :: Int64
    sreceiverloc() = new()

end

mutable struct svelocity <: Abstractvelocity

    #define contents
    period ::AbstractArray{Float64}
    v ::AbstractArray{Float64}
    npt :: Int64
    svelocity() = new()

end

mutable struct sscatterloc <: Abstractscatter

    #define contents
    x ::AbstractArray{Float64}
    y ::AbstractArray{Float64}
    strength :: AbstractArray{Float64}
    numofscatter :: Int64
    sscatterloc() = new()

end

mutable struct sattenuation<: Abstractattenuation

    #define contents
    period ::AbstractArray{Float64}
    alpha ::AbstractArray{Float64}
    npt :: Int64
    sattenuation() = new()

end

"""
initlogo()
print initial logo
"""
function initlogo(;problem_name::String)


    print("
                    _     _            _      _____      _               _      
    /\\             | |   (_)          | |    / ____|    (_)             (_)     
   /  \\   _ __ ___ | |__  _  ___ _ __ | |_  | (___   ___ _ ___ _ __ ___  _  ___ 
  / /\\ \\ | '_ ` _ \\| '_ \\| |/ _ \\ '_ \\| __|  \\___ \\ / _ \\ / __| '_ ` _ \\| |/ __|
 / ____ \\| | | | | | |_) | |  __/ | | | |_   ____) |  __/ \\__ \\ | | | | | | (__ 
/_/    \\_\\_| |_| |_|_.__/|_|\\___|_| |_|\\__| |_____/ \\___|_|___/_| |_| |_|_|\\___|
                                                                                
  _   _       _             _____                           _             
 | \\ | |     (_)           / ____|                         | |            
 |  \\| | ___  _ ___  ___  | |  __  ___ _ __   ___ _ __ __ _| |_ ___  _ __ 
 | . ` |/ _ \\| / __|/ _ \\ | | |_ |/ _ \\ '_ \\ / _ \\ '__/ _` | __/ _ \\| '__|
 | |\\  | (_) | \\__ \\  __/ | |__| |  __/ | | |  __/ | | (_| | || (_) | |   
 |_| \\_|\\___/|_|___/\\___|  \\_____|\\___|_| |_|\\___|_|  \\__,_|\\__\\___/|_|   
                                                                         
                            v1.0.0 (Last update 05/21/2019)
                                 © Kurama Okubo
                                                      
")

    println("Job start running at "*string(now()))
    println("Problem name: "*problem_name)

end



"""
    readconfig(fname)

read config locations

# Arguments
- `fname`: input file name.
"""
function readconfig(fname::String)
    
    commentmark = '#'
    
	config = sconfig()
    
    open(fname) do f
		while !eof(f)
		
		lstr = readline(f)
		val = []

			if lstr[1] == commentmark
				continue

			else
				valtemp = split(lstr, "=")
				valname = strip(valtemp[1])

				val = strip(split(valtemp[2], commentmark)[1])

				#store valname			
                if valname == "in_source"
                    config.in_source = string(val)
                elseif valname == "in_receiver"
                    config.in_receiver = string(val)
                elseif valname == "in_velocity"
                    config.in_velocity = string(val)
                elseif valname == "in_scatter"
                    config.in_scatter = string(val)
                elseif valname == "in_attenuation"
                    config.in_attenuation = string(val)
                elseif valname == "NumofTime_ID"
                    config.NumofTime_ID = parse(Int64, val)
                elseif valname == "Time_ID"
                    config.Time_ID = zeros(Int64, config.NumofTime_ID)
                    for i = 1:config.NumofTime_ID
                        config.Time_ID[i] = parse(Int64, split(val)[i])
                    end

                elseif valname == "init_day"
                    config.init_day = zeros(Int64, config.NumofTime_ID)
                    for i = 1:config.NumofTime_ID
                        config.init_day[i] = parse(Int64, split(val)[i])
                    end

                elseif valname == "end_day"
                    config.end_day = zeros(Int64, config.NumofTime_ID)
                    for i = 1:config.NumofTime_ID
                        config.end_day[i] = parse(Int64, split(val)[i])
                    end

                elseif valname == "rho"
                    config.rho = zeros(Float64, config.NumofTime_ID)
                    for i = 1:config.NumofTime_ID
                        config.rho[i] = parse(Float64, split(val)[i])
                    end

                elseif valname == "sampling_frequency"
                    config.sampling_frequency = parse(Float64, val)
                elseif valname == "minimnum_frequency"
                    config.minimnum_frequency = parse(Float64, val)
                elseif valname == "IsAttenuation"
                    config.IsAttenuation = parse(Bool, val)
                elseif valname == "IsSpectralNormalization"
                    config.IsSpectralNormalization = parse(Bool, val)
                elseif valname == "IsGaussianNoise"
                    config.IsGaussianNoise = parse(Bool, val)
                elseif valname == "Sourcetype"
                    config.Sourcetype = string(val)
                elseif valname == "Synthesize_unit_duration"
                    config.Synthesize_unit_duration = parse(Float64, val)
                elseif valname == "T_for_gf"
                    config.T_for_gf = parse(Float64, val)
                elseif valname == "Source_randomseed"
                    config.Source_randomseed = zeros(Int64, config.NumofTime_ID)
                    for i = 1:config.NumofTime_ID
                        config.Source_randomseed[i] = parse(Int64, split(val)[i])
                    end
                elseif valname == "Source_amp_mean"
                    config.Source_amp_mean = parse(Float64, val)
                elseif valname == "Source_amp_variance"
                    config.Source_amp_variance = parse(Float64, val)
                elseif valname == "Source_peakfreq_mean"
                    config.Source_peakfreq_mean = parse(Float64, val)
                elseif valname == "Source_peakfreq_variance"
                    config.Source_peakfreq_variance = parse(Float64, val)
                elseif valname == "Stacking_Frequency"
                    config.Stacking_Frequency = parse(Float64, val)
                elseif valname == "Source_average_num_per_unithour"
                    config.Source_average_num_per_unithour = parse(Int64, val)
                elseif valname == "C3_truncate_alpha"
                    config.C3_truncate_alpha = parse(Float64, val)
                elseif valname == "C3_max_time_lag"
                    config.C3_max_time_lag = parse(Float64, val)
                elseif valname == "IsSaveRawTrace"
                    config.IsSaveRawTrace = parse(Bool, val)
                elseif valname == "IsSaveCC"
                    config.IsSaveCC = parse(Bool, val)

				else
					println("valuable \""*valname*"\" does not exist. skip this valuable.")
				end
			
			end
		end
	end

    return config

end


"""
    readsource(fname)

read source locations

# Arguments
- `fname`: input file name.
"""
function readsource(fname::String)

    commentmark = '#'
    
    sourceloc = ssourceloc()
    numsource = 0

    sourceloc.x = zeros(0)
    sourceloc.y = zeros(0)

    open(fname) do f
        while !eof(f)

            lstr = readline(f)
            #println(lstr)
            if lstr[1] == commentmark
                continue

            else
                sourceloc.x = append!(sourceloc.x, parse(Float32, (split(lstr, ',')[1])))
                sourceloc.y = append!(sourceloc.y, parse(Float32, (split(lstr, ',')[2])))
            end

            numsource += 1
        end
    end

    sourceloc.numofsource = numsource

    return sourceloc    

end

"""
    readreceiver(fname)

read receiver locations

# Arguments
- `fname`: input file name.
"""
function readreceiver(fname::String)

    commentmark = '#'
    
    receiverloc = sreceiverloc()
    numreceiver = 0

    receiverloc.x = zeros(0)
    receiverloc.y = zeros(0)

    open(fname) do f
        while !eof(f)

            lstr = readline(f)
            #println(lstr)
            if lstr[1] == commentmark
                continue

            else
                receiverloc.x = append!(receiverloc.x, parse(Float64, (split(lstr, ',')[1])))
                receiverloc.y = append!(receiverloc.y, parse(Float64, (split(lstr, ',')[2])))
            end

            numreceiver += 1
        end
    end

    receiverloc.numofreceiver = numreceiver

    return receiverloc    

end

"""
    readvelocity(fname)

read dispersion curve

# Arguments
- `fname`: input file name.
"""
function readvelocity(fname::String)

    commentmark = '#'
    
    velocity = svelocity()
    npt = 0

    velocity.period = zeros(0)
    velocity.v = zeros(0)

    open(fname) do f
        while !eof(f)

            lstr = readline(f)
            #println(lstr)
            if lstr[1] == commentmark
                continue

            else
                velocity.period = append!(velocity.period, parse(Float64, (split(lstr, ',')[1])))
                velocity.v = append!(velocity.v, parse(Float64, (split(lstr, ',')[2])))
            end

            npt += 1
        end
    end

    velocity.npt = npt

    return velocity

end


"""
    readscatter(fname)

read scatter locations

# Arguments
- `fname`: input file name.
"""
function readscatter(fname::String)

    commentmark = '#'
    
    scatterloc = sscatterloc()
    numofscatter = 0

    scatterloc.x = zeros(0)
    scatterloc.y = zeros(0)
    scatterloc.strength = zeros(0)

    open(fname) do f
        while !eof(f)

            lstr = readline(f)
            #println(lstr)
            if lstr[1] == commentmark
                continue

            else
                scatterloc.x = append!(scatterloc.x, parse(Float64, (split(lstr, ',')[1])))
                scatterloc.y = append!(scatterloc.y, parse(Float64, (split(lstr, ',')[2])))
                scatterloc.strength = append!(scatterloc.strength, parse(Float64, (split(lstr, ',')[3])))
            end

            numofscatter += 1
        end
    end

    scatterloc.numofscatter = numofscatter

    return scatterloc

end



"""
    readattenuation(fname)

read attenuation

# Arguments
- `fname`: input file name.
"""
function readattenuation(fname::String)

    commentmark = '#'
    
    attenuation = sattenuation()
    npt = 0

    attenuation.period = zeros(0)
    attenuation.alpha = zeros(0)

    open(fname) do f
        while !eof(f)

            lstr = readline(f)
            #println(lstr)
            if lstr[1] == commentmark
                continue

            else
                attenuation.period = append!(attenuation.period, parse(Float64, (split(lstr, ',')[1])))
                attenuation.alpha = append!(attenuation.alpha, parse(Float64, (split(lstr, ',')[2])))
            end

            npt += 1
        end
    end

    attenuation.npt = npt

    return attenuation

end
end
