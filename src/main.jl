"""
Ambient noise generator
April 24, 2019
Kurama Okubo
"""

#---Include original libraries---#
include("./lib/readfile.jl")
include("./lib/gfandsource.jl")
include("./lib/receiverdata.jl")
include("./lib/synthesis_wave.jl")

using Revise
using Printf, LinearAlgebra, FFTW, PlotlyJS, Dierckx, Random, Distributions, ORCA, DSP,  HDF5, BenchmarkTools, ProgressMeter
import .Readfile, .GFandSource, .ReceiverData, .Synthesis_wave


function main(;problem_name::String)

    if isempty(problem_name)
        println("Enter problem name:\n")
        problem_name = readline()
    end

    Readfile.initlogo(problem_name=problem_name)

    #Open HDF5
    outputdir = problem_name*"/OUTPUT_FILES"
    if ~isdir(outputdir);  mkdir(outputdir);  end


    Figdirroot = problem_name * "/OUTPUT_FILES/figs"
    if ~isdir(Figdirroot) mkdir(Figdirroot) end

    #read config, source and receivers
    config      = Readfile.readconfig(problem_name*"/inputfiles/config.in")
    sourceloc   = Readfile.readsource(problem_name*"/inputfiles/"*config.in_source);
    receiverloc = Readfile.readreceiver(problem_name*"/inputfiles/"*config.in_receiver);

    hdf5_o_C1 = h5open(outputdir*"/CC1.h5", "w")
    hdf5_o_C2 = h5open(outputdir*"/CC2.h5", "w")
    hdf5_o_C3 = h5open(outputdir*"/CC3.h5", "w")
    hdf5_o_R = h5open(outputdir*"/Receiver.h5", "w")

    #form pair ids for C1
    NumofRecPair = Int(receiverloc.numofreceiver * (receiverloc.numofreceiver-1)/2)
    recpair = zeros(Int64,NumofRecPair, 2)
    icount = 0
    for i = 1:receiverloc.numofreceiver
        for j = i:receiverloc.numofreceiver
            if i == j
                continue;
            else
                icount += 1
                recpair[icount, 1] = i
                recpair[icount, 2] = j
            end
        end
    end

    #make receiver group for C2 and C3
    NumofRecGroup = Int((receiverloc.numofreceiver * (receiverloc.numofreceiver-1) * (receiverloc.numofreceiver-2) / 6) * 3) #nC3 * 3C1
    recgroup= zeros(Int64, NumofRecGroup, 3) #virtual source, R1, R2
    icount = 0

    for i = 1:receiverloc.numofreceiver
        for j = i:receiverloc.numofreceiver
            for k = j:receiverloc.numofreceiver
                if i == j || j==k || k ==i
                    continue;
                else
                    icount += 1
                    recgroup[icount, 1] = i
                    recgroup[icount, 2] = j
                    recgroup[icount, 3] = k

                    icount += 1
                    recgroup[icount, 1] = k
                    recgroup[icount, 2] = i
                    recgroup[icount, 3] = j
                    
                    icount += 1
                    recgroup[icount, 1] = j
                    recgroup[icount, 2] = k
                    recgroup[icount, 3] = i

                end
            end
        end
    end

    #Time info

    NumofSynthesize_unit = round(Int64, 24 * 60 * 60 / config.Synthesize_unit_duration) 

    T = config.T_for_gf # compute time for green's function between source and receivers [s]
    Synthesize_unit_duration = config.Synthesize_unit_duration; #every this unit the signal is synthesized

    fs = config.sampling_frequency
    dt = 1.0 / fs
    NT = round(Int64, T/dt)
    omega = 2*pi.*(fs*(0:(NT/2))/NT); #0 -> fnyq
    t_trace = dt .* range(0, stop=round(Int64, Synthesize_unit_duration./dt)-1)
    #define random source activate time: this is constant for all receivers
    maxT = Synthesize_unit_duration - T #avoid truncated at the end of unit duration 

    t_cc = dt .* range(-(length(t_trace)-1), stop=length(t_trace)-1)

    write(hdf5_o_C1, "Receiver_pair", recpair)
    write(hdf5_o_C1, "Lag_time", collect(t_cc))
    write(hdf5_o_C2, "Receiver_group", recgroup)
    write(hdf5_o_C2, "Lag_time", collect(t_cc))
    write(hdf5_o_C3, "Receiver_group", recgroup)
    write(hdf5_o_C3, "Lag_time", collect(t_cc))
    write(hdf5_o_R, "t", collect(t_trace))

    #loop time id
    #for Time_ID = config.Time_ID
    for Time_ID = 1
        #Time_ID = 1
        println("#-----------------#")
        println(@sprintf("Time ID%2d", Time_ID))
        println("#-----------------#")

        dname = @sprintf("TimeID_%02d", Time_ID)
        modelpath = problem_name*"/inputfiles/"*dname*"/"

        #read dispersion curve
        velocity = Readfile.readvelocity(modelpath*config.in_velocity);
        vel_omega = 2*pi ./ velocity.period;

        #read scatter location
        scatterloc = Readfile.readscatter(modelpath*config.in_scatter);

        #read attenuation curve
        attenuation = Readfile.readattenuation(modelpath*config.in_attenuation);

        #plot source, receiver and scatter
        trace1 = scatter(;x=sourceloc.x./1e3, y=sourceloc.y./1e3, mode="markers", name="source")
        textrec = []
        for i = 1:receiverloc.numofreceiver
        	push!(textrec, string(i))
        end

        trace2 = scatter(;x=receiverloc.x./1e3, y=receiverloc.y./1e3, mode="markers+text", name="receiver",
        	textposition="top center",
          	text=textrec)
        trace3 = scatter(;x=scatterloc.x./1e3, y=scatterloc.y./1e3, mode="markers", name="scatter")

        layout = Layout(width=500,height=500,
        	xaxis_range=[-300, 300], yaxis_range=[-300, 300])
        p = plot([trace1, trace2, trace3], layout)
        savefig(p, Figdirroot*"/sourceandreceiverloc.png")

        #loop day
        NumofDay_per_TimeID = config.end_day[Time_ID] - config.init_day[Time_ID] + 1;

        #Create random seed for each synthesize unit
        rng = MersenneTwister(config.Source_randomseed[Time_ID]);
        randseeds = rand(rng, 1:999999, NumofDay_per_TimeID*NumofSynthesize_unit);

        #for Day_ID = config.init_day[Time_ID]:config.end_day[Time_ID]
        for Day_ID = config.init_day[Time_ID]:2
            
            #Day_ID = 1 # day
            println("#-----------------#")
            println(@sprintf("Day %4d", Day_ID))
            println("#-----------------#")

            #------------------------#
            #---Source Randomness---#
            #------------------------#
           
            #Progress Bar
            prog = Progress(12, 1.0)

            #for Synthesize_unit_ID = 1:NumofSynthesize_unit
            for Synthesize_unit_ID = 1:12

                Random.seed!(randseeds[(Day_ID-1)*NumofSynthesize_unit+Synthesize_unit_ID]); #This should be random for each Syntesize unit duration
                source_activatetime = rand(Uniform(0, maxT), sourceloc.numofsource, config.Source_average_num_per_unithour);

                #Source fM is distributed Gaussian in log scale
                rng = Normal(log10(config.Source_peakfreq_mean), abs(0.01*config.Source_peakfreq_variance*log10(config.Source_peakfreq_mean)))
                Random.seed!(randseeds[(Day_ID-1)*NumofSynthesize_unit+Synthesize_unit_ID]);
                fMrand = 10.0.^(rand(rng, sourceloc.numofsource * config.Source_average_num_per_unithour)[:]);
                write(hdf5_o_C1, @sprintf("fMrand/TimeID%02d/Day%04d/UnitID%04d/",Time_ID, Day_ID, Synthesize_unit_ID), fMrand)

                receiverdata = Array{ReceiverData.sreceiverdata, 1}(undef, receiverloc.numofreceiver);

                for i = 1:receiverloc.numofreceiver
                    receiverdata[i] = ReceiverData.sreceiverdata();
                    #Synthesize noise field
                    Synthesis_wave.synthesize_noise_field!(receiverdata=receiverdata[i], omega=omega, t_trace=t_trace, source_activatetime=source_activatetime,
                        Receiver_ID = i, config=config, rho=config.rho[Time_ID], velocity=velocity, attenuation=attenuation, sourceloc = sourceloc,
                        receiverloc=receiverloc, scatterloc=scatterloc, fMrand=fMrand);

                    name_o = @sprintf("TimeID%02d/Day%04d/UnitID%04d/U.%02d", Time_ID, Day_ID, Synthesize_unit_ID, receiverdata[i].Receiver_ID)

                    g1 = g_create(hdf5_o_R, name_o)

                    g1["U"] = receiverdata[i].u
                    attrs(g1)["fs"] = fs
                    attrs(g1)["r1x"] = receiverdata[i].loc_x
                    attrs(g1)["r1y"] = receiverdata[i].loc_y
                    
                end

                #initialize C1 struct
                C1data =Array{ReceiverData.sC1data, 1}(undef, NumofRecPair)

                for i = 1:NumofRecPair

                	C1data[i] = ReceiverData.sC1data()

                	#define contents
                    C1data[i].Receiver_1 = recpair[i,1]
                    C1data[i].Receiver_2 = recpair[i,2]
                    C1data[i].r1_loc_x = receiverloc.x[C1data[i].Receiver_1]
                    C1data[i].r1_loc_y = receiverloc.y[C1data[i].Receiver_1]
                    C1data[i].r2_loc_x = receiverloc.x[C1data[i].Receiver_2]
                    C1data[i].r2_loc_y = receiverloc.y[C1data[i].Receiver_2]

                    C1data[i].dist     = norm([C1data[i].r1_loc_x, C1data[i].r1_loc_y] .- [C1data[i].r2_loc_x, C1data[i].r2_loc_y])
                    C1data[i].nt    = length(t_trace)

                    #compute cc1
                    #filtering traces
                    rec1_u_origin = receiverdata[C1data[i].Receiver_1].u 
                    rec2_u_origin = receiverdata[C1data[i].Receiver_2].u 


                    #--------------------------------#
                    #---Filtering original traces----#
                    #--------------------------------#

                    #bandpass_minfreq = 0.01
                    #bandpass_maxfreq = 0.2

                	#bpfilter = digitalfilter(Bandpass(bandpass_minfreq, bandpass_maxfreq; fs=fs),Butterworth(2));

                	#rec1_u_origin = filtfilt(bpfilter, rec1_u_origin)
                	#rec2_u_origin = filtfilt(bpfilter, rec2_u_origin)

                    #--------------------------------#
                    #--------------------------------#
                    rec1_u = rec1_u_origin
                    rec2_u = rec2_u_origin
                    
                    #fft

                    rec1_FU = fft(rec1_u)
                    rec2_FU = fft(rec2_u)

                	signal_magnification = 1e6
                	minplotT = -T/2
                	maxplotT = T/2


                	"""
                	Exercise 1: First order cross correlation C1:1->2
                	"""

                	signal_magnification = 1e6

                	cc1_12_pos = zeros(Complex{Float64}, C1data[i].nt)
                	cc1_12_neg = zeros(Complex{Float64}, C1data[i].nt)

                	ccid = -(C1data[i].nt-1):C1data[i].nt-1

                	#modified summation

                    if config.IsSpectralNormalization
                        #Whiten signal
                        freqmin = 0.01;
                        freqmax = 2.0;

                        #println(length(rec1_FU ))
                        #rec1_FU1 = GFandSource.whiten(real.(ifft(rec1_FU)), freqmin, freqmax, fs, pad=50)
                        #println(length(rec1_FU1))


                        #trace1 = scatter(;y = abs.(rec1_FU), mode="lines")
                        #trace2 = scatter(;y = abs.(rec1_FU1), mode="lines")

                        #p = plot([trace2])
                        #savefig(p, Figdirroot*"/whitentest.png")

                        #rec2_FU = GFandSource.whiten(real.(ifft(rec2_FU)), freqmin, freqmax, fs, pad=50)

                        for j = 1:C1data[i].nt
                                
                    	        cc1_12_pos[j] = conj(rec1_FU[j]) * rec2_FU[j] / (abs(rec1_FU[j]) * abs(rec2_FU[j]))
                    	        cc1_12_neg[j] = rec1_FU[j] * conj(rec2_FU[j]) / (abs(rec1_FU[j]) * abs(rec2_FU[j]))
                        end

                    else
                        
                        for j = 1:C1data[i].nt

                                cc1_12_pos[j] = conj(rec1_FU[j]) * rec2_FU[j] 
                                cc1_12_neg[j] = rec1_FU[j] * conj(rec2_FU[j])

                        end
                    end

                	# plot along azimuth
                	t_cc = dt .* range(-(C1data[i].nt-1), stop=C1data[i].nt-1)

                	cc1_12_temp = vcat(ifft(cc1_12_neg)[end:-1:2], ifft(cc1_12_pos)) 
                	#signal_magnification = 5*1/maximum(real.(cc1_12_temp))

                    C1data[i].lag_t = t_cc
                    C1data[i].cc1 = real.(cc1_12_temp)

                end

              
                #-----Save CC1----#
                #recpair
                for i=1:NumofRecPair
                        
                    name_o = @sprintf("TimeID%02d/Day%04d/UnitID%04d/CC1.%02d-%02d", Time_ID, Day_ID, Synthesize_unit_ID, C1data[i].Receiver_1,   C1data[i].Receiver_2)

                    g1 = g_create(hdf5_o_C1, name_o)

                    g1["CC1"] = C1data[i].cc1 

                    attrs(g1)["fs"] = fs
                    attrs(g1)["Receiver_1"] = C1data[i].Receiver_1
                    attrs(g1)["Receiver_2"] = C1data[i].Receiver_2
                    attrs(g1)["r1x"] = C1data[i].r1_loc_x
                    attrs(g1)["r1y"] = C1data[i].r1_loc_y
                    attrs(g1)["r2x"] = C1data[i].r2_loc_x
                    attrs(g1)["r2y"] = C1data[i].r2_loc_y

                end

                """
                Exercise 2: Second order cross correlation C2:1->2
                """

                #initialize C1 struct
                C2data =Array{ReceiverData.sC2data, 1}(undef, NumofRecGroup)

                #println("Compute C2")
                for i = 1:NumofRecGroup

                	C2data[i] = ReceiverData.sC2data()

                	#define contents
                    C2data[i].Receiver_v = recgroup[i,1]
                    C2data[i].Receiver_1 = recgroup[i,2]
                    C2data[i].Receiver_2 = recgroup[i,3]
                    C2data[i].rv_loc_x = receiverloc.x[C2data[i].Receiver_v]
                    C2data[i].rv_loc_y = receiverloc.y[C2data[i].Receiver_v]
                    C2data[i].r1_loc_x = receiverloc.x[C2data[i].Receiver_1]
                    C2data[i].r1_loc_y = receiverloc.y[C2data[i].Receiver_1]
                    C2data[i].r2_loc_x = receiverloc.x[C2data[i].Receiver_2]
                    C2data[i].r2_loc_y = receiverloc.y[C2data[i].Receiver_2]

                    C2data[i].dist_v_1     = norm([C2data[i].rv_loc_x, C2data[i].rv_loc_y] .- [C2data[i].r1_loc_x, C2data[i].r1_loc_y])
                    C2data[i].dist_v_2     = norm([C2data[i].rv_loc_x, C2data[i].rv_loc_y] .- [C2data[i].r2_loc_x, C2data[i].r2_loc_y])
                    C2data[i].dist_1_2     = norm([C2data[i].r1_loc_x, C2data[i].r1_loc_y] .- [C2data[i].r2_loc_x, C2data[i].r2_loc_y])
                    C2data[i].nt    = length(t_trace)

                	t_cc = dt .* range(-(C2data[i].nt-1), stop=C2data[i].nt-1)
                    C2data[i].lag_t = t_cc 

                    #compute cc1
                    #filtering traces
                    recv_u_origin = receiverdata[C2data[i].Receiver_v].u 
                    rec1_u_origin = receiverdata[C2data[i].Receiver_1].u 
                    rec2_u_origin = receiverdata[C2data[i].Receiver_2].u 

                    #--------------------------------#
                    #---Filtering original traces----#
                    #--------------------------------#

                    #bandpass_minfreq = 0.01
                    #bandpass_maxfreq = 0.2

                	#bpfilter = digitalfilter(Bandpass(bandpass_minfreq, bandpass_maxfreq; fs=fs),Butterworth(2));

                	#recv_u_origin = filtfilt(bpfilter, recv_u_origin)
                	#rec1_u_origin = filtfilt(bpfilter, rec1_u_origin)
                	#rec2_u_origin = filtfilt(bpfilter, rec2_u_origin)

                    #--------------------------------#
                    #--------------------------------#
                    recv_u = recv_u_origin
                    rec1_u = rec1_u_origin
                    rec2_u = rec2_u_origin
                    
                    #fft

                    recv_FU = fft(recv_u)
                    rec1_FU = fft(rec1_u)
                    rec2_FU = fft(rec2_u)

                	#CC1: v->1, v-2
                	cc1_v1_pos = zeros(Complex{Float64}, C2data[i].nt)
                	cc1_v1_neg = zeros(Complex{Float64}, C2data[i].nt)
                	cc1_v2_pos = zeros(Complex{Float64}, C2data[i].nt)
                	cc1_v2_neg = zeros(Complex{Float64}, C2data[i].nt)


                    if config.IsSpectralNormalization

                        for j = 1:C2data[i].nt
                            
                            cc1_v1_pos[j] = conj(recv_FU[j]) * rec1_FU[j]
                            cc1_v1_neg[j] = recv_FU[j] * conj(rec1_FU[j])

                            cc1_v2_pos[j] = conj(recv_FU[j]) * rec2_FU[j]
                            cc1_v2_neg[j] = recv_FU[j] * conj(rec2_FU[j])
                        end

                    else

                        for j = 1:C2data[i].nt
                        
                        cc1_v1_pos[j] = conj(recv_FU[j]) * rec1_FU[j]
                        cc1_v1_neg[j] = recv_FU[j] * conj(rec1_FU[j])

                        cc1_v2_pos[j] = conj(recv_FU[j]) * rec2_FU[j]
                        cc1_v2_neg[j] = recv_FU[j] * conj(rec2_FU[j])
                        end
                    end

                    cc1_v1_temp = vcat(ifft(cc1_v1_neg)[end:-1:2], ifft(cc1_v1_pos)) 
                    cc1_v2_temp = vcat(ifft(cc1_v2_neg)[end:-1:2], ifft(cc1_v2_pos)) 

                   
                    """
                    Stack positive side C2POS and negative side C2NEG
                    """

                    #CC2POS and CC2NEG
                    cc2POS_12_pos = zeros(Complex{Float64}, C2data[i].nt)
                    cc2POS_12_neg = zeros(Complex{Float64}, C2data[i].nt)

                    cc2NEG_12_pos = zeros(Complex{Float64}, C2data[i].nt)
                    cc2NEG_12_neg = zeros(Complex{Float64}, C2data[i].nt)

                    
                    if config.IsSpectralNormalization

                        for j = 1:C2data[i].nt
                            
                            cc2POS_12_pos[j] = conj(cc1_v1_pos[j]) * cc1_v2_pos[j] / (abs(cc1_v1_pos[j]) * abs(cc1_v2_pos[j]))
                            cc2POS_12_neg[j] = cc1_v1_pos[j] * conj(cc1_v2_pos[j]) / (abs(cc1_v1_pos[j]) * abs(cc1_v2_pos[j]))

                    		cc2NEG_12_pos[j] = conj(cc1_v1_neg[j]) * cc1_v2_neg[j] / (abs(cc1_v1_neg[j]) * abs(cc1_v2_neg[j]))
                            cc2NEG_12_neg[j] = cc1_v1_neg[j] * conj(cc1_v2_neg[j]) / (abs(cc1_v1_neg[j]) * abs(cc1_v2_neg[j]))
                        end

                    else
                        for j = 1:C2data[i].nt
                            
                            cc2POS_12_pos[j] = conj(cc1_v1_pos[j]) * cc1_v2_pos[j] 
                            cc2POS_12_neg[j] = cc1_v1_pos[j] * conj(cc1_v2_pos[j])

                            cc2NEG_12_pos[j] = conj(cc1_v1_neg[j]) * cc1_v2_neg[j]
                            cc2NEG_12_neg[j] = cc1_v1_neg[j] * conj(cc1_v2_neg[j])
                        end
                    end


                	cc2_POStemp = real.(vcat(ifft(cc2POS_12_neg)[end:-1:2], ifft(cc2POS_12_pos)))
                	cc2_NEGtemp = real.(vcat(ifft(cc2NEG_12_neg)[end:-1:2], ifft(cc2NEG_12_pos)))
                	cc2_POSandNEG = cc2_POStemp .+ cc2_NEGtemp

                    C2data[i].cc2_POS = cc2_POStemp
                    C2data[i].cc2_NEG = cc2_NEGtemp
                    C2data[i].cc2_ALL = cc2_POSandNEG

                    #-----Save CC2----#

                    name_o = @sprintf("TimeID%02d/Day%04d/UnitID%04d/CC2.V%02d-%02d-%02d", Time_ID, Day_ID, Synthesize_unit_ID, C2data[i].Receiver_v, C2data[i].Receiver_1,   C2data[i].Receiver_2)

                    g1 = g_create(hdf5_o_C2, name_o)

                    g1["CC2"] = C2data[i].cc2_ALL

                    attrs(g1)["fs"] = fs
                    attrs(g1)["Receiver_v"] = C2data[i].Receiver_v
                    attrs(g1)["Receiver_1"] = C2data[i].Receiver_1
                    attrs(g1)["Receiver_2"] = C2data[i].Receiver_2
                    attrs(g1)["rvx"] = C2data[i].rv_loc_x
                    attrs(g1)["rvy"] = C2data[i].rv_loc_y
                    attrs(g1)["r1x"] = C2data[i].r1_loc_x
                    attrs(g1)["r1y"] = C2data[i].r1_loc_y
                    attrs(g1)["r2x"] = C2data[i].r2_loc_x
                    attrs(g1)["r2y"] = C2data[i].r2_loc_y

                end

                """
                Exercise 3: C3 
                """
                #initialize C3 struct
                C3data =Array{ReceiverData.sC3data, 1}(undef, NumofRecGroup)
                vref = mean(velocity.v) #truncate ballistic wave with t = alpha * distance / vref

                for i = 1:NumofRecGroup

                    C3data[i] = ReceiverData.sC3data()

                    #define contents
                    C3data[i].Receiver_v = recgroup[i,1]
                    C3data[i].Receiver_1 = recgroup[i,2]
                    C3data[i].Receiver_2 = recgroup[i,3]
                    C3data[i].rv_loc_x = receiverloc.x[C3data[i].Receiver_v]
                    C3data[i].rv_loc_y = receiverloc.y[C3data[i].Receiver_v]
                    C3data[i].r1_loc_x = receiverloc.x[C3data[i].Receiver_1]
                    C3data[i].r1_loc_y = receiverloc.y[C3data[i].Receiver_1]
                    C3data[i].r2_loc_x = receiverloc.x[C3data[i].Receiver_2]
                    C3data[i].r2_loc_y = receiverloc.y[C3data[i].Receiver_2]

                    C3data[i].dist_v_1     = norm([C3data[i].rv_loc_x, C3data[i].rv_loc_y] .- [C3data[i].r1_loc_x, C3data[i].r1_loc_y])
                    C3data[i].dist_v_2     = norm([C3data[i].rv_loc_x, C3data[i].rv_loc_y] .- [C3data[i].r2_loc_x, C3data[i].r2_loc_y])
                    C3data[i].dist_1_2     = norm([C3data[i].r1_loc_x, C3data[i].r1_loc_y] .- [C3data[i].r2_loc_x, C3data[i].r2_loc_y])
                    C3data[i].nt    = length(t_trace)

                    t_cc = dt .* range(-(C3data[i].nt-1), stop=C3data[i].nt-1)
                    C3data[i].lag_t = t_cc 

                    #compute cc1
                    #filtering traces
                    recv_u_origin = receiverdata[C3data[i].Receiver_v].u 
                    rec1_u_origin = receiverdata[C3data[i].Receiver_1].u 
                    rec2_u_origin = receiverdata[C3data[i].Receiver_2].u 

                    #--------------------------------#
                    #---Filtering original traces----#
                    #--------------------------------#

                    #bandpass_minfreq = 0.01
                    #bandpass_maxfreq = 0.2

                    #bpfilter = digitalfilter(Bandpass(bandpass_minfreq, bandpass_maxfreq; fs=fs),Butterworth(2));

                    #recv_u_origin = filtfilt(bpfilter, recv_u_origin)
                    #rec1_u_origin = filtfilt(bpfilter, rec1_u_origin)
                    #rec2_u_origin = filtfilt(bpfilter, rec2_u_origin)

                    #--------------------------------#
                    #--------------------------------#
                    recv_u = recv_u_origin
                    rec1_u = rec1_u_origin
                    rec2_u = rec2_u_origin
                    
                    #fft

                    recv_FU = fft(recv_u)
                    rec1_FU = fft(rec1_u)
                    rec2_FU = fft(rec2_u)

                    #CC1: v->1, v-2
                    cc1_v1_pos = zeros(Complex{Float64}, C3data[i].nt)
                    cc1_v1_neg = zeros(Complex{Float64}, C3data[i].nt)
                    cc1_v2_pos = zeros(Complex{Float64}, C3data[i].nt)
                    cc1_v2_neg = zeros(Complex{Float64}, C3data[i].nt)


                    if config.IsSpectralNormalization

                        for j = 1:C3data[i].nt
                            
                            cc1_v1_pos[j] = conj(recv_FU[j]) * rec1_FU[j] / (abs(recv_FU[j]) * abs(rec1_FU[j]))
                            cc1_v1_neg[j] = recv_FU[j] * conj(rec1_FU[j]) / (abs(recv_FU[j]) * abs(rec1_FU[j]))

                            cc1_v2_pos[j] = conj(recv_FU[j]) * rec2_FU[j] / (abs(recv_FU[j]) * abs(rec2_FU[j]))
                            cc1_v2_neg[j] = recv_FU[j] * conj(rec2_FU[j]) / (abs(recv_FU[j]) * abs(rec2_FU[j]))

                            #cc1_v1_pos[j] = conj(recv_FU[j]) * rec1_FU[j] 
                            #cc1_v1_neg[j] = recv_FU[j] * conj(rec1_FU[j]) 

                            #cc1_v2_pos[j] = conj(recv_FU[j]) * rec2_FU[j] 
                            #cc1_v2_neg[j] = recv_FU[j] * conj(rec2_FU[j]) 

                        end

                    else
                        for j = 1:C3data[i].nt
                            
                            cc1_v1_pos[j] = conj(recv_FU[j]) * rec1_FU[j]
                            cc1_v1_neg[j] = recv_FU[j] * conj(rec1_FU[j])

                            cc1_v2_pos[j] = conj(recv_FU[j]) * rec2_FU[j]
                            cc1_v2_neg[j] = recv_FU[j] * conj(rec2_FU[j])
                        end
                    end


                    cc1_v1_temp = real.(vcat(ifft(cc1_v1_neg)[end:-1:2], ifft(cc1_v1_pos)))
                    cc1_v2_temp = real.(vcat(ifft(cc1_v2_neg)[end:-1:2], ifft(cc1_v2_pos)))

                    #--------------------------------#
                    #---Truncating ballistic wave----#
                    #--------------------------------#

                    #find truncate time for cc_v1
                    t_trancate = config.C3_truncate_alpha * C3data[i].dist_v_1 / vref
                    tid_v_1_init = findfirst(x -> x>=t_trancate, t_cc[round(Int64, length(t_cc)/2):end]);
                    tid_v_1_end = findfirst(x -> x>=config.C3_max_time_lag, t_cc[round(Int64, length(t_cc)/2):end]);

                    cc1_v1_pos_windowing = real.(ifft(cc1_v1_pos))
                    cc1_v1_neg_windowing = real.(ifft(cc1_v1_neg))

                    trace1 =  scatter(;y=cc1_v1_pos_windowing, mode="lines", name="whole CC1:POS")
                    trace2 =  scatter(;y=cc1_v1_neg_windowing, mode="lines", name="Coda CC1:NEG")
                    layout = Layout(width=800,height=300)

                    data = [trace1, trace2]
                    p = plot(data, layout)


                    cc1_v1_pos_windowing[1:tid_v_1_init] = zeros(length(1:tid_v_1_init))
                    cc1_v1_neg_windowing[1:tid_v_1_init] = zeros(length(1:tid_v_1_init))

                    cc1_v1_pos_windowing[tid_v_1_end:end] = zeros(length(cc1_v1_pos_windowing[tid_v_1_end:end]))
                    cc1_v1_neg_windowing[tid_v_1_end:end] = zeros(length(cc1_v1_neg_windowing[tid_v_1_end:end]))

                    #check windowed signal
                    cc1_v1_temp_windowed = real.(vcat(cc1_v1_neg_windowing[end:-1:2], cc1_v1_pos_windowing))

                    cc1_v1_pos = fft(cc1_v1_pos_windowing)
                    cc1_v1_neg = fft(cc1_v1_neg_windowing)

                    #find truncate time for cc_v2
                    t_trancate = config.C3_truncate_alpha * C3data[i].dist_v_2 / vref
                    tid_v_2_init = findfirst(x -> x>=t_trancate, t_cc[round(Int64, length(t_cc)/2):end]);
                    tid_v_2_end = findfirst(x -> x>=config.C3_max_time_lag, t_cc[round(Int64, length(t_cc)/2):end]);

                    cc1_v2_pos_windowing = real.(ifft(cc1_v2_pos))
                    cc1_v2_neg_windowing = real.(ifft(cc1_v2_neg))

                    cc1_v2_pos_windowing[1:tid_v_2_init] = zeros(length(1:tid_v_2_init))
                    cc1_v2_neg_windowing[1:tid_v_2_init] = zeros(length(1:tid_v_2_init))

                    cc1_v2_pos_windowing[tid_v_2_end:end] = zeros(length(cc1_v2_pos_windowing[tid_v_2_end:end]))
                    cc1_v2_neg_windowing[tid_v_2_end:end] = zeros(length(cc1_v2_neg_windowing[tid_v_2_end:end]))

                    #check windowed signal
                    cc1_v2_temp_windowed = real.(vcat(cc1_v2_neg_windowing[end:-1:2], cc1_v2_pos_windowing))
                    
                    cc1_v2_pos = fft(cc1_v2_pos_windowing)
                    cc1_v2_neg = fft(cc1_v2_neg_windowing)

                    #Check trace
                    if C3data[i].Receiver_v == 1
                        #test plot when virtual source is No. 8
                        #Plot u
                        trace1 =  scatter(;x=t_cc, y=cc1_v1_temp, mode="lines", name="whole CC1:V-R1")
                        trace2 =  scatter(;x=t_cc, y=cc1_v1_temp_windowed, mode="lines", name="Coda CC1:V-R1")
                        layout = Layout(width=800,height=300)

                        data = [trace1, trace2]
                        p = plot(data, layout)

                    end


                    """
                    Stack positive side C3POS and negative side C3NEG
                    """

                    #CC2POS and CC2NEG
                    cc3POS_12_pos = zeros(Complex{Float64}, C3data[i].nt)
                    cc3POS_12_neg = zeros(Complex{Float64}, C3data[i].nt)

                    cc3NEG_12_pos = zeros(Complex{Float64}, C3data[i].nt)
                    cc3NEG_12_neg = zeros(Complex{Float64}, C3data[i].nt)

                    if config.IsSpectralNormalization
                        for j = 1:C3data[i].nt
                            
                            cc3POS_12_pos[j] = conj(cc1_v1_pos[j]) * cc1_v2_pos[j] / (abs(cc1_v1_pos[j]) * abs(cc1_v2_pos[j]))
                            cc3POS_12_neg[j] = cc1_v1_pos[j] * conj(cc1_v2_pos[j]) / (abs(cc1_v1_pos[j]) * abs(cc1_v2_pos[j]))

                            cc3NEG_12_pos[j] = conj(cc1_v1_neg[j]) * cc1_v2_neg[j] / (abs(cc1_v1_neg[j]) * abs(cc1_v2_neg[j]))
                            cc3NEG_12_neg[j] = cc1_v1_neg[j] * conj(cc1_v2_neg[j]) / (abs(cc1_v1_neg[j]) * abs(cc1_v2_neg[j]))
                        end
                    else
                        
                        for j = 1:C3data[i].nt
                            
                            cc3POS_12_pos[j] = conj(cc1_v1_pos[j]) * cc1_v2_pos[j]
                            cc3POS_12_neg[j] = cc1_v1_pos[j] * conj(cc1_v2_pos[j])

                            cc3NEG_12_pos[j] = conj(cc1_v1_neg[j]) * cc1_v2_neg[j]
                            cc3NEG_12_neg[j] = cc1_v1_neg[j] * conj(cc1_v2_neg[j])
                        end
                    end

                    cc3_POStemp = real.(vcat(ifft(cc3POS_12_neg)[end:-1:2], ifft(cc3POS_12_pos)))
                    cc3_NEGtemp = real.(vcat(ifft(cc3NEG_12_neg)[end:-1:2], ifft(cc3NEG_12_pos)))
                    cc3_POSandNEG = cc3_POStemp .+ cc3_NEGtemp

                    C3data[i].cc3_POS = cc3_POStemp
                    C3data[i].cc3_NEG = cc3_NEGtemp
                    C3data[i].cc3_ALL = cc3_POSandNEG

                    #-----Save CC3----#
                    name_o = @sprintf("TimeID%02d/Day%04d/UnitID%04d/CC3.V%02d-%02d-%02d", Time_ID, Day_ID, Synthesize_unit_ID, C3data[i].Receiver_v, C3data[i].Receiver_1,   C3data[i].Receiver_2)

                    g1 = g_create(hdf5_o_C3, name_o)

                    g1["CC3"] = C3data[i].cc3_ALL
                    g1["CC1_v1_all"] = cc1_v1_temp
                    g1["CC1_v2_all"] = cc1_v2_temp
                    g1["CC1_v1_windowed"] = cc1_v1_temp_windowed
                    g1["CC1_v2_windowed"] = cc1_v2_temp_windowed

                    attrs(g1)["fs"] = fs
                    attrs(g1)["Receiver_v"] = C3data[i].Receiver_v
                    attrs(g1)["Receiver_1"] = C3data[i].Receiver_1
                    attrs(g1)["Receiver_2"] = C3data[i].Receiver_2
                    attrs(g1)["rvx"] = C3data[i].rv_loc_x
                    attrs(g1)["rvy"] = C3data[i].rv_loc_y
                    attrs(g1)["r1x"] = C3data[i].r1_loc_x
                    attrs(g1)["r1y"] = C3data[i].r1_loc_y
                    attrs(g1)["r2x"] = C3data[i].r2_loc_x
                    attrs(g1)["r2y"] = C3data[i].r2_loc_y


                end
                
                #progress bar
                next!(prog)

            end
        end
    end

    close(hdf5_o_C1)
    close(hdf5_o_C2)
    close(hdf5_o_C3)
    close(hdf5_o_R)

end

#Run the main program
main(problem_name="../EXAMPLE/" * "multi4_source";)


