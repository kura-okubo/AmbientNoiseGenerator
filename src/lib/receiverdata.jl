module ReceiverData

export sreceiverdata, sC1data, sC2data

abstract type Abstractreceiverdata end
abstract type AbstractC1 end
abstract type AbstractC2 end
abstract type AbstractC3 end

mutable struct sreceiverdata <: Abstractreceiverdata

    #define contents
    Receiver_ID ::Int64
    loc_x ::Float64
    loc_y ::Float64
    t :: Array{Float64}
    nt :: Int64
    u ::  Array{Float64}

    sreceiverdata() = new()
end


mutable struct sC1data <: AbstractC1

    #define contents
    Receiver_1 ::Int64
    Receiver_2 ::Int64
    r1_loc_x ::Float64
    r1_loc_y ::Float64
    r2_loc_x ::Float64
    r2_loc_y ::Float64
    dist     ::Float64
    lag_t :: Array{Float64}
    nt :: Int64
    cc1 ::  Array{Float64}

    sC1data() = new()
end

mutable struct sC2data <: AbstractC2

    #define contents
    Receiver_v ::Int64
    Receiver_1 ::Int64
    Receiver_2 ::Int64
    rv_loc_x ::Float64
    rv_loc_y ::Float64
    r1_loc_x ::Float64
    r1_loc_y ::Float64
    r2_loc_x ::Float64
    r2_loc_y ::Float64
    dist_v_1     ::Float64
    dist_v_2     ::Float64
    dist_1_2     ::Float64
    lag_t :: Array{Float64}
    nt :: Int64

    cc2_POS ::  Array{Float64}
    cc2_NEG ::  Array{Float64}
    cc2_ALL ::  Array{Float64}

    sC2data() = new()
end


mutable struct sC3data <: AbstractC3

    #define contents
    Receiver_v ::Int64
    Receiver_1 ::Int64
    Receiver_2 ::Int64
    rv_loc_x ::Float64
    rv_loc_y ::Float64
    r1_loc_x ::Float64
    r1_loc_y ::Float64
    r2_loc_x ::Float64
    r2_loc_y ::Float64
    dist_v_1     ::Float64
    dist_v_2     ::Float64
    dist_1_2     ::Float64
    lag_t :: Array{Float64}
    nt :: Int64

    CC1_v1_all ::  Array{Float64}
    CC1_v2_all ::  Array{Float64}
    CC1_v1_windowed ::  Array{Float64}
    CC1_v2_windowed ::  Array{Float64}

    cc3_POS ::  Array{Float64}
    cc3_NEG ::  Array{Float64}
    cc3_ALL ::  Array{Float64}

    sC3data() = new()
end

end
