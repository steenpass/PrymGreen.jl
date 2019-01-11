using PrymGreen
using Random
using Test

@test run_example("./data/PCNC@26431g12_2%1.sd"; print_info = false)
@test string(PrymGreen.random_PCNC(8, 2, Random.MersenneTwister(0))) * "\n" ==
        read(open("./data/random_PCNC", "r"), String)
