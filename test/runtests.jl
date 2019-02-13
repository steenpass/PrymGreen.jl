using PrymGreen
using Random
using Test

# writes a file 'pcnc_g12_submatrix_399'
@test run_example("./data/PCNC@26431g12_2%1.sd"; print_info = true)

filename = "pcnc_g12_submatrix_399"
checksum = "3b823f6d15e78a31f6e8cc047de45c07645c85edaf50203c4485e0de9d8f8bed  "
checksum *= filename * "\n"
@test read(`sha256sum $filename`, String) == checksum
rm(filename)

s1 = string(PrymGreen.random_PCNC(8, 2, Random.MersenneTwister(0))) * "\n"
s2 = read(open("./data/random_PCNC", "r"), String)
@test s1 == s2
