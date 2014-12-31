println("Loading modules")

# code for testing functions
using Images, TestImages
require("llf.jl")

img = separate(testimage("mandrill"))
I0 = convert(Array{Float64,3}, img.data)
R = I0[:,:,1]

#R2 = locallaplacianfilter(R, sigma_r = 0.1, levels=3)

imwrite(R, "input.png")
#imwrite(R2, "output.png")

#R3 = fastlocallaplacianfilter(R, sigma_r = 0.1, levels=3)
#imwrite(R3, "output_fast2.png")

img = separate(imread("/Users/loganw/Documents/Julia/local_laplacian_filtering/madison.png"))
I0 = convert(Array{Float64, 3}, img.data)
G = I0[:,:,2]

imwrite(G, "madison_i.png")
G2 = fastlocallaplacianfilter(G, sigma_r = 0.1, alpha=3, levels=5)
imwrite(G2, "madison_o.png")
