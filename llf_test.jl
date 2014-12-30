println("Loading modules")

# code for testing functions
using Images, TestImages, ImageView
require("llf.jl")

img = separate(testimage("mandrill"))
I0 = convert(Array{Float64,3}, img.data)
I0 = impyramid(I0[:,:,1])
I0 = impyramid(I0)
R = I0[:,:,1]

view(R)

R2 = locallaplacianfilter(R, 0.1, levels=3)

view(R2)
