using Plots
include("imtz.jl")
filename="rbpat15.gif"
anim = @animate for i=0:250:150000
    println(i)
    imtz(i)
end
gif(anim,filename,fps=15)
