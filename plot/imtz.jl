using DelimitedFiles
using Printf
using Plots
using Statistics

"""
    imtz(t; clims=(-0.5,0.5), cmap=:pu_or, title="", anom=false, L=2*pi)

Plots a vertical rbsolve section at step `t`. 

# Optional keywords
- `clims::Tuple{Real,Real}`: limits for the colormap
- `cmap::Symbol`: the colormap to use
- `title::String`: title of the plot
- `anom::Bool`: compute horizontal anomaly
- `L::Real`: X domain size
"""
function imtz(t; clims=(-0.5,0.5), cmap=:pu_or, title="", anom=false, L=2*pi)
    coord=readdlm("coord")
    (nz,nc)=size(coord)
    nz=convert(Int64,(nz-1)/2)+1
    coord=coord[1:nz,4]
    nz=nz-2
    filename=@sprintf("t%07d.zzz", t)
    tt=readdlm(filename)
    n=length(tt); n=n÷nz
    zz=coord[2:(nz+1)]
    xx=(0:(n-1))/n*L
    tt=reshape(tt, n, nz)
    if anom
       tt1=tt.-ones(n,1)*mean(tt,dims=1)
    else
       tt1=tt
    end
    plot(xx, zz.+0.5, tt1', st=:heatmap, c=cmap, clims=clims)
    xlabel!("x")
    ylabel!("z")
    title!(title)
end

