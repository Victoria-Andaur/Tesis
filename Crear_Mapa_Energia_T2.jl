using DataFrames
using Printf, DelimitedFiles;

Energias_T2 = Dict() 
for i in 1:10
    Energia,header_energia=readdlm("Energía T2/T2 7200 Modo "*string(i-1)*" Energia.txt",'\t', '\n',header=true);
    Energias_T2[i]= DataFrame(Energia, vec(header_energia))
end
vel=readdlm("Velocidad Viento" * ".txt",'\t', Float64, '\n',skipstart=5);
dens_obj =1.192
filas,cols=size(vel)

densidades = parse.(Float64, names(Energias_T2[1])[2:end])
_, idy = findmin(abs.(densidades .- dens_obj))
densidad_cercana = densidades[idy]
Energias_mapa_modos_T2=zeros(filas,cols,10)
for m in 1:10
    for i in 1:filas
        for j in 1:cols
            _, idx = findmin(abs.(Energias_T2[m]."Air Density [kg/m3]" .- vel[i,j]))
            Energias_mapa_modos_T2[i,j,m]=Energias_T2[m][!, Symbol(densidad_cercana)][idx]
        end
    end
end
for j in 1:10
    writedlm("Mapa Energía T2/MapaEnergia T2 Modo "*string(j-1)*".csv",Energias_mapa_modos_T2[:,:,j],',')
end
