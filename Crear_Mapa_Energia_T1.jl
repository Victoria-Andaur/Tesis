using DataFrames
using Printf, DelimitedFiles;

Energias_T1 = Dict() 
for i in 1:16
    Energia,header_energia=readdlm("Energía T1/T1 6800 Modo "*string(i-1)*" Energia.txt",'\t', '\n',header=true);
    Energias_T1[i]= DataFrame(Energia, vec(header_energia))
end
vel=readdlm("Velocidad Viento" * ".txt",'\t', Float64, '\n',skipstart=5);
dens_obj =1.192
filas,cols=size(vel)

densidades = parse.(Float64, names(Energias_T1[1])[2:end])
_, idy = findmin(abs.(densidades .- dens_obj))
densidad_cercana = densidades[idy]
Energias_mapa_modos_T1=zeros(filas,cols,16)
for m in 1:16
    for i in 1:filas
        for j in 1:cols
            _, idx = findmin(abs.(Energias_T1[m]."Air Density [kg/m3]" .- vel[i,j]))
            Energias_mapa_modos_T1[i,j,m]=Energias_T1[m][!, Symbol(densidad_cercana)][idx]
        end
    end
end
for j in 1:16
    writedlm("Mapa Energía T1/MapaEnergia T1 Modo "*string(j-1)*".csv",Energias_mapa_modos_T1[:,:,j],',')
end
