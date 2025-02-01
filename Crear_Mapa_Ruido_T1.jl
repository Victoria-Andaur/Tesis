using DataFrames
using Printf, DelimitedFiles, CSV;



vel=readdlm("Velocidad Viento" * ".txt",'\t', Float64, '\n',skipstart=5);
Rec,headerRec=readdlm("Receptores" * ".txt",'\t', Float64, '\n',skipstart=0,header=true);
Receptores=DataFrame(Rec, vec(headerRec));

dens_obj =1.192

Alturas_data=readdlm("Curvas de nivel 50m" * ".asc",' ',skipblanks=true, '\n',skipstart=6,header=false);
Alturas=convert(Matrix{Float64},Alturas_data[1:609,1:663]);

Energias_T1 = Dict() 
for i in 1:16
    Energia,header_energia=readdlm("Energía T1/T1 6800 Modo "*string(i-1)*" Energia.txt",'\t', '\n',header=true);
    Energias_T1[i]= DataFrame(Energia, vec(header_energia))
end
Modos_T1 = Dict() 
for i in 1:16
    Modo,header=readdlm("Ruido T1/T1 6800 Modo "*string(i-1)*".txt",'\t', '\n',header=true);
    Modos_T1[i]= DataFrame(Modo, vec(header))
end


J1=length(Modos_T1); #Modos de turbinas
K=length(Rec[:,1]); #Receptores

filas,cols=size(vel)
densidad=1.192

function obtener_altura(x,y,Alturas,ncols=663,nrows=609,xllcorner=-1232.0,yllcorner=-810.0,paso_z=50)
    return Alturas[round(Int,(x-xllcorner)/paso_z),round(Int,(y-yllcorner)/paso_z)]
end

function func_velocidad_1(x,y,mapa_vel=vel)
    return vel[x,y]
end

function distancia_turbina_1(x,y,k,filas=filas,cols=cols,x_val1=28600,y_val1=23760)
    pasox=x_val1/cols
    pasoy=y_val1/filas
    x2=pasox*(x-1)
    y2=pasoy*(y-1)
    z2=obtener_altura(x2,y2,Alturas)
    x1=Receptores."X [m]"[k]
    y1=Receptores."Y [m]"[k]
    z1=Receptores."Altura[m]"[k]
    dx = (x2 - x1)^2
    dy = (y2 - y1)^2
    dz = (143+z2 - z1)^2
    return sqrt(dx+dy+dz)
end

function func_ruido_turbina_octavas_T1(j)
    indice=findall(z->z==Int(round(12)), Modos_T1[j]."V[m/s]")
    return [Modos_T1[j]."31Hz"[indice[1]],Modos_T1[j]."63Hz"[indice[1]],Modos_T1[j]."125Hz"[indice[1]],Modos_T1[j]."250Hz"[indice[1]],Modos_T1[j]."500Hz"[indice[1]],Modos_T1[j]."1kHz"[indice[1]],Modos_T1[j]."2kHz"[indice[1]],Modos_T1[j]."4kHz"[indice[1]],Modos_T1[j]."8kHz"[indice[1]]]
end
function suma_ruido(vector_ruido)
    suma_parcial=0
    for db in vector_ruido
        suma_parcial+=10^(db/10)
    end
    if suma_parcial>0
        return 10*log10(suma_parcial)
    else
        return 0
    end
end

function atenuacion_atmosferica(distancia)
    #Distancia en m
    ate=[0,0.1,0.4,1.0,1.9,3.7,9.7,32.8,117]
    return ate.*distancia/1000
end

function ruido_octavas_T1(x,y,j,k)
    dist=distancia_turbina_1(x,y,k)
    if dist!=0
        #ate_div=(20*log10(dist)+11)
        ate_div=(10*log10(4*pi*dist^2))
        ate_atm=func_ruido_turbina_octavas_T1(j)-atenuacion_atmosferica(dist)
        return max(suma_ruido(ate_atm)-ate_div,0)
    else
        return func_ruido_turbina_1_T1(x,y,j)
    end
end

R=zeros(Float64,filas,cols, J1, K)
for x in 1:filas
    for y in 1:cols
        for j in 1:J1
            for k in 1:K
                R[x,y,j,k]=ruido_octavas_T1(x,y,j,k)
            end
        end
    end
end


for j in 1:J1
    for k in 1:K
        CSV.write("Mapa Ruido T1 octavas/MapaRuido_T1_Modo"*string(j-1)*"_receptor "*string(k)*"_octavas.csv", DataFrame(R[:,:,j,k],:auto))
    end
end


