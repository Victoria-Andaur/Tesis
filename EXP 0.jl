using JuMP, Gurobi, Plots
using DataFrames
using Printf, DelimitedFiles;
using CSV
# Abrir un solo ambiente de Gurobi, Jupyter se marea si se crean demasiados modelos durante en una misma celda.
const GUROBI_ENV = Gurobi.Env()

Tur,headerTur=readdlm("Turbinas" * ".txt",'\t', Float64, '\n',skipstart=0,header=true);
Turbinas=DataFrame(Tur, vec(headerTur));

Rec,headerRec=readdlm("Receptores" * ".txt",'\t', Float64, '\n',skipstart=0,header=true);
Receptores=DataFrame(Rec, vec(headerRec));

Poligono=readdlm("Poligono" * ".txt",'\t', '\n',skipstart=0,header=false);
contador_prueba=1;

Poligono_para_grafico=DataFrame("X [m]" => Float64[], "Y [m]" => Float64[]);
for val in Poligono
    if contador_prueba%4==0
        push!(Poligono_para_grafico,(parse(Float64,Poligono[contador_prueba][1:findfirst(',',Poligono[contador_prueba])-1]),parse(Float64,Poligono[contador_prueba][findfirst(',',Poligono[contador_prueba])+1:length(Poligono[contador_prueba])])));
    end
    contador_prueba+=1
end

N=length(Tur[:,1]); #numero turbinas
J1=length(Modos_T1); #Modos de turbinas T2
J2=length(Modos_T2); #Modos de turbinas T2
K=length(Rec[:,1]); #Receptores

vel=readdlm("Velocidad Viento" * ".txt",'\t', Float64, '\n',skipstart=5);
x_val0=0;
y_val0=0;
x_val1=28600;
y_val1=23760;
filas,cols=size(vel)
densidad=1.192

E1=zeros(Float64, filas,cols,J1);
for j in 1:J1
    E1[:,:,j]=readdlm("Mapa Energía T1/MapaEnergia T1 Modo "*string(j-1)*".csv",',',header=false);
end

R1=zeros(Float64, filas, cols,J1,K);
for j in 1:J1
    for k in 1:K
        filepath = "Mapa Ruido T1 octavas/MapaRuido_T1_Modo" * string(j-1) * "_receptor " * string(k) * "_octavas.csv"
        if isfile(filepath)
            try
                R1[:,:,j,k] = readdlm(filepath, ',', header=false,skipstart=1)
            catch e
                println("Error al leer el archivo $filepath: $e")
            end
        else
            println("El archivo no existe: $filepath")
        end
    end
end


E2=zeros(Float64, filas,cols,J2);
for j in 1:J2
    E2[:,:,j]=readdlm("Mapa Energía T2/MapaEnergia T2 Modo "*string(j-1)*".csv",',',header=false);
end

R2=zeros(Float64, filas, cols,J2,K);
for j in 1:J2
    for k in 1:K
        filepath = "Mapa Ruido T2 octavas/MapaRuido_T2_Modo" * string(j-1) * "_receptor " * string(k) * "_octavas.csv"
        if isfile(filepath)
            try
                R2[:,:,j,k] = readdlm(filepath, ',', header=false,skipstart=1)
            catch e
                println("Error al leer el archivo $filepath: $e")
            end
        else
            println("El archivo no existe: $filepath")
        end
    end
end

function graficar_solucion(Turbinas,Receptores,Nombre,J)
    paleta = palette([:green, :gray], J)
    plot(Poligono_para_grafico."X [m]",Poligono_para_grafico."Y [m]",seriestype=:scatter,label="Poligono",seriesalpha=0.05)
    plot!(Receptores."X [m]",Receptores."Y [m]",seriestype=:scatter, label="Receptores",markershape = :star,markercolor = :red)
    x_dummy = [-100, -100]  # Coordenadas X ficticias
    y_dummy = [-100, -100]  # Coordenadas Y ficticias
    modo_dummy = [0, J-1]  # Valores mínimos y máximos de MODO
    plot!([Turbinas."X [m]";x_dummy],[Turbinas."Y [m]";y_dummy],seriestype=:scatter,marker_z = [Turbinas."MODO";modo_dummy],colorbar_title = "Modo de ruido",  color = cgrad(paleta, J, categorical=true), label="Turbinas")
    xlims!(minimum(Poligono_para_grafico."X [m]")-300, maximum(Poligono_para_grafico."X [m]")+300)  # Cambiar los límites del eje X
    ylims!(minimum(Poligono_para_grafico."Y [m]")-300, maximum(Poligono_para_grafico."Y [m]")+300)   # Cambiar los límites del eje Y
    annotate!(Receptores."X [m]".+10,Receptores."Y [m]".-200,text.(Int.(Receptores."Receptor"),:down,9))
    annotate!(Turbinas."X [m]".+10,Turbinas."Y [m]".-200,text.(Int.(Turbinas."NUM"),:down,9))
    savefig(Nombre)
end


modelo_dummy_T1=Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)))

@variable(modelo_dummy_T1, x_T1[1:N,1:J1], Bin)
@objective(modelo_dummy_T1, Max, sum(x_T1[i,j]*E1[Int(round(Turbinas."Y [m]"[i]*filas/y_val1)),Int(round(Turbinas."X [m]"[i]*cols/x_val1)),j] for i in 1:N, j in 1:J1))
@constraint(modelo_dummy_T1, sum(sum(x_T1[i,j] for i in 1:N) for j in 1:J1) <= 30)
for i in 1:N
    @constraint(modelo_dummy_T1, sum(x_T1[i,j] for j in 1:J1)<=1)
end

suma_ruido_T1=0
for k in 1:K
    for i in 1:N
        for j in 1:J1
            suma_ruido_T1+=10^(R1[Int(round(Turbinas."Y [m]"[i]*filas/y_val1)),Int(round(Turbinas."X [m]"[i]*cols/x_val1)),j,k]/10)*x_T1[i,j] 
        end
    end
    @constraint(modelo_dummy_T1,suma_ruido_T1<=10^4)
end


elapsed_time_T1 = @elapsed optimize!(modelo_dummy_T1)


Var_T1=value.(x_T1);
Energia_T1 = objective_value(modelo_dummy_T1)

Turbinas_solucion_dummy_T1 = DataFrame("X [m]" => Float64[], "Y [m]" => Float64[], "NUM" => Int[], "MODO" => Int[], "Index_x" => Int[], "Index_y" => Int[])
cont = 1
for i in 1:N
    for j in 1:J1
        if Var_T1[i, j]!= 0 
            push!(Turbinas_solucion_dummy_T1, (Turbinas."X [m]"[i],Turbinas."Y [m]"[i],cont, j-1, Int(round(Turbinas."X [m]"[i]*cols/x_val1)), Int(round(Turbinas."Y [m]"[i]*filas/y_val1))))
            cont += 1      
        end
    end    
end

Ruido_receptores_T1=zeros(Float64,K);
for k in 1:K
    for i in 1:N
        for j in 1:J1
            Ruido_receptores_T1[k]+=10^(R1[Int(round(Turbinas."Y [m]"[i]*filas/y_val1)),Int(round(Turbinas."X [m]"[i]*cols/x_val1)),j,k]/10)*Var_T1[i,j] 
        end
    end
end
for k in 1:K
    Ruido_receptores_T1[k]=10*log10(Ruido_receptores_T1[k])
end
graficar_solucion(Turbinas_solucion_dummy_T1,Receptores,"Imágenes/0. Exp0 Turbinas/Mapa Turbinas T1.svg",J1)
CSV.write("Experimento 0/T1.csv", Turbinas_solucion_dummy_T1)
writedlm("Experimento 0/Ruido_T1.csv", Ruido_receptores_T1)





modelo_dummy_T2=Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)))

@variable(modelo_dummy_T2, x_T2[1:N,1:J2], Bin)
@objective(modelo_dummy_T2, Max, sum(x_T2[i,j]*E2[Int(round(Turbinas."Y [m]"[i]*filas/y_val1)),Int(round(Turbinas."X [m]"[i]*cols/x_val1)),j] for i in 1:N, j in 1:J2))
@constraint(modelo_dummy_T2, sum(sum(x_T2[i,j] for i in 1:N) for j in 1:J2) <= 30)
for i in 1:N
    @constraint(modelo_dummy_T2, sum(x_T2[i,j] for j in 1:J2)<=1)
end

suma_ruido_T2=0
for k in 1:K
    for i in 1:N
        for j in 1:J2
            suma_ruido_T2+=10^(R2[Int(round(Turbinas."Y [m]"[i]*filas/y_val1)),Int(round(Turbinas."X [m]"[i]*cols/x_val1)),j,k]/10)*x_T2[i,j] 
        end
    end
    @constraint(modelo_dummy_T2,suma_ruido_T2<=10^4)
end


elapsed_time_T2=@elapsed optimize!(modelo_dummy_T2)


Var_T2=value.(x_T2);
Energia_T2 = objective_value(modelo_dummy_T2)

Turbinas_solucion_dummy_T2 = DataFrame("X [m]" => Float64[], "Y [m]" => Float64[], "NUM" => Int[], "MODO" => Int[], "Index_x" => Int[], "Index_y" => Int[])
cont = 1
for i in 1:N
    for j in 1:J2
        if Var_T2[i, j]!= 0 
            push!(Turbinas_solucion_dummy_T2, (Turbinas."X [m]"[i],Turbinas."Y [m]"[i],cont, j-1, Int(round(Turbinas."X [m]"[i]*cols/x_val1)), Int(round(Turbinas."Y [m]"[i]*filas/y_val1))))
            cont += 1      
        end
    end    
end
Ruido_receptores_T2=zeros(Float64,K);
for k in 1:K
    for i in 1:N
        for j in 1:J2
            Ruido_receptores_T2[k]+=10^(R2[Int(round(Turbinas."Y [m]"[i]*filas/y_val1)),Int(round(Turbinas."X [m]"[i]*cols/x_val1)),j,k]/10)*Var_T2[i,j] 
        end
    end
end
for k in 1:K
    Ruido_receptores_T2[k]=10*log10(Ruido_receptores_T2[k])
end

graficar_solucion(Turbinas_solucion_dummy_T2,Receptores,"Imágenes/0. Exp0 Turbinas/Mapa Turbinas T2.svg",J2)
CSV.write("Experimento 0/T2.csv", Turbinas_solucion_dummy_T2)

writedlm("Experimento 0/Ruido_T2.csv", Ruido_receptores_T2)


