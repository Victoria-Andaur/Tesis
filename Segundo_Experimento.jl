import Pkg
Pkg.update("JuMP")
Pkg.add("Plots")
Pkg.status("JuMP")
Pkg.add("DataFrames")

using JuMP, Gurobi, Plots
using DataFrames
using Printf, DelimitedFiles;
# Abrir un solo ambiente de Gurobi, Jupyter se marea si se crean demasiados modelos durante en una misma celda.
const GUROBI_ENV = Gurobi.Env()

Tur,headerTur=readdlm("Turbinas" * ".txt",'\t', Float64, '\n',skipstart=0,header=true);
Turbinas=DataFrame(Tur, vec(headerTur));

Rec,headerRec=readdlm("Receptores" * ".txt",'\t', Float64, '\n',skipstart=0,header=true);
Receptores=DataFrame(Rec, vec(headerRec));

Modos_T1 = Dict() 
for i in 1:16
    Modo,header=readdlm("Ruido T1/T1 6800 Modo "*string(i-1)*".txt",'\t', '\n',header=true);
    Modos_T1[i]= DataFrame(Modo, vec(header))
end

Modos_T2 = Dict() 
for i in 1:10
    Modo,header=readdlm("Ruido T2/T2 7200 Modo "*string(i-1)*".txt",'\t', '\n',header=true);
    Modos_T2[i]= DataFrame(Modo, vec(header))
end

Energias_T1 = Dict() 
for i in 1:16
    Energia,header_energia=readdlm("Energía T1/T1 6800 Modo "*string(i-1)*" Energia.txt",'\t', '\n',header=true);
    Energias_T1[i]= DataFrame(Energia, vec(header_energia))
end
Energias_T2 = Dict() 
for i in 1:10
    Energia,header_energia=readdlm("Energía T2/T2 7200 Modo "*string(i-1)*" Energia.txt",'\t', '\n',header=true);
    Energias_T2[i]= DataFrame(Energia, vec(header_energia))
end

N=length(Tur[:,1]); #numero turbinas
J1=length(Modos_T1); #Modos de turbina T1
J2=length(Modos_T2); #Modos de turbinas T2
K=length(Rec[:,1]); #Receptores

vel=readdlm("Velocidad Viento" * ".txt",'\t', Float64, '\n',skipstart=5);
x_val0=0;
y_val0=0;
x_val1=28600;
y_val1=23760;
x_length=length(vel[1,:])
y_length=length(vel[:,1])
densidad=1.192

E1=zeros(Float64, x_length, y_length,J1);
for j in 1:J1
    E1[:,:,j]=readdlm("Mapa Energía T1/MapaEnergia T1 Modo "*string(j-1)*".csv",',',header=false);
end

E2=zeros(Float64, x_length, y_length,J2);
for j in 1:J2
    E2[:,:,j]=readdlm("Mapa Energía T2/MapaEnergia T2 Modo "*string(j-1)*".csv",',',header=false);
end

R1=zeros(Float64, x_length, y_length,J1,K);
for j in 1:J1
    for k in 1:K
        R1[:,:,j,k]=readdlm("Mapa Ruido T1 octavas/MapaRuido_T1_Modo"*string(j-1)*"_receptor "*string(k)*"_octavas.csv", ',',header=false);
    end
end

R2=zeros(Float64, x_length, y_length,J2,K);
for j in 1:J2
    for k in 1:K
        filepath = "Mapa Ruido T2 octavas/MapaRuido_T2_Modo" * string(j-1) * "_receptor " * string(k) * "_octavas.csv"
        if isfile(filepath)
            try
                R2[:,:,j,k] = readdlm(filepath, ',', header=false)
            catch e
                println("Error al leer el archivo $filepath: $e")
            end
        else
            println("El archivo no existe: $filepath")
        end
    end
end


div_x=50;
div_y=50;
x_step=floor(x_length/div_x);
y_step=floor(y_length/div_y);
modelo1 = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)))
@variable(modelo1,0<=x[1:div_x,1:div_y,1:J1]<=1 ,Bin);
@objective(modelo1, Max, sum(sum(sum(E1[round(Int,ix*x_step-x_step+1),round(Int,iy*y_step-y_step+1),j]*x[ix,iy,j] for ix in 1:div_x ) for iy in 1:div_y) for j in 1:J1));
for ix in 1:div_x
    for iy in 1:div_y
        @constraint(modelo1,sum(x[ix,iy,:])<=1)
    end
end
for j in 1:J1
    @constraint(modelo1,sum(x[:,:,j])<=30)
end


for k in 1:K
    @constraint(modelo1,sum(sum(sum(10^((abs(R1[round(Int,ix*x_step-x_step+1),round(Int,iy*y_step-y_step+1),j,k]))/10)*x[ix,iy,j] for ix in 1:div_x) for iy in 1:div_y) for j in 1:J1)<=10^4)
end

@elapsed optimize!(modelo1)

Mx=value.(x);

div_x=50;
div_y=50;
x_step=floor(x_length/div_x);
y_step=floor(y_length/div_y);
modelo2 = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)))
@variable(modelo2,0<=z[1:div_x,1:div_y,1:J2]<=1 ,Bin);
@objective(modelo2, Max, sum(sum(sum(E2[round(Int,ix*x_step-x_step+1),round(Int,iy*y_step-y_step+1),j]*z[ix,iy,j] for ix in 1:div_x ) for iy in 1:div_y) for j in 1:J2));
for ix in 1:div_x
    for iy in 1:div_y
        @constraint(modelo2,sum(z[ix,iy,:])<=1)
    end
end
for j in 1:J2
    @constraint(modelo2,sum(z[:,:,j])<=30)
end


for k in 1:K
    @constraint(modelo2,sum(sum(sum(10^(((R2[round(Int,ix*x_step-x_step+1),round(Int,iy*y_step-y_step+1),j,k]/10)))*z[ix,iy,j] for ix in 1:div_x) for iy in 1:div_y) for j in 1:J2)<=10^4)
end

@elapsed optimize!(modelo2)

Mz=value.(z);


Turbinas_solucion= DataFrame("X [m]" => Float64[], "Y [m]" => Float64[], "NUM" => Int[], "MODO" => Int[],"Index_x"=> Int[],"Index_y"=> Int[])
cont=1;
for ix in 1:div_x
    for iy in 1:div_y
        for k in 1:K
            if Mz[ix,iy,k]!=0
                push!(Turbinas_solucion, ((x_val0+((ix-1)*x_step+1)*(x_val1-x_val0)/x_length), (y_val0+((iy-1)*y_step+1)*(y_val1-y_val0)/y_length), cont,k,(ix-1)*x_step+1,(iy-1)*y_step+1))
                cont+=1
            end
        end
    end
end

function graficar_solucion(Turbinas,Receptores)
    plot(Receptores."X [m]",Receptores."Y [m]",seriestype=:scatter, label="Receptores")
    plot!(Turbinas."X [m]",Turbinas."Y [m]",seriestype=:scatter, label="Turbinas")
    annotate!(Receptores."X [m]".+10,Receptores."Y [m]".-150,text.(Int.(Receptores."Receptor"),:down,9))
    annotate!(Turbinas."X [m]".+10,Turbinas."Y [m]".-150,text.(Int.(Turbinas."NUM"),:down,9))
end

function ruido_receptores(ruido,Recep,Turb)
    for k in 1:length(Recep[:,1])
        print("El ruido que siente el receptor ")
        print(k)
        print(" es ")
        print(round(10*log10(sum(10^(ruido[Turb."Index_x"[i],Turb."Index_y"[i],Turb."MODO"[i],k]/10) for i in 1:length(Turb."Index_x"))),digits=2))
        print("\n")
    end
end


graficar_solucion(Turbinas_solucion,Receptores)
ruido_receptores(R2,Receptores,Turbinas_solucion)

# Función para generar y filtrar variaciones válidas
function valid_variations(x, y,val_sep=1)
    variations = []
    # Generar y filtrar variaciones
    for (new_x, new_y) in [(x - val_sep, y), (x, y), (x + val_sep, y),
                            (x, y - val_sep), (x, y + val_sep),
                            (x - val_sep, y - val_sep), (x + val_sep, y + val_sep),
                            (x - val_sep, y + val_sep), (x + val_sep, y - val_sep)]
        if new_x > 1 && new_x < x_length && new_y > 1 && new_y < y_length
            push!(variations, (new_x, new_y))
        end
    end
    return variations
end



Indices=[]
for aero in eachrow(Turbinas_solucion)
    index_x = aero["Index_x"]
    index_y = aero["Index_y"]
    for (new_x, new_y) in valid_variations(index_x, index_y)
        push!(Indices, (new_x, new_y))
    end
end

modelo_2 = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)))
num_var=length(Indices)
#set_silent(modelo)
@variable(modelo_2,0<=x[1:num_var,1:num_var,1:J]<=1, Bin );
@objective(modelo_2, Max, sum(sum(sum(E1[ix,iy,j]*x[index,index,j] for (index,(ix,iy)) in enumerate(Indices))) for j in 1:J1));
for valor in 1:num_var
    @constraint(modelo_2,sum(x[valor,valor,:])<=1)
end
for j in 1:J1
    @constraint(modelo_2,sum(x[:,:,j])<=30)
end

for k in 1:K
    @constraint(modelo_2,sum(sum(sum(R[ix,iy,j,k]*x[index,index,j] for (index,(ix,iy)) in enumerate(Indices))) for j in 1:J1)<=40)
end

@elapsed optimize!(modelo_2)

Mx=value.(x);

Turbinas_solucion= DataFrame("X [m]" => Float64[], "Y [m]" => Float64[], "NUM" => Int[], "MODO" => Int[],"Index_x"=> Int[],"Index_y"=> Int[])
cont=1;
for (index,(ix,iy)) in enumerate(Indices)
    for k in 1:16
        if Mx[index,index,k]!=0
            push!(Turbinas_solucion, (x_val0+ix*(x_val1-x_val0)/x_length, y_val0+iy*(y_val1-y_val0)/y_length, cont,k,ix,iy))
            cont+=1
        end
    end
end

graficar_solucion(Turbinas_solucion,Receptores)
ruido_receptores(R,Receptores,Turbinas_solucion)
