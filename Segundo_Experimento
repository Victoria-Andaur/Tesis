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

Modos = Dict() 
for i in 1:16
    Modo,header=readdlm("N175 6.8 Modo "*string(i)*".txt",'\t', '\n',header=true);
    Modos[i]= DataFrame(Modo, vec(header))
end

Energias = Dict() 
for i in 1:16
    Energia,header_energia=readdlm("Energia Modo "*string(i)*".txt",'\t', '\n',header=true);
    Energias[i]= DataFrame(Energia, vec(header_energia))
end

N=length(Tur[:,1]); #numero turbinas
J=length(Modos); #Modos de turbinas
K=length(Rec[:,1]); #Receptores

vel=readdlm("Velocidad Viento" * ".txt",'\t', Float64, '\n',skipstart=5);
x_val0=7.000210E+05;
y_val0=5.850087E+06;
x_val1=7.286210E+05;
y_val1=5.873847E+06;
x_length=length(vel[1,:])
y_length=length(vel[:,1])
densidad=1.192
E=zeros(Float64, x_length, y_length,J);
for j in 1:J
    E[:,:,j]=readdlm("MapaEnergia"*string(j)*".csv",',',header=false);
end

R=zeros(Float64, x_length, y_length,J,K);
for j in 1:J
    for k in 1:K
        R[:,:,j,k]=readdlm( "MapaRuido"*string(j)*"_receptor "*string(k)*".csv", ',',header=false);
    end
end

div_x=50;
div_y=50;
x_step=floor(x_length/div_x);
y_step=floor(y_length/div_y);
modelo = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)))
#set_silent(modelo)
@variable(modelo,0<=x[1:div_x,1:div_y,1:J]<=1 ,Bin);
@objective(modelo, Max, sum(sum(sum(E[round(Int,ix*x_step-x_step+1),round(Int,iy*y_step-y_step+1),j]*x[ix,iy,j] for ix in 1:div_x ) for iy in 1:div_y) for j in 1:J));
for ix in 1:div_x
    for iy in 1:div_y
        @constraint(modelo,sum(x[ix,iy,:])<=1)
    end
end
for j in 1:J
    @constraint(modelo,sum(x[:,:,j])<=30)
end

for k in 1:K
    @constraint(modelo,sum(sum(sum(R[round(Int,ix*x_step-x_step+1),round(Int,iy*y_step-y_step+1),j,k]*x[ix,iy,j] for ix in 1:div_x) for iy in 1:div_y) for j in 1:J)<=40)
end

@elapsed optimize!(modelo)

Mx=value.(x);

Turbinas_solucion= DataFrame("X [m]" => Float64[], "Y [m]" => Float64[], "NUM" => Int[], "MODO" => Int[],"Index_x"=> Int[],"Index_y"=> Int[])
cont=1;
for ix in 1:div_x
    for iy in 1:div_y
        for k in 1:16
            if Mx[ix,iy,k]!=0
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
        print(round(sum(ruido[Turb."Index_x"[i],Turb."Index_y"[i],Turb."MODO"[i],k] for i in 1:length(Turb."Index_x")),digits=2))
        print("\n")
    end
end


graficar_solucion(Turbinas_solucion,Receptores)
ruido_receptores(R,Receptores,Turbinas_solucion)

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
@objective(modelo_2, Max, sum(sum(sum(E[ix,iy,j]*x[index,index,j] for (index,(ix,iy)) in enumerate(Indices))) for j in 1:J));
for valor in 1:num_var
    @constraint(modelo_2,sum(x[valor,valor,:])<=1)
end
for j in 1:J
    @constraint(modelo_2,sum(x[:,:,j])<=30)
end

for k in 1:K
    @constraint(modelo_2,sum(sum(sum(R[ix,iy,j,k]*x[index,index,j] for (index,(ix,iy)) in enumerate(Indices))) for j in 1:J)<=40)
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
