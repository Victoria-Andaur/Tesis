using JuMP, Gurobi, Plots
using DataFrames
using Printf, DelimitedFiles;
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

Modos_T2 = Dict() 
for i in 1:10
    Modo,header=readdlm("Ruido T2/T2 7200 Modo "*string(i-1)*".txt",'\t', '\n',header=true);
    Modos_T2[i]= DataFrame(Modo, vec(header))
end

Energias_T2 = Dict() 
for i in 1:10
    Energia,header_energia=readdlm("Energía T2/T2 7200 Modo "*string(i-1)*" Energia.txt",'\t', '\n',header=true);
    Energias_T2[i]= DataFrame(Energia, vec(header_energia))
end

N=length(Tur[:,1]); #numero turbinas
J2=length(Modos_T2); #Modos de turbinas T2
K=length(Rec[:,1]); #Receptores

vel=readdlm("Velocidad Viento" * ".txt",'\t', Float64, '\n',skipstart=5);
x_val0=0;
y_val0=0;
x_val1=28600;
y_val1=23760;
filas,cols=size(vel);
densidad=1.192


E2=zeros(Float64, filas, cols,J2);
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


"""
    funcion_objetivo(J2, divisions_y, divisions_x, E2, x_step, y_step, binary_decision_vars, modelo2)
Función que define la función objetivo del modelo de optimización de la ubicación de las turbinas T2.
"""
function funcion_objetivo(J2, divisions_y, divisions_x, E2, x_step, y_step, binary_decision_vars, modelo2)
    objective_sum = 0
    for j in 1:J2
        for iy in 1:divisions_y
            for ix in 1:divisions_x
                objective_sum += E2[round(Int, iy * y_step - y_step + 1),round(Int, ix * x_step - x_step + 1), j] * binary_decision_vars[ix, iy, j];
            end
        end
    end
    @objective(modelo2, Max, objective_sum);
end
"""
    Un_solo_modo_por_turbina(divisions_y, divisions_x, binary_decision_vars, modelo2)

Función que agrega la restricción de que solo se puede seleccionar un modo por turbina.
"""
function Agregar_restriccion_un_solo_modo_por_turbina(divisions_y, divisions_x, binary_decision_vars, modelo2)
    for ix in 1:divisions_x
        for iy in 1:divisions_y
            @constraint(modelo2,sum(binary_decision_vars[ix,iy,:])<=1);
        end
    end
end
"""
    A_lo_más_n_turbinas(binary_decision_vars, modelo2,N)

Función que agrega la restricción de que se pueden seleccionar a lo más N turbinas.
"""
function Agregar_restriccion_a_lo_más_n_turbinas(binary_decision_vars, modelo2,N)
        @constraint(modelo2,sum(binary_decision_vars[:,:,:])<=N)
end
"""
    Restriccion_ruido(R2, divisions_y, divisions_x, K, binary_decision_vars, modelo2)

Función que agrega la restricción de que el ruido no puede superar los 40 dB.
"""
function Agregar_restriccion_ruido(R2, divisions_y, divisions_x, K, binary_decision_vars, modelo2)
    for k in 1:K
        @constraint(modelo2,(sum(sum(sum(10^(R2[round(Int,iy*y_step-y_step+1),round(Int,ix*x_step-x_step+1),j,k]/10)*binary_decision_vars[ix,iy,j] for ix in 1:divisions_x) for iy in 1:divisions_y) for j in 1:J2))<=10^4)
    end
end
function en_poligono(idx,idy,poligono)
    return "NAME=R"*string(idy)*"R_C"*string(idx)*"C" in poligono
end

function graficar_solucion(Turbinas,Receptores,Nombre,J2)
    paleta = palette([:green, :gray], J2)
    plot(Poligono_para_grafico."X [m]",Poligono_para_grafico."Y [m]",seriestype=:scatter,label="Poligono",seriesalpha=0.05)
    plot!(Receptores."X [m]",Receptores."Y [m]",seriestype=:scatter, label="Receptores",markershape = :star,markercolor = :red)
    x_dummy = [-100, -100]  # Coordenadas X ficticias
    y_dummy = [-100, -100]  # Coordenadas Y ficticias
    modo_dummy = [0, J2-1]  # Valores mínimos y máximos de MODO
    #plot!(Turbinas."X [m]",Turbinas."Y [m]",seriestype=:scatter, label="Turbinas")
    plot!([Turbinas."X [m]";x_dummy],[Turbinas."Y [m]";y_dummy],seriestype=:scatter,marker_z = [Turbinas."MODO";modo_dummy],colorbar_title = "Modo de ruido",  color = cgrad(paleta, 10, categorical=true), label="Turbinas")
    xlims!(minimum(Poligono_para_grafico."X [m]")-300, maximum(Poligono_para_grafico."X [m]")+300)  # Cambiar los límites del eje X
    ylims!(minimum(Poligono_para_grafico."Y [m]")-300, maximum(Poligono_para_grafico."Y [m]")+300)   # Cambiar los límites del eje Y
    annotate!(Receptores."X [m]".+10,Receptores."Y [m]".-200,text.(Int.(Receptores."Receptor"),:down,9))
    annotate!(Turbinas."X [m]".+10,Turbinas."Y [m]".-200,text.(Int.(Turbinas."NUM"),:down,9))
    savefig(Nombre)
end

function mostrar_ruido_receptores(ruido,Recep,Turb)
    for k in 1:length(Recep[:,1])
        print("El ruido que siente el receptor ")
        print(k)
        print(" es ")
        print(round(10*log10(sum(10^(ruido[Turb."Index_y"[i],Turb."Index_x"[i],Turb."MODO"[i],k]/10) for i in 1:length(Turb."Index_x"))),digits=2))
        print("\n")
    end
end

function guardar_resultados(filename::String, N::Int, Energia::Float64, tiempo::Float64,ruido,Recep,Turb)
    ruido_lista=[]
    for k in 1:length(Recep[:,1])
        push!(ruido_lista, round(10*log10(sum(10^(ruido[Turb."Index_y"[i],Turb."Index_x"[i],Turb."MODO"[i]+1,k]/10) for i in 1:length(Turb."Index_x"))), digits=2))
    end
    
    content = """
    \\begin{table}[ht]
    \\begin{tabular}{
    >{\\columncolor[HTML]{808080}}l ll}
    {\\color[HTML]{FFFFFF} \\textbf{Método   utlizado}} & \\multicolumn{1}{r}{\\cellcolor[HTML]{F2F2F2}{\\color[HTML]{595959} Grilla 50x50}} & \\multicolumn{1}{r}{\\cellcolor[HTML]{F2F2F2}{\\color[HTML]{595959} }} \\\\
    {\\color[HTML]{FFFFFF} \\textbf{Número de turbinas máximas}} & {\\color[HTML]{595959} $N} & {\\color[HTML]{595959} } \\\\
    {\\color[HTML]{FFFFFF} \\textbf{Número de turbinas agegadas}} & \\cellcolor[HTML]{F2F2F2}{\\color[HTML]{595959} $(size(Turb)[1])} & \\cellcolor[HTML]{F2F2F2}{\\color[HTML]{595959} } \\\\
    {\\color[HTML]{FFFFFF} \\textbf{Energía producida}} & {\\color[HTML]{595959} $Energia} & {\\color[HTML]{595959} MWh} \\\\
    {\\color[HTML]{FFFFFF} \\textbf{Tiempo de cálculo}} & \\cellcolor[HTML]{F2F2F2}{\\color[HTML]{595959} $tiempo} & \\cellcolor[HTML]{F2F2F2}{\\color[HTML]{595959} s}\\\\
    \\hline\\hline
    \\end{tabular}
    \\centering
    \\end{table}
    
    \\begin{table}[ht]
    \\begin{tabular}{cccc}
    \\rowcolor[HTML]{808080} 
    \\multicolumn{1}{c}{\\cellcolor[HTML]{808080}{\\color[HTML]{FFFFFF} \\textbf{Nº Receptor}}} & {\\color[HTML]{FFFFFF} \\textbf{Pos\\_X}} & {\\color[HTML]{FFFFFF} \\textbf{Pos\\_Y}} & {\\color[HTML]{FFFFFF} \\textbf{dB}} \\\\
    \\rowcolor[HTML]{F2F2F2} 
    {\\color[HTML]{595959} 1} & {\\color[HTML]{595959} 9005} & {\\color[HTML]{595959} 11533} & {\\color[HTML]{595959} $(ruido_lista[1])} \\\\
    {\\color[HTML]{595959} 2} & {\\color[HTML]{595959} 10026} & {\\color[HTML]{595959} 10354} & {\\color[HTML]{595959} $(ruido_lista[2])} \\\\
    \\rowcolor[HTML]{F2F2F2} 
    {\\color[HTML]{595959} 3} & {\\color[HTML]{595959} 9407} & {\\color[HTML]{595959} 12307} & {\\color[HTML]{595959} $(ruido_lista[3])} \\\\
    {\\color[HTML]{595959} 4} & {\\color[HTML]{595959} 12743} & {\\color[HTML]{595959} 12591} & {\\color[HTML]{595959} $(ruido_lista[4])} \\\\
    \\rowcolor[HTML]{F2F2F2} 
    {\\color[HTML]{595959} 5} & {\\color[HTML]{595959} 16483} & {\\color[HTML]{595959} 11370} & {\\color[HTML]{595959} $(ruido_lista[5])} \\\\
    {\\color[HTML]{595959} 6} & {\\color[HTML]{595959} 19317} & {\\color[HTML]{595959} 10737} & {\\color[HTML]{595959} $(ruido_lista[6])} \\\\
    \\rowcolor[HTML]{F2F2F2} 
    {\\color[HTML]{595959} 7} & {\\color[HTML]{595959} 18619} & {\\color[HTML]{595959} 8284} & {\\color[HTML]{595959} $(ruido_lista[7])}\\\\
    \\hline\\hline
    \\end{tabular}
    \\centering
    \\end{table}
        
    \\begin{table}[H]
    \\begin{tabular}{cccc}
    \\rowcolor[HTML]{808080} 
    \\multicolumn{1}{c}{\\cellcolor[HTML]{808080}{\\color[HTML]{FFFFFF} \\textbf{Nº Turbina}}} & 
    \\multicolumn{1}{c}{\\cellcolor[HTML]{808080}{\\color[HTML]{FFFFFF} \\textbf{Pos\\_X}}} & 
    \\multicolumn{1}{c}{\\cellcolor[HTML]{808080}{\\color[HTML]{FFFFFF} \\textbf{Pos\\_Y}}} & 
    {\\color[HTML]{FFFFFF} \\textbf{Modo}}\\\\
    """
        # Save the content to the specified file
    open(filename, "w") do file
        write(file, content)
    

        # Escribir filas de la tabla
        for i in 1:nrow(Turb)
            if i % 2 == 0
                write(file, """\\rowcolor[HTML]{F2F2F2} """)
            end
            write(file, """{\\color[HTML]{595959} $(i)} & {\\color[HTML]{595959} $(round(Turb."X [m]"[i],digits=1))} & {\\color[HTML]{595959} $(round(Turb."Y [m]"[i],digits=1))} & {\\color[HTML]{595959} $(Turb."MODO"[i])} \\\\ """)
        end

        # Escribir cierre de la tabla
        write(file, """
        \\hline\\hline
        \\end{tabular}
        \\centering
        \\end{table}
        """)
    end
    println("Table saved to $filename")
end



function resolver_y_guardar(divisions_x::Int, divisions_y::Int,N::Int, filas, cols, x_val1, y_val1, J2::Int, K::Int, E2::Array{Float64}, R2::Array{Float64}, 
    Receptores::DataFrame)

    # Paso 1: Calcular pasos en X e Y
    x_step=floor(cols/divisions_x);
    y_step=floor(filas/divisions_y);

    # Paso 2: Crear modelo de optimización
    modelo2 = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)))
    @variable(modelo2, 0 <= binary_decision_vars[1:divisions_x, 1:divisions_y, 1:J2] <= 1, Bin)

    # Paso 3: Agregar función objetivo y restricciones
    funcion_objetivo(J2, divisions_y, divisions_x, E2, x_step, y_step, binary_decision_vars, modelo2)
    Agregar_restriccion_un_solo_modo_por_turbina(divisions_y, divisions_x, binary_decision_vars, modelo2)
    Agregar_restriccion_a_lo_más_n_turbinas(binary_decision_vars, modelo2, N)
    Agregar_restriccion_ruido(R2, divisions_y, divisions_x, K, binary_decision_vars, modelo2)
    for ix in 1:divisions_x
        for iy in 1:divisions_y
            @constraint(modelo2,sum(binary_decision_vars[ix,iy,:])<=Int(en_poligono(round(Int,ix*x_step-x_step+1),round(Int,filas-(iy*y_step-y_step+1)),Poligono)))
        end
    end
    
    # Paso 4: Resolver modelo
    elapsed_time = @elapsed optimize!(modelo2)
    Energia = objective_value(modelo2)
    Mz = value.(binary_decision_vars)

    #   Paso 5: Extraer solución
    Turbinas_solucion = DataFrame("X [m]" => Float64[], "Y [m]" => Float64[], "NUM" => Int[], "MODO" => Int[], "Index_x" => Int[], "Index_y" => Int[])
    cont = 1
    for ix in 1:divisions_x
        for iy in 1:divisions_y
            for j in 1:J2
                if Mz[ix, iy, j] != 0
                    push!(Turbinas_solucion, ((((ix - 1) * x_step + 1) * (x_val1) / cols),( ((iy - 1) * y_step + 1) * (y_val1) / filas),cont, j-1, (ix - 1) * x_step + 1, (iy - 1) * y_step + 1))
                    cont += 1
                end
            end
        end
    end

    # Paso 6: Graficar solución
    graficar_solucion(Turbinas_solucion,Receptores,"Imágenes/3. Exp3 Turbinas T2/Mapa Turbinas T2_divx"*string(divisions_x)*"divy_"*string(divisions_y)*"N_"*string(N)*".svg",J2)

    # Paso 7: Guardar resultados
    guardar_resultados("Experimento 3/Turbinas T2_divx"*string(divisions_x)*"divy_"*string(divisions_y)*"N_"*string(N)*".txt",N,round(Energia,digits=2),round(elapsed_time,digits=3),R2,Receptores,Turbinas_solucion)
    return Turbinas_solucion, Energia,Mz,elapsed_time
end

divisions_x=50;
divisions_y=50;
x_step=div(cols, divisions_x);
y_step=div(filas, divisions_y);
N=1000000;

Turbinas_solucion_1,Energia_1,Mz,tiempo_1=resolver_y_guardar(
    divisions_x, divisions_y,N, filas,cols, x_val1, y_val1, J2, K, E2, R2, Receptores);

Energia_resultado_exp3=Dict()
tiempo_calculo_exp3=Dict()
Turbinas_solucion_exp3=Dict()
Ruido_exp3=Dict()
Energia_resultado_exp3[1]=Energia_1
Turbinas_solucion_exp3[1]=Turbinas_solucion_1
tiempo_calculo_exp3[1]=tiempo_1
ruido_lista=[]
for k in 1:7
    push!(ruido_lista, round(10*log10(sum(10^(R2[Turbinas_solucion_1."Index_y"[i],Turbinas_solucion_1."Index_x"[i],Turbinas_solucion_1."MODO"[i]+1,k]/10) for i in 1:length(Turbinas_solucion_1."Index_x"))), digits=2))
end
Ruido_exp3[1]=ruido_lista
# Función para generar y filtrar variaciones válidas
function valid_variations(x, y,val_sep=1)
    variations = []
    # Generar y filtrar variaciones
    for (new_x, new_y) in [(x - val_sep, y), (x, y), (x + val_sep, y),
                            (x, y - val_sep), (x, y + val_sep),
                            (x - val_sep, y - val_sep), (x + val_sep, y + val_sep),
                            (x - val_sep, y + val_sep), (x + val_sep, y - val_sep)]
        if new_x >= 1 && new_x <= cols && new_y >= 1 && new_y <= filas
            push!(variations, (new_x, new_y))
        end
    end
    return variations
end






Indices=[]
for aero in eachrow(Turbinas_solucion_1)
    index_x = aero["Index_x"]
    index_y = aero["Index_y"]
    for (new_x, new_y) in valid_variations(index_x, index_y)
        push!(Indices, (new_x, new_y))
    end
end





function funcion_objetivo_iter(J2, E2, binary_decision_vars, modelo2,Indices)
    objective_sum = 0
    for j in 1:J2
        for (index,(ix,iy)) in enumerate(Indices)
            objective_sum += E2[iy,ix,j] * binary_decision_vars[index, index, j];
        end
    end
    @objective(modelo2, Max, objective_sum);
end



"""
    Un_solo_modo_por_turbina(divisions_y, divisions_x, binary_decision_vars, modelo2)

Función que agrega la restricción de que solo se puede seleccionar un modo por turbina.
"""
function Agregar_restriccion_un_solo_modo_por_turbina_iter( binary_decision_vars, modelo,num_var)
    for valor in 1:num_var
        @constraint(modelo,sum(binary_decision_vars[valor,valor,:])<=1)
    end
end
"""
    A_lo_más_n_turbinas(binary_decision_vars, modelo2,N)

Función que agrega la restricción de que se pueden seleccionar a lo más N turbinas.
"""
function Agregar_restriccion_a_lo_más_n_turbinas(binary_decision_vars, modelo2,N)
    @constraint(modelo2,sum(binary_decision_vars[:,:,:])<=N)
end
"""
    Restriccion_ruido(R2, divisions_y, divisions_x, K, binary_decision_vars, modelo2)

Función que agrega la restricción de que el ruido no puede superar los 40 dB.

"""
function Agregar_restriccion_ruido_iter(R2, Indices, K, binary_decision_vars, modelo_iter,J2)
    for k in 1:K
        @constraint(modelo_iter,sum(sum(10^(R2[iy,ix,j,k]/10)*binary_decision_vars[index,index,j] for (index,(ix,iy)) in enumerate(Indices)) for j in 1:J2)<=10^4)
    end
end
function en_poligono(idx,idy,poligono)
    return "NAME=R"*string(idy)*"R_C"*string(idx)*"C" in poligono
end

function Agregar_restriccion_valid_variations(binary_decision_vars,modelo_iter,Indices)
    for (valor,(ix,iy)) in enumerate(Indices)
        for (valor2, (ix2,iy2)) in enumerate(Indices)
            if (ix2,iy2) in valid_variations(ix,iy)
                if valor2!=valor
                    @constraint(modelo_iter,sum(binary_decision_vars[valor,valor,j]+binary_decision_vars[valor2,valor2,j] for j in 1:J2)<=1)
                end
            end
        end
    end
end


function Agregar_restriccion_distancia_entre_turbinas_iter(binary_decision_vars,modelo,Indices)
    for (valor,(ix,iy)) in enumerate(Indices)
        for (valor2, (ix2,iy2)) in enumerate(Indices)
            if valor2!=valor
                if abs(ix-ix2)<=13 && abs((iy-iy2))<=30
                    @constraint(modelo,sum(binary_decision_vars[valor,valor,j]+binary_decision_vars[valor2,valor2,j] for j in 1:J2)<=1)
                end
            end
        end
    end
end









modelo_iter = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)))
num_var=length(Indices)
@variable(modelo_iter,0<=binary_decision_vars_iter[1:num_var,1:num_var,1:J2]<=1, Bin );

funcion_objetivo_iter(J2, E2, binary_decision_vars_iter, modelo_iter,Indices)
#Agregar_restriccion_un_solo_modo_por_turbina_iter(binary_decision_vars_iter,modelo_iter,num_var)
for (valor,(ix,iy)) in enumerate(Indices)
    @constraint(modelo_iter,sum(binary_decision_vars_iter[valor,valor,:])<=Int(en_poligono(ix,filas-iy,Poligono)))
end
Agregar_restriccion_a_lo_más_n_turbinas(binary_decision_vars_iter, modelo_iter, N);
Agregar_restriccion_ruido_iter(R2, Indices, K, binary_decision_vars_iter, modelo_iter,J2)
Agregar_restriccion_distancia_entre_turbinas_iter(binary_decision_vars_iter,modelo_iter,Indices)

for (valor,(ix,iy)) in enumerate(Indices)
    @constraint(modelo_iter,sum(binary_decision_vars_iter[valor,valor,:])<=Int(en_poligono(ix,filas-iy,Poligono)))
end
    




@elapsed optimize!(modelo_iter)

Mx=value.(binary_decision_vars_iter);
Energia = objective_value(modelo_iter)

Turbinas_solucion_iter = DataFrame("X [m]" => Float64[], "Y [m]" => Float64[], "NUM" => Int[], "MODO" => Int[], "Index_x" => Int[], "Index_y" => Int[])
cont = 1
for (index,(ix,iy)) in enumerate(Indices)
    for j in 1:J2
        if Mx[index, index, j] != 0
            push!(Turbinas_solucion_iter, ((ix) * (x_val1)/cols,(iy) * y_val1/filas,cont, j-1, ix, iy))
            cont += 1
        end
    end
end

graficar_solucion(Turbinas_solucion_iter,Receptores,"Imágenes/3. Exp3 Turbinas T2/Mapa Turbinas T2_iter_divx"*string(divisions_x)*"divy_"*string(divisions_y)*"N_"*string(N)*".svg",J2)


function iteraciones(Turbinas_solucion,Poligono,J2,K,R2,N,x_step,y_step,cols,filas,x_val1,y_val1,E2,iter)
    Indices=[]
    for aero in eachrow(Turbinas_solucion)
        index_x = aero["Index_x"]
        index_y = aero["Index_y"]
        for (new_x, new_y) in valid_variations(index_x, index_y)
            push!(Indices, (new_x, new_y))
        end
    end
    modelo_iter = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)))
    num_var=length(Indices)
    @variable(modelo_iter,0<=binary_decision_vars_iter[1:num_var,1:num_var,1:J2]<=1, Bin );

    funcion_objetivo_iter(J2, E2, binary_decision_vars_iter, modelo_iter,Indices)
    #Agregar_restriccion_un_solo_modo_por_turbina_iter(binary_decision_vars_iter,modelo_iter,num_var)
    Agregar_restriccion_a_lo_más_n_turbinas(binary_decision_vars_iter, modelo_iter, N);
    Agregar_restriccion_ruido_iter(R2, Indices, K, binary_decision_vars_iter, modelo_iter,J2)
    Agregar_restriccion_distancia_entre_turbinas_iter(binary_decision_vars_iter,modelo_iter,Indices)

    for (valor,(ix,iy)) in enumerate(Indices)
        @constraint(modelo_iter,sum(binary_decision_vars_iter[valor,valor,:])<=Int(en_poligono(ix,filas-iy,Poligono)))
    end
    elapsed_time=@elapsed optimize!(modelo_iter)
    Mx=value.(binary_decision_vars_iter);
    Energia = objective_value(modelo_iter)

    Turbinas_solucion_iter = DataFrame("X [m]" => Float64[], "Y [m]" => Float64[], "NUM" => Int[], "MODO" => Int[], "Index_x" => Int[], "Index_y" => Int[])
    cont = 1
    for (index,(ix,iy)) in enumerate(Indices)
        for j in 1:J2
            if Mx[index, index, j] != 0
                push!(Turbinas_solucion_iter, ((ix) * (x_val1)/cols,(iy) * y_val1/filas,cont, j-1, ix, iy))
                cont += 1
            end
        end
    end
    ruido_lista=[]
    for k in 1:7
        push!(ruido_lista, round(10*log10(sum(10^(R2[Turbinas_solucion_iter."Index_y"[i],Turbinas_solucion_iter."Index_x"[i],Turbinas_solucion_iter."MODO"[i]+1,k]/10) for i in 1:length(Turbinas_solucion_iter."Index_x"))), digits=2))
    end

    graficar_solucion(Turbinas_solucion_iter,Receptores,"Imágenes/3. Exp3 Turbinas T2 3D 7D/Mapa Turbinas T2_3D_7D_iter_divx"*string(divisions_x)*"divy_"*string(divisions_y)*"N_"*string(N)*"_iter"*string(iter)*".svg",J2)
    return Turbinas_solucion_iter,Energia, elapsed_time,ruido_lista,size(Turbinas_solucion_iter)[1]

end

turbinas_agregadas=Dict()
turbinas_agregadas[1]=size(Turbinas_solucion_exp3[1])[1]
for i in 21:1001
    Turbinas_solucion_exp3[i],Energia_resultado_exp3[i],tiempo_calculo_exp3[i],Ruido_exp3[i],turbinas_agregadas[i]=iteraciones(Turbinas_solucion_exp3[i-1],Poligono,J2,K,R2,N,x_step,y_step,cols,filas,x_val1,y_val1,E2,i)
    if Turbinas_solucion_exp3[i-1]==Turbinas_solucion_exp3[i]
        println(i)
        break
    end
end
Guardar_resultado_50_50=DataFrame("Iteración"=>Int[],"Turbinas agregadas"=>Int[],"Energia producida [MW]"=>Float64[],"Tiempo de cálculo [s]"=>Float64[],"Ruido receptor 1 [dB]"=>Float64[],"Ruido receptor 2 [dB]"=>Float64[],"Ruido receptor 3 [dB]"=>Float64[],"Ruido receptor 4 [dB]"=>Float64[],"Ruido receptor 5 [dB]"=>Float64[],"Ruido receptor 6 [dB]"=>Float64[],"Ruido receptor 7 [dB]"=>Float64[])
for i in 1:21
    push!(Guardar_resultado_50_50,(i,turbinas_agregadas[i],round(Energia_resultado_exp3[i]*8766/1000000,digits=1),tiempo_calculo_exp3[i],Ruido_exp3[i]...))
end

