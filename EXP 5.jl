#import Pkg
#Pkg.update("JuMP")
#Pkg.add("Plots")
#Pkg.status("JuMP")
#Pkg.add("DataFrames")


using JuMP, Gurobi, Plots
using DataFrames
using Printf, DelimitedFiles;
using CSV;
# Abrir un solo ambiente de Gurobi, Jupyter se marea si se crean demasiados modelos durante en una misma celda.
const GUROBI_ENV = Gurobi.Env()

Tur,headerTur=readdlm("Turbinas" * ".txt",'\t', Float64, '\n',skipstart=0,header=true);
Turbinas=DataFrame(Tur, vec(headerTur));

Rec,headerRec=readdlm("Receptores" * ".txt",'\t', Float64, '\n',skipstart=0,header=true);
Receptores=DataFrame(Rec, vec(headerRec));

Poligono=readdlm("Poligono" * ".txt",'\t', '\n',skipstart=0,header=false);
Poligono_2=readdlm("Poligono 2" * ".txt",'\t', '\n',skipstart=0,header=false);
Poligono_3=readdlm("Poligono 3" * ".txt",'\t', '\n',skipstart=0,header=false);
contador_prueba=1;

Poligono_para_grafico=DataFrame("X [m]" => Float64[], "Y [m]" => Float64[]);
for val in Poligono
    if contador_prueba%4==0
        push!(Poligono_para_grafico,(parse(Float64,Poligono[contador_prueba][1:findfirst(',',Poligono[contador_prueba])-1]),parse(Float64,Poligono[contador_prueba][findfirst(',',Poligono[contador_prueba])+1:length(Poligono[contador_prueba])])));
    end
    contador_prueba+=1
end
contador_prueba_2=1
Poligono_para_grafico_2=DataFrame("X [m]" => Float64[], "Y [m]" => Float64[]);
for val in Poligono_2
    if contador_prueba_2%4==0
        push!(Poligono_para_grafico_2,(parse(Float64,Poligono_2[contador_prueba_2][1:findfirst(',',Poligono_2[contador_prueba_2])-1]),parse(Float64,Poligono_2[contador_prueba_2][findfirst(',',Poligono_2[contador_prueba_2])+1:length(Poligono_2[contador_prueba_2])])));
    end
    contador_prueba_2+=1
end
contador_prueba_3=1
Poligono_para_grafico_3=DataFrame("X [m]" => Float64[], "Y [m]" => Float64[]);
for val in Poligono_3
    if contador_prueba_3%4==0
        push!(Poligono_para_grafico_3,(parse(Float64,Poligono_3[contador_prueba_3][1:findfirst(',',Poligono_3[contador_prueba_3])-1]),parse(Float64,Poligono_3[contador_prueba_3][findfirst(',',Poligono_3[contador_prueba_3])+1:length(Poligono_3[contador_prueba_3])])));
    end
    contador_prueba_3+=1
end


Modos_T2 = Dict() 
for i in 1:10
    Modo,header=readdlm("Ruido T2/T2 7200 Modo "*string(i-1)*".txt",'\t', '\n',header=true);
    Modos_T2[i]= DataFrame(Modo, vec(header))
end
Modos_T1 = Dict() 
for i in 1:16
    Modo,header=readdlm("Ruido T1/T1 6800 Modo "*string(i-1)*".txt",'\t', '\n',header=true);
    Modos_T1[i]= DataFrame(Modo, vec(header))
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
J2=length(Modos_T2); #Modos de turbinas T2
J1=length(Modos_T1); #Modos de turbinas T1
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
E2=zeros(Float64, filas,cols,J2);
for j in 1:J2
    E2[:,:,j]=readdlm("Mapa Energía T2/MapaEnergia T2 Modo "*string(j-1)*".csv",',',header=false);
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
function funcion_objetivo(J, divisions_y, divisions_x, E, x_step, y_step, binary_decision_vars, modelo2)
    objective_sum = 0
    for j in 1:J
        for iy in 1:divisions_y
            for ix in 1:divisions_x
                objective_sum += E[round(Int, iy * y_step - y_step + 1),round(Int, ix * x_step - x_step + 1), j] * binary_decision_vars[ix, iy, j];
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
function Agregar_restriccion_ruido(R, divisions_y, divisions_x, K, binary_decision_vars, modelo2,J)
    for k in 1:K
        @constraint(modelo2,(sum(sum(sum(10^(R[round(Int,iy*y_step-y_step+1),round(Int,ix*x_step-x_step+1),j,k]/10)*binary_decision_vars[ix,iy,j] for ix in 1:divisions_x) for iy in 1:divisions_y) for j in 1:J))<=10^4)
    end
end
function en_poligono(idx,idy,poligono)
    return "NAME=R"*string(idy)*"R_C"*string(idx)*"C" in poligono
end

function graficar_solucion(Turbinas,Receptores,Nombre,J,Poligono_para_grafico)
    paleta = palette([:green, :gray], J)
    plot(Poligono_para_grafico."X [m]",Poligono_para_grafico."Y [m]",seriestype=:scatter,label="Poligono",seriesalpha=0.05)
    plot!(Receptores."X [m]",Receptores."Y [m]",seriestype=:scatter, label="Receptores",markershape = :star,markercolor = :red)
    x_dummy = [-100, -100]  # Coordenadas X ficticias
    y_dummy = [-100, -100]  # Coordenadas Y ficticias
    modo_dummy = [0, J-1]  # Valores mínimos y máximos de MODO
    #plot!(Turbinas."X [m]",Turbinas."Y [m]",seriestype=:scatter, label="Turbinas")
    plot!([Turbinas."X [m]";x_dummy],[Turbinas."Y [m]";y_dummy],seriestype=:scatter,marker_z = [Turbinas."MODO";modo_dummy],colorbar_title = "Modo de ruido",  color = cgrad(paleta, J, categorical=true), label="Turbinas")
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
        print(round(10*log10(sum(10^(ruido[Turb."Index_y"[i],Turb."Index_x"[i],Turb."MODO"[i]+1,k]/10) for i in 1:length(Turb."Index_x"))),digits=2))
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


function resolver_y_guardar(divisions_x::Int, divisions_y::Int,N::Int, filas, cols, 
    x_val1, y_val1, J::Int, K::Int, E::Array{Float64}, R::Array{Float64}, 
    Receptores::DataFrame,Poligono_para_grafico::DataFrame,Poligono)

    # Paso 1: Calcular pasos en X e Y
    x_step=floor(cols/divisions_x);
    y_step=floor(filas/divisions_y);

    # Paso 2: Crear modelo de optimización
    modelo2 = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)))
    @variable(modelo2, 0 <= binary_decision_vars[1:divisions_x, 1:divisions_y, 1:J] <= 1, Bin)

    # Paso 3: Agregar función objetivo y restricciones
    funcion_objetivo(J, divisions_y, divisions_x, E, x_step, y_step, binary_decision_vars, modelo2)
    Agregar_restriccion_un_solo_modo_por_turbina(divisions_y, divisions_x, binary_decision_vars, modelo2)
    Agregar_restriccion_a_lo_más_n_turbinas(binary_decision_vars, modelo2, N)
    Agregar_restriccion_ruido(R, divisions_y, divisions_x, K, binary_decision_vars, modelo2,J)
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
            for j in 1:J
                if Mz[ix, iy, j] != 0
                    push!(Turbinas_solucion, ((((ix - 1) * x_step + 1) * (x_val1) / cols),( ((iy - 1) * y_step + 1) * (y_val1) / filas),cont, j-1, (ix - 1) * x_step + 1, (iy - 1) * y_step + 1))
                    cont += 1
                end
            end
        end
    end

    # Paso 6: Graficar solución
    graficar_solucion(Turbinas_solucion,Receptores,"Imágenes/5. Exp5 Turbinas T1 Poligono 2/Mapa Turbinas T1_divx"*string(divisions_x)*"divy_"*string(divisions_y)*"N_"*string(N)*".svg",J,Poligono_para_grafico)

    # Paso 7: Guardar resultados
    guardar_resultados("Experimento 5/Poligono 2_Turbinas T1_divx"*string(divisions_x)*"divy_"*string(divisions_y)*"N_"*string(N)*".txt",N,round(Energia,digits=2),round(elapsed_time,digits=3),R,Receptores,Turbinas_solucion)
    return Turbinas_solucion, Energia*8766/1000000,Mz,elapsed_time
end
divisions_x=50;
divisions_y=20;
x_step=div(cols, divisions_x);
y_step=div(filas, divisions_y);
N=700;

Turbinas_solucion_1_Exp5,Energia_1_exp5,Mz_exp5,tiempo_1_exp5=resolver_y_guardar(
    divisions_x, divisions_y,N, filas,cols, x_val1, y_val1, J1, K, E1, R1, Receptores,Poligono_para_grafico_2,Poligono_2);

Energia_resultado_exp5=Dict()
tiempo_calculo_exp5=Dict()
Turbinas_solucion_exp5=Dict()
Ruido_exp5=Dict()
Energia_resultado_exp5[1]=Energia_1_exp5
Turbinas_solucion_exp5[1]=Turbinas_solucion_1_Exp5
tiempo_calculo_exp5[1]=tiempo_1_exp5
ruido_lista_Exp5=[]
for k in 1:7
    push!(ruido_lista_Exp5, round(10*log10(sum(10^(R1[Turbinas_solucion_1_Exp5."Index_y"[i],Turbinas_solucion_1_Exp5."Index_x"[i],Turbinas_solucion_1_Exp5."MODO"[i]+1,k]/10) for i in 1:length(Turbinas_solucion_1_Exp5."Index_x"))), digits=2))
end
Ruido_exp5[1]=ruido_lista_Exp5

# Función para generar y filtrar variaciones válidas
function valid_variations(x, y,val_sep=1)
    variations = []
    # Generar y filtrar variaciones
    for (new_x, new_y) in [(x - val_sep, y), (x, y), (x + val_sep, y),
                            (x, y - val_sep), (x, y + val_sep),
                            (x - val_sep, y - val_sep), (x + val_sep, y + val_sep),
                            (x - val_sep, y + val_sep), (x + val_sep, y - val_sep)]
        if new_x >= 1 && new_x <= filas && new_y >= 1 && new_y <= cols
            push!(variations, (new_x, new_y))
        end
    end
    return variations
end






Indices=[]
for aero in eachrow(Turbinas_solucion_1_Exp5)
    index_x = aero["Index_x"]
    index_y = aero["Index_y"]
    for (new_x, new_y) in valid_variations(index_x, index_y)
        push!(Indices, (new_x, new_y))
    end
end





function funcion_objetivo_iter_exp4(J, E, binary_decision_vars, modelo2,Indices)
    objective_sum = 0
    for (index,(ix,iy)) in enumerate(Indices)
        for j in 2:J
            objective_sum += E[iy,ix,j] * (binary_decision_vars[index, index, j]-binary_decision_vars[index, index, j-1]);
        end
        objective_sum += E[iy,ix,1]*binary_decision_vars[index, index, 1]
    end
    @objective(modelo2, Max, objective_sum);
end



"""
    Un_solo_modo_por_turbina(divisions_y, divisions_x, binary_decision_vars, modelo2)

Función que agrega la restricción de que solo se puede seleccionar un modo por turbina.
"""
function Agregar_restriccion_un_solo_modo_por_turbina_iter_exp4( binary_decision_vars, modelo,num_var,J)
    for valor in 1:num_var
        for j in 1:(J-1)
            @constraint(modelo,(binary_decision_vars[valor,valor,j]-binary_decision_vars[valor,valor,j+1])<=0)
        end
    end
end
"""
    A_lo_más_n_turbinas(binary_decision_vars, modelo2,N)

Función que agrega la restricción de que se pueden seleccionar a lo más N turbinas.
"""
function Agregar_restriccion_a_lo_más_n_turbinas_exp4(binary_decision_vars, modelo,N,J)
    suma=0
    for valor in 1:num_var
        suma+=binary_decision_vars[valor,valor,J]
    end
    @constraint(modelo,suma<=N)
end
"""
    Restriccion_ruido(R2, divisions_y, divisions_x, K, binary_decision_vars, modelo2)

Función que agrega la restricción de que el ruido no puede superar los 40 dB.

"""
function Agregar_restriccion_ruido_iter_exp4(R, Indices, K, binary_decision_vars, modelo_iter,J)
    for k in 1:K
        suma_ruido=0
        for (index, (ix,iy)) in enumerate(Indices)
            for j in 2:J
                suma_ruido+=10^(R[iy,ix,j,k]/10)*(binary_decision_vars[index,index,j]-binary_decision_vars[index,index,j-1])
            end
            suma_ruido+=10^(R[iy,ix,1,k]/10)*binary_decision_vars[index,index,1]
        end
        @constraint(modelo_iter,suma_ruido<=10^4)
    end
end
function en_poligono(idx,idy,poligono)
    return "NAME=R"*string(idy)*"R_C"*string(idx)*"C" in poligono
end

function Agregar_restriccion_valid_variations_exp4(binary_decision_vars,modelo_iter,Indices,J)
    for (valor,(ix,iy)) in enumerate(Indices)
        for (valor2, (ix2,iy2)) in enumerate(Indices)
            if (ix2,iy2) in valid_variations(ix,iy)
                if valor2!=valor
                    @constraint(modelo_iter,binary_decision_vars[valor,valor,J]+binary_decision_vars[valor2,valor2,J]<=1)
                end
            end
        end
    end
end





function Agregar_restriccion_distancia_entre_turbinas_iter_exp4(binary_decision_vars,modelo,Indices,J)
    for (valor,(ix,iy)) in enumerate(Indices)
        for (valor2, (ix2,iy2)) in enumerate(Indices)
            if valor2!=valor
                if abs(ix-ix2)<=13 && abs((iy-iy2))<=30
                    @constraint(modelo,binary_decision_vars[valor,valor,J]+binary_decision_vars[valor2,valor2,J]<=1)
                end
            end
        end
    end
end






function iteraciones_exp5(Turbinas_solucion,Poligono,J,K,R,N,filas,cols,x_val1,y_val1,E,iter,Poligono_para_grafico)
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
    @variable(modelo_iter,0<=binary_decision_vars_iter[1:num_var,1:num_var,1:J]<=1, Bin );
    #set_optimizer_attribute(modelo_iter, "TimeLimit", 600) 
    funcion_objetivo_iter_exp4(J, E, binary_decision_vars_iter, modelo_iter,Indices)
    Agregar_restriccion_un_solo_modo_por_turbina_iter_exp4(binary_decision_vars_iter,modelo_iter,num_var,J)
    #Agregar_restriccion_a_lo_más_n_turbinas_exp4(binary_decision_vars_iter, modelo_iter, N);
    Agregar_restriccion_ruido_iter_exp4(R, Indices, K, binary_decision_vars_iter, modelo_iter,J)
    Agregar_restriccion_distancia_entre_turbinas_iter_exp4(binary_decision_vars_iter,modelo_iter,Indices,J)
    for (valor,(ix,iy)) in enumerate(Indices)
        @constraint(modelo_iter,binary_decision_vars_iter[valor,valor,J]<=Int(en_poligono(ix,filas-iy,Poligono)))
    end
    elapsed_time=@elapsed optimize!(modelo_iter)
    Mx=value.(binary_decision_vars_iter);
    Energia = objective_value(modelo_iter)

    Turbinas_solucion_iter = DataFrame("X [m]" => Float64[], "Y [m]" => Float64[], "NUM" => Int[], "MODO" => Int[], "Index_x" => Int[], "Index_y" => Int[])
    cont = 1
    for (index,(ix,iy)) in enumerate(Indices)
        for j in 0:(J-2)
            if (Mx[index, index, J-j]- Mx[index, index, J-1-j])!= 0 
                push!(Turbinas_solucion_iter, ((ix) * (x_val1)/cols,(iy) * y_val1/filas,cont, J-j-1, ix, iy))
                cont += 1
            end
        end
        if Mx[index, index, 1]!= 0
            push!(Turbinas_solucion_iter, ((ix) * (x_val1)/cols,(iy) * y_val1/filas,cont, 0, ix, iy))
            cont += 1
        end
    end
    ruido_lista_Exp5=[]
    for k in 1:7
        push!(ruido_lista_Exp5, round(10*log10(sum(10^(R[Turbinas_solucion_iter."Index_y"[i],Turbinas_solucion_iter."Index_x"[i],Turbinas_solucion_iter."MODO"[i]+1,k]/10) for i in 1:length(Turbinas_solucion_iter."Index_x"))), digits=2))
    end
    println("Iteración $iter")
    #graficar_solucion(Turbinas_solucion_iter,Receptores,"Imágenes/5. Exp5 Turbinas T1 Poligono 2/Mapa Turbinas T1_3D_7D_iter_divx"*string(divisions_x)*"divy_"*string(divisions_y)*"N_"*string(N)*"_iter"*string(iter)*".svg",J,Poligono_para_grafico)
    return Turbinas_solucion_iter,Energia, elapsed_time,ruido_lista_Exp5,size(Turbinas_solucion_iter)[1] 
end



turbinas_agregadas_exp5=Dict()
turbinas_agregadas_exp5[1]=size(Turbinas_solucion_exp5[1])[1]
for i in 2:1001
    Turbinas_solucion_exp5[i],Energia_resultado_exp5[i],tiempo_calculo_exp5[i],Ruido_exp5[i],turbinas_agregadas_exp5[i]=iteraciones_exp5(Turbinas_solucion_exp5[i-1],Poligono_2,J1,K,R1,N,filas,cols,x_val1,y_val1,E1,i,Poligono_para_grafico_2)
    if Turbinas_solucion_exp5[i-1]==Turbinas_solucion_exp5[i]
        println(i)
        break
    end
end

Guardar_resultado_50_50_exp5_T1=DataFrame("Iteración"=>Int[],"Turbinas agregadas"=>Int[],"Energia producida [MW]"=>Float64[],"Tiempo de cálculo [s]"=>Float64[],"Ruido receptor 1 [dB]"=>Float64[],"Ruido receptor 2 [dB]"=>Float64[],"Ruido receptor 3 [dB]"=>Float64[],"Ruido receptor 4 [dB]"=>Float64[],"Ruido receptor 5 [dB]"=>Float64[],"Ruido receptor 6 [dB]"=>Float64[],"Ruido receptor 7 [dB]"=>Float64[])
for i in 1:1001
    push!(Guardar_resultado_50_50_exp5_T1,(i,turbinas_agregadas_exp5[i],Energia_resultado_exp5[i],tiempo_calculo_exp5[i],Ruido_exp5[i]...))
end
CSV.write("Experimento 5/Resultados 1001 50x20 Exp5 T1 Poligono 2.csv",Guardar_resultado_50_50_exp5_T1)
CSV.write("Experimento 5/Turbinas 50x20 Exp5 T1 100 Poligono 2.csv",Turbinas_solucion_exp5[101])
CSV.write("Experimento 5/Turbinas 50x20 Exp5 T1 1000 Poligono 2.csv",Turbinas_solucion_exp5[1001])
CSV.write("Experimento 5/Turbinas 50x20 Exp5 T1 0 Poligono 2.csv",Turbinas_solucion_exp5[1])
CSV.write("Experimento 5/Turbinas 50x20 Exp5 T1 1 Poligono 2.csv",Turbinas_solucion_exp5[2])
CSV.write("Experimento 5/Turbinas 50x20 Exp5 T1 50 Poligono 2.csv",Turbinas_solucion_exp5[51])
CSV.write("Experimento 5/Turbinas 50x20 Exp5 T1 700 Poligono 2.csv",Turbinas_solucion_exp5[701])

for i in [1,2,10,11,50,51,100,101,1001]
    graficar_solucion(Turbinas_solucion_exp5[i],Receptores,"Imágenes/5. Exp5 Turbinas T1 Poligono 2/Mapa Turbinas T1_3D_7D_iter_divx"*string(divisions_x)*"divy_"*string(divisions_y)*"N_"*string(N)*"_iter"*string(i)*".svg",J1,Poligono_para_grafico_2)
end