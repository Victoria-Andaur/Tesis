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
Modos_T1 = Dict() 
for i in 1:16
    Modo,header=readdlm("Ruido T1/T1 6800 Modo "*string(i-1)*".txt",'\t', '\n',header=true);
    Modos_T1[i]= DataFrame(Modo, vec(header))
end


Energias_T2 = Dict() 
for i in 1:10
    Energia,header_energia=readdlm("Energía T2/T2 7200 Modo "*string(i-1)*" Energia.txt",'\t', '\n',header=true);
    Energias_T2[i]= DataFrame(Energia, vec(header_energia))
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
filas,cols=size(vel);
densidad=1.192

E1=zeros(Float64, filas, cols,J1);
for j in 1:J1
    E1[:,:,j]=readdlm("Mapa Energía T1/MapaEnergia T1 Modo "*string(j-1)*".csv",',',header=false);
end



E2=zeros(Float64, filas, cols,J2);
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
    funcion_objetivo(J, divisions_y, divisions_x, E, x_step, y_step, binary_decision_vars, modelo2)
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
    Restriccion_ruido(R, divisions_y, divisions_x, K, binary_decision_vars, modelo2)

Función que agrega la restricción de que el ruido no puede superar los 40 dB.
"""
function Agregar_restriccion_ruido(R, divisions_y, divisions_x, K, binary_decision_vars, modelo2,x_step,y_step,J)
    for k in 1:K
        @constraint(modelo2,(sum(sum(sum(10^(R[round(Int,iy*y_step-y_step+1),round(Int,ix*x_step-x_step+1),j,k]/10)*binary_decision_vars[ix,iy,j] for ix in 1:divisions_x) for iy in 1:divisions_y) for j in 1:J))<=10^4)
    end
end


function en_poligono(idx,idy,poligono)
    return "NAME=R"*string(idy)*"R_C"*string(idx)*"C" in poligono
end

function graficar_solucion(Turbinas,Receptores,Nombre,J)
    paleta = palette([:green, :gray], J)
    plot(Poligono_para_grafico."X [m]",Poligono_para_grafico."Y [m]",seriestype=:scatter,label="Poligono",seriesalpha=0.05)
    plot!(Receptores."X [m]",Receptores."Y [m]",seriestype=:scatter, label="Receptores",markershape = :star,markercolor = :red)
    x_dummy = [-100, -100]  # Coordenadas X ficticias
    y_dummy = [-100, -100]  # Coordenadas Y ficticias
    modo_dummy = [0, J-1]  # Valores mínimos y máximos de MODO
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




function resolver_y_guardar(divisions_x::Int, divisions_y::Int,N::Int, filas, cols, x_val1, y_val1, J::Int, K::Int, E::Array{Float64}, R::Array{Float64}, 
    Receptores::DataFrame)

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
    Agregar_restriccion_ruido(R, divisions_y, divisions_x, K, binary_decision_vars, modelo2,x_step,y_step,J)
    for ix in 1:divisions_x
        for iy in 1:divisions_y
            @constraint(modelo2,sum(binary_decision_vars[ix,iy,:])<=Int(en_poligono(round(Int,ix*x_step-x_step+1),round(Int,filas-(iy*y_step-y_step+1)),Poligono)))
        end
    end
    #Agregar_restriccion_distancia_entre_turbinas(binary_decision_vars,modelo2,x_step,y_step,divisions_x,divisions_y)
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
    graficar_solucion(Turbinas_solucion,Receptores,"Imágenes/1. Exp1 Turbinas T1/Mapa Turbinas T1_divx"*string(divisions_x)*"divy_"*string(divisions_y)*"N_"*string(N)*".svg",J)

    # Paso 7: Guardar resultados
    guardar_resultados("Experimento 1/Turbinas T1_divx"*string(divisions_x)*"divy_"*string(divisions_y)*"N_"*string(N)*".txt",N,round(Energia,digits=2),round(elapsed_time,digits=3),R,Receptores,Turbinas_solucion)

    return Turbinas_solucion, Energia,Mz,elapsed_time
end

divisions_x=50;
divisions_y=50;
#x_step=div(cols, divisions_x);
#x_step_m=div(x_val1,divisions_x);
#y_step=div(filas, divisions_y);
#y_step_m=div(y_val1,divisions_y);

Energia_resultado_exp1=Dict()
Ruido_receptores=Dict()
tiempo_calculo=Dict()
[5,10,25,30,50,80,100,125,150,500,625,100,1000]
for n in [5,10,25,30,50,80,100,125,150,500,625,1000]
    N=n;
    Turbinas_solucion,Energia_resultado,Mz,elapsed_time=resolver_y_guardar(
    divisions_x, divisions_y,N, filas, cols, 
     x_val1, y_val1, J1, K, E1, R1, Receptores)
    ruido_lista=[]
    Energia_resultado_exp1[N]=Energia_resultado
    for k in 1:7
        push!(ruido_lista, round(10*log10(sum(10^(R[Turbinas_solucion."Index_y"[i],Turbinas_solucion."Index_x"[i],Turbinas_solucion."MODO"[i]+1,k]/10) for i in 1:length(Turbinas_solucion."Index_x"))), digits=2))
    end
    Ruido_receptores[N]=ruido_lista
    tiempo_calculo[N]=elapsed_time
end
tiempo_calculo3030=[]
energia3030=[]
rr=DataFrame("1" => Float64[], "2" => Float64[], "3" => Float64[], "4" => Float64[], "5" => Float64[], "6" => Float64[],"7" => Float64[])
for n in [5,10,25,30,50,80,100,125,150,500,625]
    push!(rr,(Ruido_receptores[n][1],Ruido_receptores[n][2],Ruido_receptores[n][3],Ruido_receptores[n][4],Ruido_receptores[n][5],Ruido_receptores[n][6],Ruido_receptores[n][7]))
    push!(tiempo_calculo3030,tiempo_calculo[n])
    push!(energia3030,Energia_resultado_exp1[n])
end
writedlm("Experimento 2/Turbinas T2 tiempo_calculo3030.csv", tiempo_calculo3030)
writedlm("Experimento 2/Turbinas T2 Energia_3030.csv", energia3030)
CSV.write("Experimento 2/Turbinas T2 ruido_receptores3030.csv", rr)