using JuMP, Gurobi, Plots
using DataFrames
using Printf, DelimitedFiles;
# Abrir un solo ambiente de Gurobi, Jupyter se marea si se crean demasiados modelos durante en una misma celda.
const GUROBI_ENV = Gurobi.Env()

Tur,headerTur=readdlm("Turbinas" * ".txt",'\t', Float64, '\n',skipstart=0,header=true);
Turbinas=DataFrame(Tur, vec(headerTur));

Rec,headerRec=readdlm("Receptores" * ".txt",'\t', Float64, '\n',skipstart=0,header=true);
Receptores=DataFrame(Rec, vec(headerRec));

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
J2=length(Modos_T2); #Modos de turbinas T2
K=length(Rec[:,1]); #Receptores

vel=readdlm("Velocidad Viento" * ".txt",'\t', Float64, '\n',skipstart=5);
x_val0=0;
y_val0=0;
x_val1=28600;
y_val1=23760;
filas,cols=size(vel)
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
                objective_sum += E2[round(Int, iy * y_step - y_step + 1),round(Int, ix * x_step - x_step + 1) , j] * binary_decision_vars[iy, ix, j];
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
            @constraint(modelo2,sum(binary_decision_vars[iy,ix,:])<=1);
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
        suma_ruido_T2=0
        for ix in 1:divisions_x
            for iy in 1:divisions_y
                for j in 1:J2
                    suma_ruido_T2+=(10^((R2[round(Int,iy*y_step-y_step+1),round(Int,ix*x_step-x_step+1),j,k])/10))*binary_decision_vars[iy,ix,j] 
                end
            end
        end
        @constraint(modelo2,suma_ruido_T2<=10^4)
    end
end

function graficar_solucion(Turbinas,Receptores,Nombre)
    plot(Poligono_para_grafico."X [m]",Poligono_para_grafico."Y [m]",seriestype=:scatter,label="Poligono",seriesalpha=0.01)
    plot!(Receptores."X [m]",Receptores."Y [m]",seriestype=:scatter, label="Receptores",markershape = :star,markercolor = :red)
    plot!(Turbinas."X [m]",Turbinas."Y [m]",seriestype=:scatter, label="Turbinas")
    xlims!(x_val0-500, x_val1+500)  # Cambiar los límites del eje X
    ylims!(y_val0-500, y_val1+500)   # Cambiar los límites del eje Y
    annotate!(Receptores."X [m]".+10,Receptores."Y [m]".-1000,text.(Int.(Receptores."Receptor"),:down,14))
    savefig(Nombre)
end

function mostrar_ruido_receptores(ruido,Recep,Turb)
    for k in 1:length(Recep[:,1])
        print("El ruido que siente el receptor ")
        print(k)
        print(" es ")
        suma_ruido=0
        for i in 1:length(Turb."Index_x")
            ruido_tur=ruido[Turb."Index_y"[i],Turb."Index_x"[i],Turb."MODO"[i]+1,k]
            suma_ruido+=(10^(ruido_tur/10))
        end
        print(round(10*log10(suma_ruido),digits=2))
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
    x_val1, y_val1, J2::Int, K::Int, E2::Array{Float64}, R2::Array{Float64}, 
    Receptores::DataFrame)
    # Paso 1: Calcular pasos en X e Y
    x_step = Int(round(cols/divisions_x))
    y_step = Int(round(filas/divisions_y))
    # Paso 2: Crear modelo de optimización
    modelo2 = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)))
    @variable(modelo2, 0 <= binary_decision_vars[1:divisions_y, 1:divisions_x, 1:J2] <= 1, Bin)
    # Paso 3: Agregar función objetivo y restricciones
    funcion_objetivo(J2, divisions_y, divisions_x, E2, x_step, y_step, binary_decision_vars, modelo2)
    Agregar_restriccion_un_solo_modo_por_turbina(divisions_y, divisions_x, binary_decision_vars, modelo2)
    Agregar_restriccion_a_lo_más_n_turbinas(binary_decision_vars, modelo2, N)
    suma_ruido_2=Dict()
    for k in 1:K
        suma_ruido_T2=0
        for ix in 1:divisions_x
            for iy in 1:divisions_y
                for j in 1:J2
                    suma_ruido_T2+=(10^(R2[round(Int,iy*y_step-y_step+1),round(Int,ix*x_step-x_step+1),j,k]/10))*binary_decision_vars[iy,ix,j] 
                end
            end
        end
        suma_ruido_2[k]=suma_ruido_T2
        @constraint(modelo2,suma_ruido_2[k]<=10000)
    end
    # Paso 4: Resolver modelo
    elapsed_time = @elapsed optimize!(modelo2)
    Energia = objective_value(modelo2)
    Mz = value.(binary_decision_vars)
    Ruido_por_receptor=Dict()
    for k in 1:K
        Ruido_por_receptor[k]=value(suma_ruido_2[k])
    end
    #   Paso 5: Extraer solución
    Turbinas_solucion = DataFrame("X [m]" => Float64[], "Y [m]" => Float64[], "NUM" => Int[], "MODO" => Int[], "Index_x" => Int[], "Index_y" => Int[])
    cont = 1
    for ix in 1:divisions_x
        for iy in 1:divisions_y
            for j in 1:J2
                if Mz[ix, iy, j] != 0
                    push!(Turbinas_solucion, ((((ix - 1) * x_step + 1) * (x_val1) / cols),(((iy - 1) * y_step + 1) * (y_val1) / filas),cont, j-1, round(Int,ix*x_step-x_step+1), round(Int,iy*y_step-y_step+1)))
                    cont += 1
                end
            end
        end
    end

    # Paso 6: Graficar solución
    graficar_solucion(Turbinas_solucion,Receptores,"Imágenes/1. Exp1 Turbinas T2/Mapa Turbinas T2_divx"*string(divisions_x)*"divy_"*string(divisions_y)*"N_"*string(N)*".svg")

    # Paso 7: Guardar resultados
    guardar_resultados("Experimento 1/Turbinas T2_divx"*string(divisions_x)*"divy_"*string(divisions_y)*"N_"*string(N)*".txt",N,round(Energia,digits=2),round(elapsed_time,digits=3),R2,Receptores,Turbinas_solucion)
    return Turbinas_solucion, Ruido_por_receptor
end

divisions_x=50;
divisions_y=50;
N=2500;

Turbinas_solucion,ruido_receptor=resolver_y_guardar(
    divisions_x, divisions_y,N, filas, cols, x_val1, y_val1, J2, K, E2, R2, Receptores);









