using CSV
using DataFrames
include("./evaluarTensiones.jl")
include("./evaluarFlujos.jl")
include("./calcularAdmitancias.jl")

mutable struct ParticleHibrida
    nGeneradores::Int
    
    # Estado de los generadores (u)
    position_u::Array{Float64, 1}  # Valores continuos [0,1]
    velocity_u::Array{Float64, 1}
    pBest_u::Array{Float64, 1}
    lBest_u::Array{Float64, 1}
    
    # Potencias de generadores (PG)
    position_pg::Array{Float64, 2}  # Matriz nx2 (P y Q)
    velocity_pg::Array{Float64, 2}
    pBest_pg::Array{Float64, 2}
    lBest_pg::Array{Float64, 2}
    
    # Valores de fitness
    fitValue::Float64
    fitpBest::Float64
    fitlBest::Float64
    
    nFitEval::Int
    
    # Función para inicializar las partículas
    function ParticleHibrida(nGeneradores::Int, datosGenerador::DataFrame) 
        # Inicialización de estados (u)
        position_u = rand(nGeneradores)  # Valores continuos entre 0 y 1
        velocity_u = rand(nGeneradores) .- 0.5  # Velocidades iniciales centradas en 0
        pBest_u = copy(position_u)
        lBest_u = copy(position_u)
        
        # Inicialización de potencias, matrices a 0
        position_pg = zeros(Float64, nGeneradores, 2) # Creación de una matriz de ceros con 1 columna para P y 1 columna para Q
        velocity_pg = zeros(Float64, nGeneradores, 2)
        
        # Inicializar potencias dentro de límites para todos los generadores
        # Aquí se dan valores aleatorios dentro de los límites de potencia de los generadores 
        # para P y Q y se almacenan en la matriz position_pg inicializada previamente 
        for i in 1:nGeneradores
            # Potencia activa
            Pmin = datosGenerador.P_MIN[i]
            Pmax = datosGenerador.P_MAX[i]
            position_pg[i,1] = Pmin + rand()*(Pmax - Pmin) # Fila i, columna 1
            
            # Potencia reactiva
            Qmin = datosGenerador.Q_MIN[i]
            Qmax = datosGenerador.Q_MAX[i]
            position_pg[i,2] = Qmin + rand()*(Qmax - Qmin) # Fila i, columna 2
        end
        
        pBest_pg = copy(position_pg) # La mejor posición de la partícula es la posición actual
        lBest_pg = copy(position_pg) # La mejor posición de los vecinos es la posición actual
        
        # Se dan valores infinitos a los fitness para tratar de minimizarlos
        fitValue = Inf
        fitpBest = Inf
        fitlBest = Inf
        
        # Se inicializa el número de evaluaciones de fitness a 0
        nFitEval = 0
        
        # Se crea el objeto ParticleHibrida con los valores iniciales
        new(nGeneradores, position_u, velocity_u, pBest_u, lBest_u,
            position_pg, velocity_pg, pBest_pg, lBest_pg,
            fitValue, fitpBest, fitlBest, nFitEval)
    end       
end

# Función para actualizar la posición de las partículas
function updatePosition!(p::ParticleHibrida, w::Float64, c1::Float64, c2::Float64, datos::Tuple)
    datosGenerador = datos[2]
    
    # Actualizar velocidades de u
    p.velocity_u = w * p.velocity_u + 
                  c1 * rand() * (p.pBest_u - p.position_u) + 
                  c2 * rand() * (p.lBest_u - p.position_u)
    
    # Actualizar posiciones de u
    p.position_u += p.velocity_u
    p.position_u = clamp.(p.position_u, 0.0, 1.0)
    
    # Actualizar velocidades y posiciones de PG según máscara
    for i in 1:p.nGeneradores
        if p.position_u[i] >= 0.5
            # Actualizar velocidad y posición solo si el generador está activo
            p.velocity_pg[i,:] = w * p.velocity_pg[i,:] + 
                                c1 * rand() * (p.pBest_pg[i,:] - p.position_pg[i,:]) + 
                                c2 * rand() * (p.lBest_pg[i,:] - p.position_pg[i,:])
            
            p.position_pg[i,:] += p.velocity_pg[i,:]
            
            # Limitar potencias a sus rangos. 
            # Clamp es una función que limita un valor a un rango
            # Si el valor está por debajo del mínimo, se pone el mínimo. 
            # Si está por encima del máximo, se pone el máximo.
            p.position_pg[i,1] = clamp(p.position_pg[i,1], 
                                     datosGenerador.P_MIN[i], 
                                     datosGenerador.P_MAX[i])
            p.position_pg[i,2] = clamp(p.position_pg[i,2], 
                                     datosGenerador.Q_MIN[i], 
                                     datosGenerador.Q_MAX[i])
        else
            # Si el generador está inactivo, solo poner velocidad a cero para que no se mueva la partícula
            p.velocity_pg[i,:] .= 0.0
            # Las potencias mantienen sus valores anteriores
        end
    end
end

mutable struct SwarmHibrido
    fitFunc::Function
    nGeneradores::Int
    datos::Tuple
    
    nParticle::Int
    nNeibor::Int
    nInter::Int
    
    c1::Float
    c2::Float
    
    wMax::Float
    wMin::Float
    w::Float # Valor de inercia
    
    gBest_u::Array{Float64, 1}
    gBest_pg::Array{Float64, 2}
    fitgBest::Float64
    
    particles::Array{ParticleHibrida, 1}
    
    nFitEvals::Int
    
    # Función para inicializar el enjambre (partículas, estado de los generadores y potencias)
    function SwarmHibrido(fitFunc::Function, nGeneradores::Int, datos::Tuple;
            nParticle::Int=3, nNeibor::Int=3, nInter::Int=2000,
            c1::Float=2.0, c2::Float=2.0,
            wMax::Float=0.9, wMin::Float=0.4)
        
        if nNeibor > nParticle
            error("El número de partículas en un grupo local no debe exceder el número total de partículas")
        end
        
        w = wMax # Valor de inercia inicial
        println("w: ", w)
        datosGenerador = datos[2] # Datos de los generadores sacado del csv
        
        # Inicializar partículas con la función ParticleHibrida 	
        particles = [ParticleHibrida(nGeneradores, datosGenerador) for i in 1:nParticle]
        
        # Inicializar mejores posiciones globales
        gBest_u = rand(nGeneradores) # Vector de 0 y 1 
        gBest_pg = zeros(Float64, nGeneradores, 2) # Matriz de nGeneradores x 2 de ceros 
        fitgBest = Inf # Valor infinito para tratar de minimizarlo
        
        nFitEvals = 0 # Número de evaluaciones de fitness
        
        # Se crea el objeto SwarmHibrido con los valores iniciales
        new(fitFunc, nGeneradores, datos, nParticle, nNeibor, nInter,
            c1, c2, wMax, wMin, w, gBest_u, gBest_pg, fitgBest,
            particles, nFitEvals)
    end
end

# Función para inicializar el fitness de la partícula
function initFitValue!(p::ParticleHibrida, fitFunc::Function, datos::Tuple, log_file::Union{IOStream, Nothing}=nothing, log_enabled::Bool=false)
    println("initFitValue!")
    p.fitValue, _ = fitFunc(p, datos, log_file, log_enabled)
    p.nFitEval += 1
    nothing
end

# Función para evaluar el fitness de las partículas. Se usa durante la optimización para actualizar las partículas
function evaluate!(p::ParticleHibrida, fitFunc::Function, datos::Tuple, log_file::Union{IOStream, Nothing}, log_enabled::Bool)
    p.fitValue, _ = fitFunc(p, datos, log_file, log_enabled) 
    p.nFitEval += 1 # Aumenta el número de evaluaciones de fitness
    
    # Actualizar mejor personal si corresponde
    if p.fitValue < p.fitpBest # Si el fitness de la partícula es menor que el mejor fitness personal, se actualiza
        p.fitpBest = p.fitValue
        p.pBest_u = copy(p.position_u)
        p.pBest_pg = copy(p.position_pg)
    end
    
    # Se escribe el valor del fitness en el archivo de log que se crea en la carpeta logs
    if log_enabled
        log_to_file(log_file, "fitValue: $(p.fitValue)", log_enabled)
    end
    
    return p
end

# Función para actualizar el mejor fitness de la partícula si corresponde
function updatepBestAndFitpBest!(p::ParticleHibrida)
    if p.fitValue < p.fitpBest
        p.fitpBest = p.fitValue
        p.pBest_u = copy(p.position_u)
        p.pBest_pg = copy(p.position_pg)
    end
    nothing
end

# Función para actualizar la inercia del enjambre
function updateInertia!(s::SwarmHibrido)
    dw = (s.wMax - s.wMin)/s.nInter
    s.w -= dw
    nothing
end

# Función para actualizar la mejor posición global y el mejor fitness global
function updategBestAndFitgBest!(s::SwarmHibrido)
    gFits = [particle.fitValue for particle in s.particles]
    fitgBest, index = findmin(gFits) # Selecciona el menor fitness de todas las partículas y su índice
    
    # Si el fitness de la partícula seleccionada es menor que el mejor fitness global, se actualiza
    if fitgBest < s.fitgBest
        s.gBest_u = copy(s.particles[index].position_u)
        s.gBest_pg = copy(s.particles[index].position_pg)
        s.fitgBest = fitgBest
    end
    nothing
end

# Función para actualizar la mejor posición local y el mejor fitness local
function updatelBestAndFitlBest!(s::SwarmHibrido)
    for i in 1:s.nParticle
        neiborIds = neiborIndices(i, s.nNeibor, s.nParticle)
        neiborFits = [s.particles[Id].fitValue for Id in neiborIds]
        fitlBest, index = findmin(neiborFits)
        
        if fitlBest < s.particles[i].fitlBest
            s.particles[i].lBest_u = copy(s.particles[neiborIds[index]].position_u)
            s.particles[i].lBest_pg = copy(s.particles[neiborIds[index]].position_pg)
            s.particles[i].fitlBest = fitlBest
        end
    end
    nothing
end

# Función para inicializar el enjambre
function initialize!(s::SwarmHibrido, log_file::Union{IOStream, Nothing}=nothing, log_enabled::Bool=false)
    for particle in s.particles
        initFitValue!(particle, s.fitFunc, s.datos, log_file, log_enabled)
        updatepBestAndFitpBest!(particle)
    end
    
    updatelBestAndFitlBest!(s)
    updategBestAndFitgBest!(s)
    
    return s
end

function optimize!(s::SwarmHibrido, log_file::Union{IOStream, Nothing}, log_enabled::Bool)
    log_to_file(log_file, "\nIniciando PSO híbrido", log_enabled)
    mejor_fitness_historico = Inf
    iteraciones_sin_mejora = 0
    
    for i in 1:s.nInter
        log_to_file(log_file, "\n=== Iteración $i ===", log_enabled)
        log_to_file(log_file, "Inercia actual: $(s.w)", log_enabled)
        
        for (j, p) in enumerate(s.particles)
            log_to_file(log_file, "\nActualizando partícula $j", log_enabled)
            updatePosition!(p, s.w, s.c1, s.c2, s.datos)
            evaluate!(p, s.fitFunc, s.datos, log_file, log_enabled)
        end
        
        updatelBestAndFitlBest!(s)
        updategBestAndFitgBest!(s)
        
        if s.fitgBest < mejor_fitness_historico
            mejora = mejor_fitness_historico - s.fitgBest
            log_to_file(log_file, "\nMejora encontrada: $mejora", log_enabled)
            mejor_fitness_historico = s.fitgBest
            iteraciones_sin_mejora = 0
        else
            iteraciones_sin_mejora += 1
        end
        
        updateInertia!(s)
        log_to_file(log_file, "\nMejor fitness actual: $(s.fitgBest)", log_enabled)
        log_to_file(log_file, "Iteraciones sin mejora: $iteraciones_sin_mejora", log_enabled)
    end
    
    return s.gBest_u, s.gBest_pg, s.fitgBest
end

# Función que corresponde a la fitFunc, evalúa la partícula
function evaluarParticula(p::ParticleHibrida, datos::Tuple, log_file::Union{IOStream, Nothing}=nothing, log_enabled::Bool=false)
    datosLinea, datosGenerador, datosNodo, nNodos, nLineas, bMVA = datos
    
    # Calcular matriz de admitancias y valores relacionados
    Y_sparse, Y_serie, Y_sh = calcularAdmitancias(datosLinea, nNodos, nLineas)
    Y = Matrix(Y_sparse)
    
    # Convertir las admitancias serie y shunt a vectores
    y_series = [1.0 / (datosLinea.R[i] + im * datosLinea.X[i]) for i in 1:nLineas]
    y_shunts = [0.5 * (im * datosLinea.BSh[i]) for i in 1:nLineas]
    
    # Extraer potencias activas y reactivas
    potencias_P = p.position_pg[:,1]
    potencias_Q = p.position_pg[:,2]
    estados_u = p.position_u  # Estado de los generadores
    
    # Calcular potencia total generada (damos un valor inicial por si todos los generadores están apagados)
    potencia_total_generada = sum(potencias_P[i] for i in 1:p.nGeneradores if estados_u[i] >= 0.5; init=0.0)
    demanda_total = sum(datosNodo.PD)
    
    # Si no se cubre la demanda, retornar coste infinito
    if potencia_total_generada < demanda_total
        return Inf, zeros(nNodos)
    end
    
    # Evaluar tensiones y violaciones
    println("evaluarTensiones!")
    V, violaciones_tension = evaluarTensiones(datosLinea, datosGenerador, datosNodo,
                                            nNodos, nLineas, float(bMVA), 
                                            potencias_P, potencias_Q, estados_u,
                                            Y, log_file, log_enabled)
    
    # Evaluar flujos y sus violaciones
    println("evaluarFlujos!")
    violaciones_flujo = evaluarFlujos(datosLinea, Y, V, float(bMVA),
                                     y_series, y_shunts,  # Ahora son vectores
                                     log_file, log_enabled)
    
    # Calcular coste de generación    
    coste = 0.0
    for i in 1:p.nGeneradores
        if p.position_u[i] >= 0.5 # Usar position_u para determinar si el generador está activo
            coste += datosGenerador.P_COSTE1[i] * potencias_P[i]
        end
    end
    println("coste: ", coste)
    
    # Penalización por violaciones de tensión
    coste_total = coste + 1000 * violaciones_tension + 1000 * violaciones_flujo
    
    return coste_total, V
end

function neiborIndices(i::Int, nNeibor::Int, nParticle::Int)
    nNeibor = max(3, nNeibor)
    nLeft = (nNeibor - 1) ÷ 2
    startIndex = (i - nLeft)
    endIndex = startIndex + nNeibor - 1
    indices = collect(startIndex:endIndex)
    
    for i in 1:nNeibor
        if indices[i] < 1
            indices[i] += nParticle
        elseif indices[i] > nParticle
            indices[i] -= nParticle
        end
    end
    
    indices
end

function log_to_file(log_file::Union{IOStream, Nothing}, message::String, log_enabled::Bool)
    if log_enabled && !isnothing(log_file)
        write(log_file, message * "\n")
    end
    println(message)
end

function initialize_log(caso_estudio::String, log_enabled::Bool)
    if log_enabled
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        log_filename = "logs/PSO_$(caso_estudio)_$(timestamp).log"
        mkpath("logs")  # Crear directorio si no existe
        return open(log_filename, "w")
    end
    return nothing
end

function runPSOHibrido(datos::Tuple, nParticle::Int, nInter::Int, log_enabled::Bool)
    log_file = initialize_log(datos[end], log_enabled)  # datos[end] debería ser caso_estudio
    
    nGeneradores = size(datos[2], 1)
    log_to_file(log_file, "Iniciando PSO Híbrido con $nParticle partículas y $nInter iteraciones", log_enabled)
    
    # Crear enjambre
    swarm = SwarmHibrido(evaluarParticula, nGeneradores, datos, 
                        nParticle=nParticle, nInter=nInter)
    
    # Inicializar y ejecutar
    initialize!(swarm)
    gBest_u, gBest_pg, fitgBest = optimize!(swarm, log_file, log_enabled)
    
    if log_enabled
        close(log_file)
    end
    
    return gBest_u, gBest_pg, fitgBest
end 