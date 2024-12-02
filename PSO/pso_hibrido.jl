using CSV
using DataFrames
include("./evaluarTensiones.jl")

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
    
    function ParticleHibrida(nGeneradores::Int, datosGenerador::DataFrame) 
        # Inicialización de estados (u)
        position_u = rand(nGeneradores)  # Valores continuos entre 0 y 1
        velocity_u = rand(nGeneradores) .- 0.5  # Velocidades iniciales centradas en 0
        pBest_u = copy(position_u)
        lBest_u = copy(position_u)
        
        # Inicialización de potencias
        position_pg = zeros(Float64, nGeneradores, 2)
        velocity_pg = zeros(Float64, nGeneradores, 2)
        
        # Inicializar potencias dentro de límites para todos los generadores
        for i in 1:nGeneradores
            # Potencia activa
            Pmin = datosGenerador.P_MIN[i]
            Pmax = datosGenerador.P_MAX[i]
            position_pg[i,1] = Pmin + rand()*(Pmax - Pmin)
            
            # Potencia reactiva
            Qmin = datosGenerador.Q_MIN[i]
            Qmax = datosGenerador.Q_MAX[i]
            position_pg[i,2] = Qmin + rand()*(Qmax - Qmin)
        end
        
        pBest_pg = copy(position_pg)
        lBest_pg = copy(position_pg)
        
        fitValue = Inf
        fitpBest = Inf
        fitlBest = Inf
        
        nFitEval = 0
        
        new(nGeneradores, position_u, velocity_u, pBest_u, lBest_u,
            position_pg, velocity_pg, pBest_pg, lBest_pg,
            fitValue, fitpBest, fitlBest, nFitEval)
    end       
end

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
            
            # Limitar potencias a sus rangos
            p.position_pg[i,1] = clamp(p.position_pg[i,1], 
                                     datosGenerador.P_MIN[i], 
                                     datosGenerador.P_MAX[i])
            p.position_pg[i,2] = clamp(p.position_pg[i,2], 
                                     datosGenerador.Q_MIN[i], 
                                     datosGenerador.Q_MAX[i])
        else
            # Si el generador está inactivo, solo poner velocidad a cero
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
    w::Float
    
    gBest_u::Array{Float64, 1}
    gBest_pg::Array{Float64, 2}
    fitgBest::Float64
    
    particles::Array{ParticleHibrida, 1}
    
    nFitEvals::Int
    
    function SwarmHibrido(fitFunc::Function, nGeneradores::Int, datos::Tuple;
            nParticle::Int=3, nNeibor::Int=3, nInter::Int=2000,
            c1::Float=2.0, c2::Float=2.0,
            wMax::Float=0.9, wMin::Float=0.4)
        
        if nNeibor > nParticle
            error("El número de partículas en un grupo local no debe exceder el número total de partículas")
        end
        
        w = wMax
        println("w: ", w)
        datosGenerador = datos[2]
        # println("datosGenerador: ", datosGenerador)
        
        # Inicializar partículas
        particles = [ParticleHibrida(nGeneradores, datosGenerador) for i in 1:nParticle]
        
        # Inicializar mejores globales
        gBest_u = rand(nGeneradores)
        gBest_pg = zeros(Float64, nGeneradores, 2)
        fitgBest = Inf
        
        nFitEvals = 0
        
        new(fitFunc, nGeneradores, datos, nParticle, nNeibor, nInter,
            c1, c2, wMax, wMin, w, gBest_u, gBest_pg, fitgBest,
            particles, nFitEvals)
    end
end

function evaluate!(p::ParticleHibrida, fitFunc::Function, datos::Tuple)
    p.fitValue, _ = fitFunc(p, datos)
    p.nFitEval += 1
    
    # Actualizar mejor personal si corresponde
    if p.fitValue < p.fitpBest
        p.fitpBest = p.fitValue
        p.pBest_u = copy(p.position_u)
        p.pBest_pg = copy(p.position_pg)
    end
    return p
end

function updateBest!(p::ParticleHibrida)
    if p.fitValue < p.fitpBest
        p.fitpBest = p.fitValue
        p.pBest_u = copy(p.position_u)
        p.pBest_pg = copy(p.position_pg)
    end
    return p
end

function initFitValue!(p::ParticleHibrida, fitFunc::Function, datos::Tuple)
    println("initFitValue!")
    p.fitValue, _ = fitFunc(p, datos)
    p.nFitEval += 1
    nothing
end

function updatepBestAndFitpBest!(p::ParticleHibrida)
    if p.fitValue < p.fitpBest
        p.fitpBest = p.fitValue
        p.pBest_u = copy(p.position_u)
        p.pBest_pg = copy(p.position_pg)
    end
    nothing
end

function updateInertia!(s::SwarmHibrido)
    dw = (s.wMax - s.wMin)/s.nInter
    s.w -= dw
    nothing
end

function updategBestAndFitgBest!(s::SwarmHibrido)
    gFits = [particle.fitValue for particle in s.particles]
    fitgBest, index = findmin(gFits)
    
    if fitgBest < s.fitgBest
        s.gBest_u = copy(s.particles[index].position_u)
        s.gBest_pg = copy(s.particles[index].position_pg)
        s.fitgBest = fitgBest
    end
    nothing
end

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

function initialize!(s::SwarmHibrido)
    for particle in s.particles
        initFitValue!(particle, s.fitFunc, s.datos)
        updatepBestAndFitpBest!(particle)
    end
    
    updatelBestAndFitlBest!(s)
    updategBestAndFitgBest!(s)
    
    return s
end

function optimize!(s::SwarmHibrido)
    println("\nIniciando PSO híbrido")
    mejor_fitness_historico = Inf
    iteraciones_sin_mejora = 0
    
    for i in 1:s.nInter
        println("\n=== Iteración $i ===")
        println("Inercia actual: ", s.w)
        
        for p in s.particles
            println("UpdatePosition!")
            # updatePosition!()
            updatePosition!(p, s.w, s.c1, s.c2, s.datos)
            evaluate!(p, s.fitFunc, s.datos)
        end
        
        updatelBestAndFitlBest!(s)
        updategBestAndFitgBest!(s)
        
        if s.fitgBest < mejor_fitness_historico
            mejora = mejor_fitness_historico - s.fitgBest
            mejor_fitness_historico = s.fitgBest
            iteraciones_sin_mejora = 0
        else
            iteraciones_sin_mejora += 1
        end
        
        updateInertia!(s)
        
        # if i % 100 == 0
        #     println("\n=== Resumen iteración $i ===")
        #     println("Mejor fitness actual: ", s.fitgBest)
        #     println("Iteraciones sin mejora: ", iteraciones_sin_mejora)
        # end
    end
    
    # println("\n=== Optimización completada ===")
    # println("Mejor fitness encontrado: ", s.fitgBest)
    
    return s.gBest_u, s.gBest_pg, s.fitgBest
end

function evaluarParticula(p::ParticleHibrida, datos::Tuple)
    datosLinea, datosGenerador, datosNodo, nNodos, nLineas, bMVA = datos
    
    # Extraer potencias activas y reactivas
    potencias_P = p.position_pg[:,1]
    potencias_Q = p.position_pg[:,2]
    estados_u = p.position_u  # Estado de los generadores
    
    # Evaluar tensiones y violaciones
    println("evaluarTensiones!")
    V_mag, violaciones = evaluarTensiones(datosLinea, datosGenerador, datosNodo,
                                         nNodos, nLineas, float(bMVA), 
                                         potencias_P, potencias_Q, estados_u)
    
    # Calcular coste de generación    
    coste = 0.0
    for i in 1:p.nGeneradores
        if p.position_u[i] >= 0.5  # Usar position_u para determinar si el generador está activo
            coste += datosGenerador.P_COSTE1[i] * potencias_P[i]
        end
    end
    println("coste: ", coste)
    
    # Penalización por violaciones de tensión
    coste_total = coste + 1000 * violaciones
    
    return coste_total, V_mag
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


function runPSOHibrido(datos::Tuple, nParticle::Int, nInter::Int)
    nGeneradores = size(datos[2], 1)
    
    # Crear enjambre
    swarm = SwarmHibrido(evaluarParticula, nGeneradores, datos, 
                        nParticle=nParticle, nInter=nInter)
    
    # Inicializar y ejecutar
    initialize!(swarm)
    gBest_u, gBest_pg, fitgBest = optimize!(swarm)
    return gBest_u, gBest_pg, fitgBest
end 