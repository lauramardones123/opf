using CSV
using DataFrames
include("./evaluarTensiones.jl")

mutable struct ParticleHibrida
    nGeneradores::Int
    
    # Posiciones binarias (on/off de generadores)
    position_bin::Array{Float64, 1}
    velocity_bin::Array{Float64, 1}
    pBest_bin::Array{Float64, 1}
    lBest_bin::Array{Float64, 1}
    
    # Posiciones continuas (potencias de generadores)
    position_pot::Array{Float64, 2}  # Matriz nx2 (P y Q para cada generador)
    velocity_pot::Array{Float64, 2}
    pBest_pot::Array{Float64, 2}
    lBest_pot::Array{Float64, 2}
    
    # Valores de fitness
    fitValue::Float64
    fitpBest::Float64
    fitlBest::Float64
    
    nFitEval::Int
    
    function ParticleHibrida(nGeneradores::Int, datosGenerador::DataFrame) 
        # Inicialización binaria
        position_bin = rand(0:1.0, nGeneradores)
        velocity_bin = rand(nGeneradores) - position_bin
        pBest_bin = copy(position_bin)
        lBest_bin = copy(position_bin)
        
        # Inicialización de potencias
        position_pot = zeros(Float64, nGeneradores, 2)
        velocity_pot = zeros(Float64, nGeneradores, 2)
        
        # Inicializar potencias dentro de límites
        for i in 1:nGeneradores
            # Potencia activa
            Pmin = datosGenerador.P_MIN[i]
            Pmax = datosGenerador.P_MAX[i]
            position_pot[i,1] = Pmin + rand()*(Pmax - Pmin)
            
            # Potencia reactiva
            Qmin = datosGenerador.Q_MIN[i]
            Qmax = datosGenerador.Q_MAX[i]
            position_pot[i,2] = Qmin + rand()*(Qmax - Qmin)
        end
        
        pBest_pot = copy(position_pot)
        lBest_pot = copy(position_pot)
        
        fitValue = Inf
        fitpBest = Inf
        fitlBest = Inf
        
        nFitEval = 0
        
        new(nGeneradores, position_bin, velocity_bin, pBest_bin, lBest_bin,
            position_pot, velocity_pot, pBest_pot, lBest_pot,
            fitValue, fitpBest, fitlBest, nFitEval)
        # println("ParticleHibrida creada")
    end       

    function updatePosition!(p::ParticleHibrida, w::Float, c1::Float, c2::Float, datos::Tuple)
        # Actualizar velocidades
        p.velocity_bin = w * p.velocity_bin + 
                        c1 * rand() * (p.pBest_bin - p.position_bin) + 
                        c2 * rand() * (p.lBest_bin - p.position_bin)
        
        p.velocity_pot = w * p.velocity_pot + 
                        c1 * rand() * (p.pBest_pot - p.position_pot) + 
                        c2 * rand() * (p.lBest_pot - p.position_pot)
        
        # Actualizar posiciones
        p.position_bin += p.velocity_bin
        p.position_pot += p.velocity_pot
        
        # Discretizar binarias
        for i in 1:p.nGeneradores
            p.position_bin[i] = p.position_bin[i] >= 0.5 ? 1.0 : 0.0
        end
        
        # Limitar potencias
        datosGenerador = datos[2]
        for i in 1:p.nGeneradores
            p.position_pot[i,1] = clamp(p.position_pot[i,1], 
                                      datosGenerador.P_MIN[i], 
                                      datosGenerador.P_MAX[i])
            p.position_pot[i,2] = clamp(p.position_pot[i,2], 
                                      datosGenerador.Q_MIN[i], 
                                      datosGenerador.Q_MAX[i])
        end
        return p
    end

    function evaluate!(p::ParticleHibrida, fitFunc::Function, datos::Tuple)
        p.fitValue, _ = fitFunc(p, datos)
        p.nFitEval += 1
        
        # Actualizar mejor personal si corresponde
        if p.fitValue < p.fitpBest
            p.fitpBest = p.fitValue
            p.pBest_bin = copy(p.position_bin)
            p.pBest_pot = copy(p.position_pot)
        end
        return p
    end

    function updateBest!(p::ParticleHibrida)
        if p.fitValue < p.fitpBest
            p.fitpBest = p.fitValue
            p.pBest_bin = copy(p.position_bin)
            p.pBest_pot = copy(p.position_pot)
        end
        return p
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
    
    gBest_bin::Array{Float64, 1}    
    gBest_pot::Array{Float64, 2}
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
        gBest_bin = rand(0:1.0, nGeneradores)
        gBest_pot = zeros(Float64, nGeneradores, 2)
        fitgBest = Inf
        
        nFitEvals = 0
        
        new(fitFunc, nGeneradores, datos, nParticle, nNeibor, nInter,
            c1, c2, wMax, wMin, w, gBest_bin, gBest_pot, fitgBest,
            particles, nFitEvals)
    end
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
        p.pBest_bin = copy(p.position_bin)
        p.pBest_pot = copy(p.position_pot)
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
        s.gBest_bin = copy(s.particles[index].position_bin)
        s.gBest_pot = copy(s.particles[index].position_pot)
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
            s.particles[i].lBest_bin = copy(s.particles[neiborIds[index]].position_bin)
            s.particles[i].lBest_pot = copy(s.particles[neiborIds[index]].position_pot)
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
        
        for (j, p) in enumerate(s.particles)
            updateParticle!(p, s.w, s.c1, s.c2, s.datos)
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
    
    return s.gBest_bin, s.gBest_pot, s.fitgBest
end

function evaluarParticula(p::ParticleHibrida, datos::Tuple)
    datosLinea, datosGenerador, datosNodo, nNodos, nLineas, bMVA = datos
    
    # Extraer potencias activas y reactivas
    potencias_P = p.position_pot[:,1]
    potencias_Q = p.position_pot[:,2]
    
    # Evaluar tensiones y violaciones
    println("evaluarTensiones!")
    V_mag, violaciones = evaluarTensiones(datosLinea, datosGenerador, datosNodo,
                                         nNodos, nLineas, float(bMVA), p.position_bin, potencias_P, potencias_Q)
    # V_mag, violaciones = evaluarTensiones(  # Usar el nombre completo del módulo
    #                                     datosLinea, datosGenerador, datosNodo, 
    #                                     nNodos, nLineas, bMVA,
    #                                     p.position_bin, potencias_P, potencias_Q)  
    
    # Calcular coste de generación    
    coste = 0.0
    for i in 1:p.nGeneradores
        if p.position_bin[i] == 1
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
    optimize!(swarm)
end 