function selectCaracteristicas()
    include("./PSO/extraerDatos_PSO.jl")
    include("./PSO/selectEstudio_PSO.jl")

    casoEstudio, s = selectEstudio_PSO()
    datos = extraerDatos_PSO(casoEstudio)

    # El número de partículas será el número de generadores
    nDim = nrow(datos[2])
    nParticle = 0
    nInter = 0

    while true
    
        while true

            # Entra en un bloque try-catch para poder manejar las entradas que provocan excepciones en el sistema
            try 

                limpiarTerminal()

                # Elegir el número de soluciones que explora cada partícula
                println("\nNúmero de generadores = ", nDim)
                println("\nEscriba el número de soluciones que quiere que explore cada partícula (debe ser mayor al número de generadores)")
                nParticle = parse(Int, readline())

                if nParticle > nDim
                    break

                # En caso de que el número de soluciones exploradas sea menor al número de partículas
                else 
                    println("\n ERROR: número de soluciones a explorar menor o igual al número de generadores.")
                    sleep(2)
                    limpiarTerminal()
                    #println("\nNúmero de generadores = ", nDim)
                    #println("Escriba el número de soluciones que quiere que explore cada partícula (debe ser mayor al número de generadores)")
                    #nParticle = parse(Int, readline())
                end

                # En caso de que la entrada cause una excepción, 
                # por ejemplo introduciendo una letra al cual no se puede convertir en un int
            catch

                # Limpia el terminal
                limpiarTerminal()

                # El mensaje se muestra en pantalla por 2 segundos
                println("Entrada no válida. Por favor, introduzca un número.")
                sleep(2)
                continue

            end
        end

        while true

            try
                # Elegir el número de iteraciones del PSO
                println("\nEscriba el número de veces que quieres que se itere el algoritmo PSO: ")
                nInter = parse(Int, readline())

                break # Salgo del bucle 

                # En caso de que la entrada cause una excepción, 
                # por ejemplo introduciendo una letra al cual no se puede convertir en un int
            catch
                # Limpia el terminal
                limpiarTerminal()

                # El mensaje se muestra en pantalla por 2 segundos
                println("Entrada no válida. Por favor, introduzca un número.")
                sleep(2)
                continue

            end
        end

        # Resumen de los datos elegidos 
        println("Resumen:")
        println("Número de generadores ----- ", nDim)
        println("Número soluciones --------- ", nParticle)
        println("Número de iteraciones ----- ", nInter)

        println("\nPulsa la tecla ENTER para continuar o cualquier otra entrada para volver a elegir.")
        respuesta = readline()

        if respuesta == ""
            return nDim, nParticle, nInter, casoEstudio, s
            break 

        else
            continue      
        end
    end
end
