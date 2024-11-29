function selectEstudio_PSO()

    while true
        # Casos de estudio
        # Carga en el vector "caso" la lista de carpetas que hay en la carpeta de "Casos"
        listaCasos = readdir("Casos")
        # Se carga la lista y el nombre a la función de elegir opción
        casoEst = elegirOpcion(listaCasos, "caso")

        # Se elige el solver con el que se quiere resolver el AC OPF
        listaACSolvers = ["Ipopt", "Couenne"]
            s = elegirOpcion(listaACSolvers, "solver")

        # Limpieza del terminal
        limpiarTerminal()

        # Impirmir en terminal el resumen de todos las opciones elegidas
        println("Resumen:")
        println("Caso de estudio ----- ", casoEst)
        println("Optimizador --------- ", s)

        # Pregunta al usuario si las opciones listados anteriormente concuerdan con lo que quiere resolver, 
        # en caso negativo puede volver a seleccionar las opciones 
        println("\nPulsa la tecla ENTER para continuar o cualquier otra entrada para volver a elegir.")
        respuesta = readline()
        
        # Si la respuesta es un "ENTER" procede a continuar y devolver dichas opciones
        if respuesta == ""
            return casoEst,s
            break

        # En caso de introducir cualquier entrada, procede a cancelar y volver a seleccionar las opciones
        else
            continue
            
        end
    
    end

end