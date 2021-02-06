# Draft of Social Service Project 0.

# Posibles títulos.

- Estrategias de recolección eficiente de ensembles farmacóforicos en 3 dimensiones extraidos de la simulación de dinámica molecular con el objectivo de incrementar el éxito de un cribado high-throuput screening.
- Implementación de una librería de código abierto para la extracción de ensembles farmacofóricos 3D de la simulación dinámica molecular para el cribado virtual en la búsqueda de nuevos ligandos.
- Implementación de una librería de código abierto para la extracción de ensembles farmacofóricos 3D mediante la simulación dinámica molecular.

# Introducción

Un farmacóforo se usa como filtro de un cribado virtual.

# Unidad didáctica

## Objectivos
- Elaboración de jupyter notebooks como documento resumen que será colgado en web de librería.
- Elaboración de slides para charla.

## Contenidos
- El alumno debe, mediante el estudio guiado de la bibliografía científica, generar un discurso coherente, consistente y autocontenido, a modo de unidad didáctica, que contenga las respuestas a las siguientes cuestiones:

- ¿Qué es un modelo farmacofórico 3D y cuál es su utilidad en el campo del diseño y descubribimiento de pequeñas moléculas con potencial farmacológico?
- ¿Qué propiedades debe incluir un modelo farmacofórico 3D?, o ¿qué es un punto farmacofórico?.
- Qué diferencias hay entre el modelado farmacofórico basado en estructura y el basado en ligando?
- Qué es una estructura lead o un lead.
- Que software y librerías hay para crear modelos 3D farmacofóricos, libre y propietario.
- Que sets de compuestos, bencharmks o sistemas se usan para el testeo y comparación de software.
- ¿En qué consisten las estrategias de obtención de puntos farmacofóricos mediante la simulación de dinámica molecular del blanco junto con una alta concentración de pequeñas moléculas orgánicas como sondas químicas?
- Dado un sitio de unión en un blanco y varios ligandos que se unen, qué es el farmacóforo consenso o farmacóforo de propiedades compartidas. Cómo se obtiene?
- Puede un mismo ligando sobre distintos blancos generar un farmacóforo consenso? Cual sería su sentido y utilidad?
- Puede el farmacóforo consenso provenir de varias conformaciones de un mismo sitio de unión? Cúal sería su sentido y utilidad?
- Qué son los algoritmos de clustering basados en farmacóforos y para qué se usan? Y los algoritmos de alineamiento?
- Siguiendo el abordaje tradicional del diseño de ligandos basado en estructura, dado un blanco, el farmacóforo usado para el cribado de una quimioteca se obtiene de las estructuras resueltas experimentalente, por difracción de rayos X o por NMR en su mayoría. Qué sucedería si el cribado se llevara acabo con un farmacóforo consendo o un conjunto de farmacóforos extraidos de la dinámica molecular del blanco? Hay evidencias de que esta estrategía modifique el resultado del cribado? Cómo lo modificará? En que sistemas de test se ha probado este abordaje.
- Qué estrategias de construcción de "ensembles" de farmacóforos extraidos de la dinámica molecular de un sistem receptor-ligando, o de un receptor, encontramos en la literatura reciente.
- Hay margen para estrategias nuevas de exploración eficiente mediante, por ejemplo técnicas de sampleo adaptativo de dinámica molecular?
- Ejemplos en la literatura del uso exitoso de modelos farmacofóricos para la predicción de estructuras lead en proyectos de diseño o búsqueda de ligandos con potencial farmacológico.
- En el caso de que a través de la dinámica molecular podamos construir una red cinética o un modelo de estados de markov de los modelos farmacofóricos observados en el sistema. Qué nos interesa más... usar los farmacóforos más probables (con menor energía libre) o aquellos cuya meta-estabilidad (o tiempo promedio) de vida sea mayor? Por qué?
- Cómo podemos enriquecer, haciendo uso del simulación de dinámica molecular, un modelo farmacofórico en el caso de que sea posible la tautomerización de un grupo químico en la región que define dicho modelo?
- Qué información nos daría el estudio de la familia de modelos farmacofóricos obtenidos de realizar diferentes mutaciones puntuales en la región del modelo? Para qué nos puede servir esta información? Qué relación puede tener esto con el diseño sobre el conjunto de isformas de un blanco -si las hubiere-?

# Librería
- De la estructura química al farmacóforo
- Elaboración de documentación en web

El resultado será, dado un blanco o un sistema receptor-ligando/s. Obtendremos la red cinética o modelo de estados de markov de los modelos farmacofóricos. Lo cual nos permitirá obtener los distintos modelos farmcofóricos que el sistema visita, así como la probabilidad, sus tiempos de residencia y sus tasas de intercambio.

# Milestones
- Charla en el 3 mes con slides de unidad didáctica.
- Charla en el 6 mes con extensión de unidad didáctica a herramientas desarrolladas en la librería.

# Entregables

- Documento en Latex (unidades didácticas) con anexo de slides de las charlas.
- Librería en GitHub, donde se puede seguir la implementación así como las discusiones.
- Deployment de la librería.
- Documentación de la librería en web.
- Reporte del tutor de actividades adicionales como la asistencia a seminarios del laboratorio.

# Obligaciones adicionales
- Asistencia a seminarios virtuales del laboratorio sobre el trabajo realizado por los investigadores y estudiantes del laboratorio.
- Asistencia a talleres virtuales sobre cuestiones de interés para su formación. 

# Calendario
- 240 horas para la unidad didáctica en el primer trimestre
- 240 horas para la elaboración de la librería en el segundo trimestre

# Qué aprendera el alumno

## El contexto científico del proyecto

- Herramientas y estratégias de búsqueda de bibliografía científica.
- Herramientas de gestión biobliográfica.
- Se hará uso de la librería OpenMM para llevar a cabo una pequeña simulación de dinámica molecular con la que demostrar el uso de la librería implementada. La función de este experimento no será la investigación, sino simplemente el obtener una o varias trayectorias sobre las que demostrar que la librería obtiene el ensemble relevante de modelos farmacofóricos.

## El trabajo en un laboratorio computacional de investigación

- Linux.
- Python y librerías científicas.
- Conda y Anaconda.
- Git, control de versiones.
- GitHub, desarrollo en colaboración de herramientas científicas.
- Jupyter lab and notebooks.
- Latex.
- Deployment of packages with conda.
- Creación de documentación científica con sphinx.
- Lenguages decorados: Markdown y rst.
