Un **sistema de información y alerta hidrológico** para la detección temprana de inundaciones y sequías es una **tecnología implementada**, en principio, a fin de **incrementar el tiempo de preparación para la ejecución de una respuesta apropiada (planificada)** frente a la **posible o inminente ocurrencia del peligro**. Luego, en cualquier caso debe producir información que permita responder a los siguientes **interrogantes mínimos**: 

- ¿Cuál es el estado actual de los sistemas hídricos bajo vigilancia? (monitoreo)

- ¿Cuál es la previsión de su evolución en el corto plazo? (pronóstico operativo)

Asimismo, y sobre todo si los tiempos de respuesta o la memoria del sistema son considerables, también debe proveer información, al menos orientativa, como para responder al interrogante: 

- ¿Cuál es la previsión de su evolución en el mediano y largo plazo? (  subestacional a estacional)


 >Las proyecciones hidrológicas a escala climática exceden el alcance del ámbito operativo y corresponden más bien a una componente de análisis de escenario en otros tipos de Sistemas de Soporte a la Toma de Decisión, en todo caso, del ámbito estratégico. Esto es, el requerimiento puede presentarse en el marco de proyecto de asistencia técnica (eventual) y, de manera menos usual, podría formularse el requerimiento de generación de escenarios de forma sistemática (operativa), con periodicidad relativamente constante. Por tanto, en una primera aproximación pueden obviarse. En todo caso, los resultados obtenidos en las proyecciones hidrológicas de escala climática, en conjunto con distintos escenarios territoriales proyectados, pueden utilizarse para el diseño de líneas de desarrollo y acción en los sistemas de información y alerta para la detección temprana de sequías e inundaciones.      

Para dar cuenta del **estado** de los **sistemas hídricos** en un instante específico, el sistema de información debe brindar,_estimaciones_ sobre los **principales flujos y almacenamientos** involucrados en la **fase terrestre del ciclo hidrológico**, como por ejemplo:

- Precipitación Antecedente 
- Nivel 
- Área anegada
- Almacenamiento superficial 
    - Embalses Naturales (lagos, lagunas)
    - Embalses Artificiales (presas)
    - Nieve/Hielo (cumbres, glaciares)
- Humedad del suelo
- Almacenamiento Subterráneo
- Caudal 

A fin de evaluar un mínimo de variables clave en la instancia de implementación de un sistema de información y alerta hidrológico, es conveniente elaborar un modelo conceptual de ciclo hidrológico de considerando los elementos más relevantes en las cuencas, áreas de aporte, lagos o lagunas bajo vigilancia. Así, por ejemplo, si el interés está centrado en la previsión del nivel o el caudal, mínimamente deberá contarse con información sobre precipitación antecedente, nivel y caudal. 

El **pronóstico operativo y la perspectiva subestacional a estacional** pueden realizarse sobre la base de 2 tipos de aproximación al problema:

- Procedimientos **basados en observaciones/datos**, generalmente sobre la base de **analogías en las series disponibles** (modelos empíricos o estadísticos, inteligencia artificial)

- Procedimientos **basados en modelos dinámicos del ciclo hidrológico** (modelos de transformación de precipitación en escorrentía, modelos de embalses, modelos de tránsito hidrológico o hidrodinámico).

> En caso que se especifique la generación de proyecciones a escala climática, por lo general se asume un enfoque dinámico (debido a que el alcance temporal excede el radio de correlación máximo con las observaciones hidrológicas disponibles)

Los **procedimientos operativos** deben brindar información sobre:

- La **posible o inminente ocurrencia** de un evento peligroso (alerta/aviso)

- El **tiempo de arribo** de los valores extremos del evento (pico en crecida, mínimo en estiaje)

- El **pico de intensidad** del evento (valor de pico, valor de mínimo)

- La **duración** del peligro o afectación temporal(tiempo de base por encima/por debajo de valores umbrales)

- La **extensión espacial** del peligro o afectación espacial 

Esto es, se emitirá un **alerta** si se detecta que es **posible la ocurrencia de un evento peligroso**. Luego se procederá a **establecer el tiempo disponible** en caso de que ocurra. A la vez, la **respuesta en territorio** muy posiblemente será **distinta de acuerdo a la intensidad** prevista del evento. Los primeros dos movimientos son mínimos e indispensables de acuerdo al objetivo principal de un Sistema de Alerta (incrementar el tiempo de preparación para la ejecución de una apropiada respuesta territorial). En pocas palabras, **conforme se desarrollo un sistema de alerta debe poder comunicar la severidad prevista**. Por tanto, diríamos que el tercer movimiento también es indispensable, si bien la implementación de un sistema con los primeros dos funcionales constituye un buen inicio de tareas operativas. Ciertamente, es deseable que el **esquema** utilizado para la **categorización** de la **severidad** sea **fácilmente comunicable**, y más aun que esté en acuerdo con la semiología utilizada por los agentes destinados al manejo del riesgo. 

Efectivamente, también se desea saber **cuánto durará el peligro**, puesto que esto permite estimar el daño al territorio. Y, por último, la previsión de la **superficie afectada** permite aun más precisar esta estimación. En pocas palabras, la operación del sistema debe facilitar una apropiada catarcterización del riesgo probable, y del tiempo disponible para su atenuación dada la posible o inminente ocurrencia del evento peligroso. En caso que el **peligro** sea **inminente**, se emitirá un **aviso**.  

> - La finalidad del **alerta** es brindar conocimiento sobre la **posible ocurrencia de un evento peligroso**

> - La finalidad de un **aviso** es dar cuenta de la **inminente ocurrencia**

En una primera aproximación, **la emisión de un aviso o un alerta requiere el establecimiento de valores umbrales**, a los que se denominará **umbrales territoriales**. En pocas palabras, _la inundación o la sequía es el significado social de una crecida o una bajante_. Esto es, son **eventos que estresan una estructura territorial**. Y, por tanto, dependen de umbrales propios de tal estructura. Ejemplos son la  cota mínima de actividades residenciales, productivas o de salud, la cota mínima de la infraestructura vial, la profundidad mínima para navegación, o el caudal o almacenamiento mínimo necesario para las actividades de riego o abastecimiento. Así, el _significado_ puede **trasladarse a un valor de ordenada en un hidrograma** (representación gráfica de la evolución temporal de un elemento hidrológico específico observado). Consecuentemente, **mediante el monitoreo y la previsión de hidrogramas, y a partir de la comparación con los valores umbrales, pueden emitirse alertas y avisos**.

Esto último impone **funcionalidades mínimas deseables** para el diseño de componentes de un sistema de información y alerta hidrológico. Por ejemplo, el sistema tiene que poder representar hidrogramas de las variables evaluadas del ciclo hidrológico, al menos para el monitoreo y en _tiempo operativo_. Asimismo, tiene que **facilitar la representación** de los **hidrogramas previstos**, de manera que sea relativamente sencilla su **superposición a los hidrogramas observados, incluyendo los valores umbrales**. 

Para esto, deberán implementarse procedimientos para la generación, la evaluación y ajuste de previsiones, los cuales serán de los 2 tipos ya mencionados. Asimismo, resulta necesario que el **servicio de pronóstico** sea **confiable**. Por tanto, será necesario **establecer y comunicar** la **confiabilidad** en la **emisión de alertas**, tanto como en la **categorización de la severidad** prevista. 

Definiremos intutitivamente a la **confiabilidad** como el **grado de confianza que un usuario debe asignar a un mensaje**. Analíticamente, como _la probabilidad de que el mensaje sea el correcto_. Por tanto, la _confiabilidad_ debe analizarse de acuerdo al _punto de vista del usuario_. Esto es, si bien depende de la eficiencia de los procedimientos implementados para el monitoreo y la previsión, **de acuerdo al análisis desde este punto de vista podrán ajustarse parámetros del sistema**, como por ejemplo valores umbrales, a fin de **incrementar** lo primero y, de ahí, la **confiabilidad** (por ejemplo en procedimientos sesgados).    

Los procedimientos de **pronóstico subestacional a estacional** deberán brindar información sobre ... (continuar tomando línea argumental como lo desarrollado de 35 a 47)

