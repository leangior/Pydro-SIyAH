# Modelos Matriciales de Tránsito Hidrológico
Notas Técnicas en Hidrología Operativa

Dr. Leandro Giordano

SIyAH, Instituto Nacional del Agua

## Modelo lineal para un tramo compuesto por N-subtramos en serie

Supóngase un tramo fluvial compuesto por N sub-tramos fluviales. Los sub-tramos están definidos entre secciones transversales, asociadas al ingreso desde aguas arriba y a la descarga hacia aguas abajo, dentro del tramo. Se asumirá que el sistema es abierto, esto es: con seguridad, el primer tramo recibirá aporte desde un borde superior y, además, se acepta que cada sub-tramo puede recibir aporte lateral. Asimismo, la descarga del sistema aguas abajo, a través de su frontera, será igual a la descarga del último sub-tramo. Consecuentemente, por conservación de volumen, la descarga de cada sub-tramo se encontrará en función del ingreso desde aguas arriba (sub-tramos antecedentes o borde superior), del aporte lateral y del almacenamiento inicial en estos. La Figura 1 muestra la topología impuesta y su representación como sistema de reservorios lineales. 

![Topología de tramo fluvial con $N$-subtramos en serie. Representación del tramo como sistema de reservorios lineales en serie. $K$ es la proporción de agua libre que permanece en cada subtramo y $(1-K)$ la proporción que constituye descarga aguas abajo](SistemaMatricialReservoriosv2.png)

Partiremos de la simple asunción que la descarga $Q_i(j)$ del $i$-ésimo sub-tramo, durante el $j$-ésimo intervalo regular de tiempo de duración $\Delta t$, puede asociarse de forma lineal al almacenamiento $X_i(j)$. Por otro lado, el sistema es abierto. Esto es, puede presentar bordes laterales en cada tramo y, con seguridad, presenta al menos un borde superior (tramo inicial). Esto es, presenta aportes variables en el tiempo $I_i(j)$ .  

Así, imponiendo linealidad e incluyendo los aportes variables por la frontera del sistema en cada sub-tramo, se asumirá que la descarga durante cada intervalo será igual a una proporción de la suma del almacenamiento inicial y de los aportes de los bordes 

$$Q_i(j)=(1-K_i)[X_i(j)+I_i(j)]$$

En donde $K_i$ es la proporción del agua libre que permanece en el $i$-ésimo sub- tramo durante el intervalo. En principio, esta proporción se asumirá constante entre los intervalos de tiempo (imposición de linealidad) y se considerará que $0<K_i<=1$ (conservación del volumen). Ciertamente, en estas condiciones podrá asumirse que el vector de parámetros $\mathbf{K}$ está asociado al tiempo de residencia medio en cada sub-tramo. 

Así, el modelo queda definido por tres elementos básicos:

- Un vector de parámetros $\mathbf{K}=(K_1,K_2,...,K_N)$
- Un vector de condiciones iniciales $\mathbf{X}=(X_1(j),X_2(j),..,X_N(j))$.
- Una matriz de condiciones de borde/forzantes $\mathbf{I}$, en donde cada entrada $I_{ij}=I_i(j)$ se corresponderá al valor del aporte de borde al $i$-ésimo tramo para el $j$-ésimo intervalo (matriz de forzantes)

Luego si se desarrolla la ecuación de conservación para cada tramo, primeramente podrá observarse que debe cumplirse (Fig. 1):

$$X_{i}(j+1)-X_{i}(j)=Q_{i-1}(j)+I_i(j)-Q_i(j)$$

En donde el subíndice $i$ indica el número y orden de cada sub-tramo, desde aguas arriba a aguas abajo. Luego, desarrollando $Q_i$ y $Q_{i-1}$, de acuerdo a su definición, y realizando las operaciones corrrespondientes, se obtiene:
$$X_{i}(j+1)=(1-K_{i-1})[X_{i-1}(j)+I_{i-1}(j)]+K_i[X_i(j)+I_i(j)]$$
Que para $i=1$ se reduce a:
$$X_{1}(j+1)=K[X_i(j)+I_i(j)]$$
puesto que es el tramo inicial, por lo que se asume que $Q_0(j)=0$. 

Por tanto podrá notarse que, para cada $j$-ésimo intervalo de cálculo, la operación del sistema, que afecta los estados $\mathbf{X}$,  quedará definida por la transformación lineal:
$$
% \begin{align}
\mathbf{A}\,\mathbf{[X}(j)+\mathbf{I}(j)\mathbf{]}=\mathbf{[X}(j+1){]}
% \end{align}
$$
En donde $\mathbf{A}$ es la matriz de transformación del sistema dinámico.  Para el caso de un sistema sub-tramos o reservorios lineales en serie (Fig. 1), cada uno con constante de prorateo $K_i$, con $i=1,2,...,N$, esto es:
$$
% \begin{align*}
\mathbf{A=}
    \begin{bmatrix}
    K_1 & 0 & 0 & .& .& . & 0 \\
    1 - K_1 & K_2 & 0 & .& .& . & 0 \\
    0 & 1-K_2 & K_3 & .& .& . & 0 \\
    . & . & . & .& .& . & . \\
    . & . & . & .& .& . & . \\
    . & . & . & .& .& . & . \\
    . & . & . & .& .& . & . \\
    . & . & . & .& .& . & . \\
    0 & . & . & .& .& 1-K_{N-1} & K_N \\
    \end{bmatrix}
% \end{align*}
$$
Por tanto, para el caso evalaudo el modelo dinámico es:
$$
% \begin{align*}
    \begin{bmatrix}
    K_1 & 0 & 0 & .& .& . & 0 \\
    1 - K_1 & K_2 & 0 & .& .& . & 0 \\
    0 & 1-K_2 & K_3 & .& .& . & 0 \\
    . & . & . & .& .& . & . \\
    . & . & . & .& .& . & . \\
    . & . & . & .& .& . & . \\
    . & . & . & .& .& . & . \\
    . & . & . & .& .& . & . \\
    0 & . & . & .& .& 1-K_{N-1} & K_N \\
    \end{bmatrix}
    \begin{bmatrix}
    X_1(j) + I_1(j) \\
    X_2(j) + I_2(j) \\
    X_3(j) + I_3(j) \\
    . \\
    . \\
    . \\
    . \\
    . \\
    X_N(j) + I_N(j) \\
    \end{bmatrix}
    =
    \begin{bmatrix}
    X_1(j+1) \\
    X_2(j+1)  \\
    X_3(j+1)  \\
    . \\
    . \\
    . \\
    . \\
    . \\
    X_N(j+1)  \\
    \end{bmatrix}
% \end{align*}
$$ 

>
>

El subíndice $i$ indica el orden topológico de cada sub-tramo en la serie, desde aguas arriba hacia aguas abajo. El subíndice $j$ indica el orden temporal, para pasos de cálculo de longitud $\Delta t$, duración que debe ser significativamente mayor al tiempo de residencia medio. La matriz $\mathbf{I}$ está compuesta por información de borde producto de procesos exógenos que afectan al sistema desde su frontera (e.g. afluencia desde aguas arriba, afluencia lateral). 

Las descargas $Q_i(j)$ se asocian linealmente a los estados del sistema $X_i(j)$ y la información de borde $I_i(j)$, por lo que constituye un modelo lineal simple de tránsito distribuido en un tramo, con ingreso $I_1(j)$ y descarga $Q(j)=Q_N(j)$. Así, nótese que para el caso más simple con $N=1$, el tránsito es agregado y se reduce a:

$$Q(j)=(1-K)[X(j)+I(j)]$$
 
En donde $I(j)$ es el hidrograma de aportes al tramo (afluencia por sección de entrada más los aportes laterales) y $Q(j)$ es el hidrograma transitado, en la sección ubicada aguas abajo. 

Si $N>1$, y más aun si $K_i \neq K_u$ con $u \neq i$ ($u=1,2,...,N$), el tránsito será distribuido. Por otro lado, en caso que  $K_i=K_u=K$ con $\forall u \neq i$ ($u=1,2,...,N$) el tránsito será semejante a un tránsito agregado utilizando una aproximación de $N$ reservorios lineales con constante de residencia $K$ para la formulación de una función de distribución y la ejecución de la posterior convolución. Este modelo se desarrolla posteriormente. 

Procedimentalmente, se debe especificar el conjunto de parámetros $\mathbf{K}$, para construir la matriz de transformación $\mathbf{A}$, siguiendo el orden topológico propuesto. Consecuentemente, esta matriz resume la información topológica del sistema (indica la tasa de intercambio desde el sub-tramo 'columna' al sub-tramo 'fila', luego cada _entrada_ distinta de cero indica conexión entre _nodos_). Luego, dado un vector inicial $\mathbf{X}(j)$ y un vector de condiciones iniciales $\mathbf{I}(j)$, la operación es directa y se obtiene $\mathbf{X}(j+1)$. A partir de la ecuación $Q_i(j)=(1-K_i()[X_i(j)+I_i(j)]$ se obtiene la descarga de cada sub-tramo y, de ahí, la descarga aguas abajo, mediante $Q(j)=Q_N(j)$.    

Por último, bajo las condiciones impuestas podrá mostrarse que el vector $\mathbf{K}$ es efectivamente el vector de autovalores del sistema $\mathbf{A}\mathbf{X}=\mathbf{X}$ (aspecto que puede entenderse intuitivamente, asumiendo que cada $K_i$ es la constante de _extinción_ del almacenamiento en cada sub-tramo).

## Modelo lineal para una red de drenaje compuesta por N tramos (paralelos o en serie)

Es evidente que mediante las mismas asunciones el modelo lineal puede extenderse a una topología de red más compleja, por ejemplo compuesta por tramos en serie y paralelos. Esto modificaría la matriz $\mathbf{A}$, puesto que es la que resume la información topológica (conectividad entre tramos y tasas de intercambio). Diremos que los tramos que confluyen para dar lugar a un nuevo tramo, serán _paralelos_. Por otro lado, al igual que en el caso anterior, los tramos también podrán encontrase en serie. La figura 2 muestra el caso más simple de 2 tramos paralelos (2 cursos de orden 1) que confluyen en una sección (la confluencia se denota mediante el símbolo de suma en nodo).

![Topología de red fluvial con 2 tramos paralelos. Representación de la red como sistema de reservorios lineales. $K$ es la proporción de agua libre que permanece en cada tramo y $(1-K)$ la proporción que constituye descarga aguas abajo](Reservorios_Paralelos_min.png)


Ciertamente, siguiendo el mismo razonamiento que en el caso precedente y de acuerdo a la información provista por la Fig. 2, podrá notarse que la matriz de transición $\mathbf{A}$ del modelo dinámico, para este caso, es:

$$\mathbf{A}=
\scriptsize
 \begin{bmatrix}
    K_1 & 0  \\
    0 & K_2   \\ 
\end{bmatrix}$$

y, ademaś:

$$Q(j)=Q_1(j)+Q_2(j)=(1-K_1)[X_1(j)+I_1(j)]+(1-K_2)[X_2(j)+I_2(j)]$$


Luego, teniendo en cuenta esto si consideramos el caso simple de 2 tramos _paralelos_ que descargan sobre un último tramo en _serie_ (2 cursos de orden 1 y 1 curso de orden 2, Fig. 3), la matriz de transición quedará definida por:

y, además:

$$\mathbf{A}=
\scriptsize
 \begin{bmatrix}
    K_1 & 0 & 0  \\
    0 & K_2 & 0  \\
    (1-K_1) & (1-K_2) & K_3  \\ 
\end{bmatrix}$$

$$Q(j)=Q_3(j)=(1-K_3)[X_3(j)+I_3(j)]$$

>
>

![Topología de red fluvial con 2 tramos paralelos. Representación de la red como sistema de reservorios lineales. $K$ es la proporción de agua libre que permanece en cada tramo y $(1-K)$ la proporción que constituye descarga aguas abajo](Reservorios_Paralelos_Serie.png)

 Si continuamos desarrollando nuestro razonamiento y consideramos un caso más extenso que involucre tanto tramos _paralelos_ como en _serie_, podríamos considerar una pequeña _red dendítrica_ . Esto se ejemplifica en la Fig. 4, para un sistema con $N=9$ tramos (4 cursos de orden 1, 2 cursos de orden 2, 2 cursos de orden 3 y 1 curso de orden 4).

>
>

 ![Topología de red fluvial con $N=9$ tramos paralelos y en serie. Representación de la red como sistema de reservorios lineales. $K$ es la proporción de agua libre que permanece en cada tramo y $(1-K)$ la proporción que constituye descarga aguas abajo](ModeloLinealRed.png)

 Así, al desarrollar el modelo dinámico:

$$\mathbf{A}[\mathbf{x}(j)+\mathbf{I}(j)]=\mathbf{x}(j+1)$$

de acuerdo a las conexiones y tasas de intercambio descritas en el esquema de sistema propuesto en la Fig. 4, la matriz de transformación $\mathbf{A}$ queda definida mediante:



$$\mathbf{A}=
\scriptsize
 \begin{bmatrix}
    K_1 & 0 & 0 & 0& 0& 0 & 0 & 0 & 0 \\
    0 & K_2 & 0 & 0& 0& 0 & 0 & 0 & 0 \\
    0 & 0 & K_3 & 0& 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & K_4 & 0 & 0 & 0 & 0 & 0  \\
    (1 - K_1) & (1-K_2) & 0 & 0 & K_5 & 0 & 0 & 0 & 0  \\
    0 & 0 & (1-K_3) & (1-K_4) & 0 & K_6 & 0 & 0 & 0  \\
    0 & 0 & 0 & 0 & (1-K_5) & 0 & K_7 & 0 & 0  \\
    0 & 0 & 0 & 0 & 0 & (1-K_6) & 0 & K_8 & 0 \\
    0 & 0 & 0 & 0 & 0 & 0 & (1-K_7) & (1-K_8) & K_9  \\
\end{bmatrix}$$

Así, nuevamente de acuerdo a la información topológica del sistema se obtendrá $\mathbf{A}$. Consecuentemente, a partir de un vector de condiciones iniciales $\mathbf{X}(j)$ y sobre la base de las condiciones de borde $\mathbf{I}(j)$ la operación es directa, obteniéndose $\mathbf{X}(j+1)$. 

En suma, en principio si se asumen tasas de intercambio constantes, siempre podrá desarrollarse una aplicación lineal  incluyendo la información topológica en la matriz de transformación $\mathbf{A}$, más allá del tamaño o complejidad de la red, mediante un razonamiento semejante al realizado para el caso más simple y generalizable de un tramo mediante sub-tramos en serie. Asimismo, como en apra el caso de los tramos en serie, el vector $\mathbf{K}$ será el conjutno de autovalores del sistema lineal $\mathbf{AX=X}$.  

## Modelo de tránsito en un tramo por convolución lineal 

Modelo de tránsito lineal mediante convolución de $N$ hidrogramas de entrada a un tramo definido entre un borde superior y un borde inferior. 
