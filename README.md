# Pydro-SIyAH
Librería de procedimientos numéricos para la generación de previsiones hidrológicas, en proceso de codificación en lenguaje Python. Objetos y métodos/funciones para la simulación del proceso de transformación de precipitación en escorrentía y del tránsito de caudales (traslado y suma) 

## instalación

    python -m venv .
    source bin/activate
    pip install -r requirements.txt


## configuración

   mkdir config
   nano config/config.yml
   # ingresar parámetros de conexión a api:
   api:
     url: host:port/path
     token: my_token
   
