Controllo Robusto di un Elicottero 2-DOF

Questo repository contiene il codice sorgente, i modelli di simulazione MATLAB/Simulink e la documentazione sviluppati per il progetto del corso di Controllo Robusto del corso di Laurea Magistrale in Ingegneria Robotica e dell'Automazione dell'Università di Pisa (A.A. 2025/2026).

Panoramica del Progetto: L'obiettivo del progetto è lo studio, la progettazione e la validazione di diverse architetture di controllo per un elicottero a due gradi di libertà (2-DOF), caratterizzato da dinamiche non lineari e accoppiamenti geometrici e aerodinamici tra gli assi di pitch e yaw.

Il sistema è modellato attraverso una rappresentazione non lineare completa, comprensiva delle dinamiche degli attuatori e di un modello sensoriale basato su IMU. A partire dalla linearizzazione attorno a un punto di equilibrio nominale, vengono sviluppate e confrontate diverse strategie di controllo, con particolare attenzione alle prestazioni, alla reiezione dei disturbi e alla robustezza rispetto alle incertezze parametriche e dinamiche del sistema.

Il progetto si articola nelle seguenti attività:

Modellazione e linearizzazione della dinamica non lineare dell'elicottero, includendo attuatori, sensori, rumori di misura e disturbi esterni.

Progettazione e confronto di controllori LQG, con diverse architetture di feedback, precompensazione e azione integrale, e analisi della relazione con il controllo H2.

Sviluppo di controllori robusti H∞ mediante tecniche di mixed-sensitivity, utilizzando le sintesi mixsyn, hinfsyn e hinfstruct.

Implementazione della µ-sintesi attraverso la D-K Iteration, modellando esplicitamente le incertezze parametriche e dinamiche degli attuatori.

Analisi della Regione di Attrazione (RAS) dei diversi controllori sul modello non lineare, per confrontarne la capacità di riportare il sistema al punto di equilibrio a partire da condizioni iniziali perturbate.

Validazione Monte Carlo sul modello incerto, per verificare il mantenimento delle prestazioni di tracking, della reiezione dei disturbi e dei limiti operativi degli attuatori.
