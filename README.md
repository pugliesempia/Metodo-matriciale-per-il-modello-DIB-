# Descrizione sintetica di quanto implementato per il modello DIB 

Questo script implementa un metodo numerico per esplorare il comportamento del modello DIB nel tempo e nello spazio.

Inizialmente, vengono definiti un set di parametri (denominato 'par') e due funzioni (F1 e F2) che descrivono le equazioni di Turing specifiche per la dinamica del modello DIB.
Dopo la discretizzazione del dominio spaziale (x, y) e con la conseguente costruzione della griglia spaziale, vengono creati operatori come la matrice tridiagonale T
necessari per il modello.
Successivamente, viene eseguita la discretizzazione temporale per consentire l'iterazione nel tempo al fine di risolvere le equazioni di Sylvester.

Per visualizzare i risultati, il programma genera un grafico bidimensionale che mostra la distribuzione spaziale di U e V. 
 Inoltre, vengono visualizzati due subplot: uno per rappresentare il valore medio spaziale di U nel tempo (un indicatore del suo comportamento medio su tutto il dominio) 
 e l'altro per mostrare l'incremento in norma di Frobenius (utilizzato per valutare la convergenza della soluzione nel tempo).
