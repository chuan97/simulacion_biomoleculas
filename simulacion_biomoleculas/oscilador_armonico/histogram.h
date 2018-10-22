//
//  histogram.h
//  simulacion_biomoleculas
//
//  Created by Juan Román Roche on 06/10/2018.
//  Copyright © 2018 Grupo_5. All rights reserved.
//

#ifndef histogram_h
#define histogram_h
//calcula el histograma de una distribucion de valores

//data -> entrada: data con los que calcula el histograma
//Hist -> salida: histograma que devuelve
//N_datos -> entrada: numero de casillas en el array data
//intervalos -> entrada: numero de casillas en el array Hist
//delta -> salida: tamaño de los intervalos
//m -> salida: valor Minimo de la distribucion
//M -> salida: valor Maximo de la distribucion


void histograma (double * data, double * Hist, int N_datos, int intervalos, double * delta, double * m, double * M)
{
    int j, indice;
    double span, Min, Max;
    
    for (j=0; j<intervalos; j++)         //inicializa a cero el histograma
        Hist[j]=0;
    
    Min=Max=data[0];
    
    for (j=1; j<N_datos; j++)         //calcula el Max y Min del intervalo
    {
        if (data[j]<Min)
            Min=data[j];
        
        if (data[j]>Max)
            Max=data[j];
    }
    
    span=(Max-Min)/intervalos;
    
    if (span<=0)
    {
        printf("ERROR");
        exit(1);
    }
    
    for (j=0; j<N_datos; j++)           // nucleo del programa, asigna a cada valor de data el intervalo que le corresponde
    {                                   //despues calcula el indice correspondiente en el array Hist
        indice=((data[j]-Min)/span);    //y contabiliza ese valor
        Hist[indice]++;
    }                                   //****** fin del nucleo ******
    
    *delta=span;
    *m=Min;
    *M=Max;
    
    for (j=0; j<intervalos; j++)      //normaliza el histograma para que la integral del la densidad de probabilidad sea uno
        Hist[j]*=1/(span*N_datos);
}

#endif /* histogram_h */
