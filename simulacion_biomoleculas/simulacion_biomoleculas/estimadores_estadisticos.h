//
//  estimadores_estadisticos.h
//  simulacion_biomoleculas
//
//  Created by Juan Román Roche on 09/10/2018.
//  Copyright © 2018 Grupo_5. All rights reserved.
//

#ifndef estimadores_estadisticos_h
#define estimadores_estadisticos_h

void estimadores_estadisticos(double * datos, int Numero_datos, double * media, double * varianza)
{
    int i;//voy a utilizar las variables media y varianza a lo largo del algoritmo para almacenar valores intermedios
    
    *media=0;
    *varianza=0;
    
    for (i=0; i<Numero_datos; i++)                              //****** nucleo ******
    {
        *media+=datos[i];
        *varianza+=datos[i]*datos[i];
    }
    
    *media=(*media)/Numero_datos;
    *varianza=(Numero_datos)/(Numero_datos-1)*(((*varianza)/Numero_datos)-(*media)*(*media));     //*** fin del nucleo ***
    
    if ((*varianza)<0) //reviso posibles errores de redondeo y los corrijo avisando al usuario
    {
        *varianza=0;
        printf("\nWARNING: la varianza salia ligeramente negativa, se ha entendido como un error de redondeo y se ha igualado a cero\n");
    }
}


#endif /* estimadores_estadisticos_h */
