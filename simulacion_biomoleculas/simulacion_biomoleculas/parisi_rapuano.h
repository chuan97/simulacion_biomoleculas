//
//  parisi_rapuano.h
//  simulacion_biomoleculas
//
//  Created by Juan Román Roche on 08/10/2018.
//  Copyright © 2018 Grupo_5. All rights reserved.
//

#ifndef parisi_rapuano_h
#define parisi_rapuano_h
#define NormRANu (2.3283063671E-10F)
//este generador de numeros aleatorios genera un double en [0, 1)
//para ello se debe inicializar ejecutando la funcion semilla_parisi_rapuano
//una vez inicializado se puede llama a la funcion rand_parisi_rapuano cada vez que se desee un nuevo numero aleatorio

unsigned int irr[256], ir1; //funciona, pero falta hacer un test y ver que genera un histograma plano
unsigned char ind_ran, ig1, ig2, ig3;

void semilla_parisi_rapuano(int seed)
{
    int i;
    
    if (seed==0)
        srand((unsigned int) time(NULL));
    else
        srand(seed);
    
    for (i=0; i<256; i++)
        irr[i]=(rand()<<16)+rand();
    
    ind_ran=ig1=ig2=ig3=0;
}

double rand_parisi_rapuano(void)
{
    ig1=ind_ran-24;
    ig2=ind_ran-55;
    ig3=ind_ran-61;
    irr[ind_ran]=irr[ig1]+irr[ig2];
    ir1=(irr[ind_ran]^irr[ig3]);
    ind_ran++;
    
    return ir1*NormRANu;
}

double interval_rand_parisi_rapuano(double a, double b)
{
    return (a+(b-a)*rand_parisi_rapuano());
}


#endif /* parisi_rapuano_h */
