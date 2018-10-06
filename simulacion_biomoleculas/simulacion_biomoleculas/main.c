//
//  main.c
//  oscilador armónico estocástico
//
//  Created by Juan Román Roche & Fernando Lorén Mastral on 02/10/2018.
//  Copyright © 2018 Grupo_5. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

# define MAX_INT (double) (pow(2, 31) - 1)

//he definido unas constantes pero tengo infinitas dudas con esto
#define h 0.001
#define m 1
#define T 1
#define mu 1
#define kb 1
#define c0 2 * mu * kb * T

//genera numeros aleatorios en dist plana [0, 1)
float rdm(void){
    float r = rand() / MAX_INT; //el -1 es por convención, intervalo abierto
    return r;
}

//algoritmo de box-muller, si vemos que no usamos los dos numeros podemos quitar uno y así optimizamos y lo sacamos con un return
void gauss(float * g1, float * g2){
    float root = sqrtf(-1 *2 * logf(rdm()));
    float arg = 2 * M_PI * rdm();
    
    *g1 = -1 * root * cosf(arg);
    *g2 = -1 * root * sinf(arg);
}

//funcion fuerza, esta aparte para poder irla cambiando
float force(float x, float t){
    return -1 * x;
}

//integración por euler maruyama
float euler_maru(float t0, float x0){
    float g1, g2;
    
    gauss(&g1, &g2);
    
    return x0 + h * force(x0, t0) + sqrtf(c0 * h) * g1;
}


int main(int argc, const char * argv[]) {
    srand((unsigned int) time(NULL));
    FILE *f;
    
    f = fopen("debug.out", "w");
    
    printf("%f", euler_maru(0.0, 0.0));
        
    
    return 0;
}
