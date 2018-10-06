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
#define k 1
#define m 1
#define T 1
#define mu 0.0000001
#define kb 1
#define c0 2 * mu * kb * T
#define x0 0.01
#define n_steps 1000000

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
    return -1 * k * x;
}

//integración por euler maruyama
float euler_maru(float t_prev, float x_prev){
    float g1, g2;
    
    gauss(&g1, &g2);
    
    return x_prev + h * force(x_prev, t_prev) + sqrtf(c0 * h) * g1;
}

void save_trajectory(float * x, float * t){
    FILE *f;
    int i;
    
    f = fopen("trajectory.out", "w");
    
    for (i = 0; i < n_steps; i++){
        fprintf(f, "%f %f\n", x[i], t[i]);
    }
    
    fclose(f);
}


int main(int argc, const char * argv[]) {
    srand((unsigned int) time(NULL));
    int i = 0;
    float x[n_steps];
    float t[n_steps];
    
    t[0] = 0;
    x[0] = x0;
    
    for (i = 1; i < n_steps; i++){
        x[i] = euler_maru(t[i - 1], x[i - 1]);
        t[i] = t[i - 1] + h;
    }
        
    save_trajectory(x, t);
    return 0;
}
