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

# define MAX_INT (double) (pow(2, 31) - 1) //el -1 es por convención, intervalo abierto

//he definido unas constantes pero tengo infinitas dudas con esto
#define h 0.1
#define k 1
#define m 1
#define T 1
#define nu 10
#define kb 1
#define c0 2 * nu * kb * T
#define x0 100
#define v0 0
#define n_steps 1000

//Parámetros del Runge Kutta
#define A1 0.5
#define A2 0.5
#define beta 1
#define lambda0 1
#define lambda1 1 // también se puede usar lambda1=0 y lambda2=1
#define lambda2 0

//genera numeros aleatorios en dist plana [0, 1)
float rdm(void){
    float r = rand() / MAX_INT;
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
float Euler_maru(float t_prev, float x_prev){
    float g1, g2;

    gauss(&g1, &g2);

    return x_prev + h * force(x_prev, t_prev) + sqrtf(c0 * h) * g1;
}

//Integración por Runge-Kutta(2o orden)
float Runge_kutta2(float t_prev, float x_prev){
    float g1, g2, Z1, Z2;

    gauss(&Z1, &Z2);

    g1=force(x_prev + sqrtf(c0 * h) * lambda1 * Z1, t_prev);
    g2=force(x_prev + beta * h * g1 + sqrtf(c0 * h) * lambda2 * Z1, t_prev);

    return x_prev + h * ( A1 * g1 + A2 * g2 ) + sqrtf(c0 * h) * lambda0 * Z1;
}

//Integración por Verlet Explicito
void Verlet_exp(float t_prev, float x_prev, float v_prev, float * t_next, float * x_next, float * v_next){
    float g1, g2, a, b;
    * t_next = t_prev + h;
    gauss(&g1, &g2);

    a = (1 - 0.5 * nu * h / m ) / (1 + 0.5 * nu * h / m) ;
    b = 1.0 / (1 + 0.5 * nu * h / m) ; // no estoy seguro si el .0 hace falta

    * x_next = x_prev + b * h * v_prev + 0.5 * b * h * h * force(x_prev, t_prev) / m + 0.5 * b * h * g1 / m ;
    * v_next = a * v_prev + 0.5 * h * ( a * force(x_prev, t_prev) + force(* x_next, * t_next)) / m +  b * g1 / m ; // no tengo claro si el número aleatorio tiene que ser el mismo o no

}


void save_trajectory(float * t, float * x, float * v){
    FILE *f;
    int i;

    f = fopen("trajectoryVerlet_10.out", "w");

    /// Guardar trayectoria para Euler-Maruyama y RK2
    /*
    for (i = 0; i < n_steps; i++){
        fprintf(f, "%f %f\n", t[i], x[i]);
    }
    */
    /// Guardar trayectoria (con velocidades) para Verlet explicito

    for (i = 0; i < n_steps; i++){
        fprintf(f, "%f %f %f\n", t[i], x[i], v[i]);
    }


    fclose(f);
}


int main(int argc, const char * argv[]) {
    srand((unsigned int) time(NULL));
    int i = 0;
    float x[n_steps];
    float t[n_steps];
    float v[n_steps];

    t[0] = 0;
    x[0] = x0;
    v[0] = v0;

    /// Calcula la trayectoria de Euler-Maruyama o RK2
    /*
    for (i = 1; i < n_steps; i++){
        x[i] = Runge_kutta2(t[i - 1], x[i - 1]);
        t[i] = t[i - 1] + h;
    }
    */
    /// Calcula la trayectoria (con velocidades) para Verlet explicito
    
    for (i = 1; i < n_steps; i++){
        Verlet_exp(t[i-1], x[i-1], v[i-1], &t[i], &x[i], &v[i]);
    }


    save_trajectory(t, x, v);
    return 0;
}
