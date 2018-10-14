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
#include "parisi_rapuano.h"
#include "estimadores_estadisticos.h"

#define h 0.0001
#define k 1.0
#define m 1.0
#define T 1.0
#define nu 0.0
#define kb 1.0
#define c0 (2.0 * nu * kb * T)
#define x0 1.0
#define v0 0.0
#define n_steps 1000000

//Parámetros del Runge Kutta
#define A1 0.5
#define A2 0.5
#define beta 1.0
#define lambda0 1.0
#define lambda1 1.0 // también se puede usar lambda1=0 y lambda2=1
#define lambda2 0.0

void gauss(double* g1, double* g2);
double force(double x, double t);
void Euler_maru(double t_prev, double x_prev, double v_prev, double* t_next, double* x_next, double* v_next);
double Runge_kutta2(double t_prev, double x_prev);
void Verlet_exp(double t_prev, double x_prev, double v_prev, double* t_next, double* x_next, double* v_next);
void save_trajectory(double* t, double* x, double* v, double* E_pot, double* E_kin);


int main(int argc, const char* argv[]) {
    semilla_parisi_rapuano(0);
    int i = 0;
    double* x = malloc(n_steps * sizeof(double));
    double* t = malloc(n_steps * sizeof(double));
    double* v = malloc(n_steps * sizeof(double));
    double* E_pot = malloc(n_steps * sizeof(double));
    double* E_kin = malloc(n_steps * sizeof(double));

    t[0] = 0.0;
    x[0] = x0;
    v[0] = v0;

    for (i = 1; i < n_steps; i++){
        Euler_maru(t[i-1], x[i-1], v[i-1], &t[i], &x[i], &v[i]);
    }

    for (i = 0; i < n_steps; i++){
        E_pot[i] = 0.5 * k * x[i] * x[i];
        E_kin[i] = 0.5 * m * v[i] * v[i];
    }
    
    save_trajectory(t, x, v, E_pot, E_kin);

    free(x);
    free(v);
    free(t);
    free(E_pot);
    free(E_kin);
    
    return 0;
}

//algoritmo de box-muller, si vemos que no usamos los dos numeros podemos quitar uno y así optimizamos y lo sacamos con un return
void gauss(double* g1, double* g2){
    double root = sqrt(-2.0 * log(rand_parisi_rapuano()));
    double arg = 2.0 * M_PI * rand_parisi_rapuano();

    *g1 = -root * cos(arg);
    *g2 = -root * sin(arg);
}

//funcion fuerza, está aparte para poder irla cambiando
double force(double x, double t){
    return -k * x;
}

//integración por euler maruyama
void Euler_maru(double t_prev, double x_prev, double v_prev, double* t_next, double* x_next, double* v_next){
    double g1, g2;
    gauss(&g1, &g2);
    
    *t_next = t_prev + h;
    *x_next = x_prev + h * v_prev;
    *v_next = v_prev + (h * (-nu * v_prev + force(x_prev, t_prev)) + sqrt(c0 * h) * g1) / m;
}

//Integración por Runge-Kutta(2o orden)
double Runge_kutta2(double t_prev, double x_prev){
    double g1, g2, Z1, Z2;

    gauss(&Z1, &Z2);
    
    g1 = force(x_prev + sqrt(c0 * h) * lambda1 * Z1, t_prev);
    g2 = force(x_prev + beta * h * g1 + sqrt(c0 * h) * lambda2 * Z1, t_prev);
    
    return x_prev + h * ( A1 * g1 + A2 * g2 ) + sqrt(c0 * h) * lambda0 * Z1;
}

//Integración por Verlet Explicito
void Verlet_exp(double t_prev, double x_prev, double v_prev, double* t_next, double* x_next, double* v_next){
    double g1, g2, a, b;
    *t_next = t_prev + h;

    gauss(&g1, &g2);

    a = (1.0 - 0.5 * nu * h / m ) / (1.0 + 0.5 * nu * h / m) ;
    b = 1.0 / (1.0 + 0.5 * nu * h / m) ;
    
    *x_next = x_prev + b * h * v_prev + 0.5 * b * h * h * force(x_prev, t_prev) / m + 0.5 * b * sqrt(h * c0) * g1 / m ;
    *v_next = a * v_prev + 0.5 * h * ( a * force(x_prev, t_prev) + force(*x_next, *t_next)) / m + b * sqrt(h * c0) * g1 / m ;
}


void save_trajectory(double* t, double* x, double* v, double* E_pot, double* E_kin){
    FILE* f;
    int i;

    f = fopen("prueba.out", "w");

    for (i = 0; i < n_steps; i++){
        fprintf(f, "%lf %lf %lf %lf %lf %lf\n", t[i], x[i], v[i], E_pot[i], E_kin[i], E_kin[i] + E_pot[i]);
    }

    fclose(f);
}

void save_trajectory_Euler_RK(double* t, double* x){
    FILE* f;
    int i;

    f = fopen("trajectoryRK_nu10.0_h0.01.out", "w");

    for (i = 0; i < n_steps; i++){
        fprintf(f, "%lf %lf\n", t[i], x[i]);
    }

    fclose(f);
}
