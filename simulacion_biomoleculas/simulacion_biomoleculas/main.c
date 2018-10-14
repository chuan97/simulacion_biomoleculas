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
#include "histogram.h"

#define h 0.0001
#define k 1.0
#define m 1.0
#define T 1.0
#define nu 1.0
#define kb 1.0
#define c0 (2.0 * nu * kb * T)
#define x0 1.0
#define v0 0.0
#define n_steps 10000000
#define n_term (n_steps / 10) //cantidad arbitraria

//Parámetros del Runge Kutta
#define A1 0.5
#define A2 0.5
#define beta 1.0
#define lambda0 1.0
#define lambda1 1.0 // también se puede usar lambda1=0 y lambda2=1
#define lambda2 0.0

double force(double x, double t);
void Euler_maru(double t_prev, double x_prev, double v_prev, double* t_next, double* x_next, double* v_next);
void Runge_kutta2(double t_prev, double x_prev, double v_prev, double* t_next, double* x_next, double* v_next);
void Verlet_exp(double t_prev, double x_prev, double v_prev, double* t_next, double* x_next, double* v_next);
void save_trajectory(double* t, double* x, double* v, double* E_kin, double* E_pot);
void control_parameters(double* x, double* v, double* E_kin, double* E_pot);
double gauss(void);

int main(int argc, const char* argv[]) {
    semilla_parisi_rapuano(0);
    int i = 0;
    double* x = malloc(n_steps * sizeof(double)); //declarados en el heap en vez del stack para tener más espacio
    double* t = malloc(n_steps * sizeof(double));
    double* v = malloc(n_steps * sizeof(double));
    double* E_pot = malloc(n_steps * sizeof(double));
    double* E_kin = malloc(n_steps * sizeof(double));
    
    //termalización
    t[0] = 0.0;
    x[0] = x0;
    v[0] = v0;
    
    for (i = 1; i < n_term; i++){
        Euler_maru(t[i-1], x[i-1], v[i-1], &t[i], &x[i], &v[i]);
    }
    //fin de termalización, se toman la última posición y velocidad como nuevos parámetros de inicio
    
    t[0] = 0.0;
    x[0] = x[n_term - 1];
    v[0] = v[n_term - 1];

    for (i = 1; i < n_steps; i++){
        Verlet_exp(t[i-1], x[i-1], v[i-1], &t[i], &x[i], &v[i]);
    }

    for (i = 0; i < n_steps; i++){
        E_pot[i] = 0.5 * k * x[i] * x[i];
        E_kin[i] = 0.5 * m * v[i] * v[i];
    }
    
    save_trajectory(t, x, v, E_kin, E_pot);
    control_parameters(x, v, E_kin, E_pot);

    free(x); //liberando el heap
    free(v);
    free(t);
    free(E_pot);
    free(E_kin);
    
    return 0;
}

//algoritmo de box-muller, optimizado
double gauss(void){
    static int phase = 1;
    static double g1, g2;
    
    if (phase == 1){
        phase = 2;
        
        double root = sqrt(-2.0 * log(rand_parisi_rapuano()));
        double arg = 2.0 * M_PI * rand_parisi_rapuano();
        
        g1 = -root * cos(arg);
        g2 = -root * sin(arg);
        
        return g1;
    }
    
    else{
        phase = 1;
        
        return g2;
    }
}


//funcion fuerza, está aparte para poder irla cambiando
double force(double x, double t){
    return -k * x;
}

//integración por euler maruyama
void Euler_maru(double t_prev, double x_prev, double v_prev, double* t_next, double* x_next, double* v_next){
    *t_next = t_prev + h;
    *x_next = x_prev + h * v_prev;
    *v_next = v_prev + (h * (-nu * v_prev + force(x_prev, t_prev)) + sqrt(c0 * h) * gauss()) / m;
}

//Integración por Runge-Kutta(2o orden)
void Runge_kutta2(double t_prev, double x_prev, double v_prev, double* t_next, double* x_next, double* v_next){
    double g1, g2, Z1;
    
    Z1 = gauss();

    g1 = -nu * v_prev + force(x_prev + sqrt(c0 * h) * lambda1 * Z1, t_prev);
    g2 = -nu * v_prev + force(x_prev + beta * h * g1 + sqrt(c0 * h) * lambda2 * Z1, t_prev);
    
    *t_next = t_prev + h;
    *x_next = x_prev + h * v_prev;
    *v_next = v_prev + (h * ( A1 * g1 + A2 * g2 ) + sqrt(c0 * h) * lambda0 * Z1) / m;
}

//Integración por Verlet Explicito
void Verlet_exp(double t_prev, double x_prev, double v_prev, double* t_next, double* x_next, double* v_next){
    double g1, a, b;
   
    g1 = gauss();

    a = (1.0 - 0.5 * nu * h / m ) / (1.0 + 0.5 * nu * h / m) ;
    b = 1.0 / (1.0 + 0.5 * nu * h / m) ;
    
    *t_next = t_prev + h;
    *x_next = x_prev + b * h * v_prev + 0.5 * b * h * h * force(x_prev, t_prev) / m + 0.5 * b * h * sqrt(h * c0) * g1 / m ;
    *v_next = a * v_prev + 0.5 * h * ( a * force(x_prev, t_prev) + force(*x_next, *t_next)) / m + b * sqrt(h * c0) * g1 / m ;
}


void save_trajectory(double* t, double* x, double* v, double* E_kin, double* E_pot){
    FILE* f;
    int i;

    f = fopen("prueba.out", "w");
    
    int increase = n_steps / 100000; // limita los puntos a plotear a 100000, para no saturar el plotter

    for (i = 0; i < n_steps; i += increase){
        fprintf(f, "%lf %lf %lf %lf %lf %lf\n", t[i], x[i], v[i], E_pot[i], E_kin[i], E_kin[i] + E_pot[i]);
    }

    fclose(f);
}

void control_parameters(double* x, double* v, double* E_kin, double* E_pot){
    double mean_kin, var_kin, mean_pot, var_pot;
    
    estimadores_estadisticos(E_kin, n_steps, &mean_kin, &var_kin);
    estimadores_estadisticos(E_pot, n_steps, &mean_pot, &var_pot);
    
    printf("Energía cinética media: %lf\nEnergía potencial media: %lf\nEnergía térmica: %lf\n", mean_kin, mean_pot, 0.5 * kb * T);
}
