//
//  main.c
//  doble_pozo
//
//  Created by Juan Román Roche on 22/10/2018.
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
#define B 1.0
#define m 1.0
#define T 0.2
#define eta 0.1
#define kb 1.0
#define c0 (2.0 * eta * kb * T)
#define x0 2.0
#define v0 0.0
#define n_steps 1000000
#define n_term (n_steps / 10) //cantidad arbitraria
#define n_intervalos 100

double force(double x, double t);
void Verlet_exp(double t_prev, double x_prev, double v_prev, double* t_next, double* x_next, double* v_next);
void save_trajectory(double* t, double* x, double* v, double* E_kin, double* E_pot);
void control_parameters(double* x, double* v, double* E_kin, double* E_pot);
double gauss(void);

int main(int argc, const char* argv[]) {
    semilla_parisi_rapuano(0);
    int i = 0;
    double* x = malloc(n_steps * sizeof(double)); //declarados en el heap en vez del stack para tener m·s espacio
    double* t = malloc(n_steps * sizeof(double));
    double* v = malloc(n_steps * sizeof(double));
    double* E_pot = malloc(n_steps * sizeof(double));
    double* E_kin = malloc(n_steps * sizeof(double));
    
    //termalizaciÛn
    t[0] = 0.0;
    x[0] = x0;
    v[0] = v0;
    
    for (i = 1; i < n_term; i++){
        Verlet_exp(t[i-1], x[i-1], v[i-1], &t[i], &x[i], &v[i]);
    }
    //fin de termalizaciÛn, se toman la ˙ltima posiciÛn y velocidad como etaevos par·metros de inicio
    
    t[0] = 0.0;
    x[0] = x[n_term - 1];
    v[0] = v[n_term - 1];
    
    for (i = 1; i < n_steps; i++){
        Verlet_exp(t[i-1], x[i-1], v[i-1], &t[i], &x[i], &v[i]);
    }
    
    for (i = 0; i < n_steps; i++){
        E_pot[i] = 0.5 * B * (x[i] * x[i] - 1) * (x[i] * x[i] - 1);
        E_kin[i] = 0.5 * m * v[i] * v[i];
    }
    
    save_trajectory(t, x, v, E_kin, E_pot);
    //control_parameters(x, v, E_kin, E_pot);
    
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

//funcion fuerza del potencial de doble pozo
double force(double x, double t){
    return -2.0 * B * (x * x * x - x);
}

//IntegraciÛn por Verlet Explicito
void Verlet_exp(double t_prev, double x_prev, double v_prev, double* t_next, double* x_next, double* v_next){
    double g1, a, b;
    
    g1 = gauss();
    
    a = (1.0 - 0.5 * eta * h / m ) / (1.0 + 0.5 * eta * h / m) ;
    b = 1.0 / (1.0 + 0.5 * eta * h / m) ;
    
    *t_next = t_prev + h;
    *x_next = x_prev + b * h * v_prev + 0.5 * b * h * h * force(x_prev, t_prev) / m + 0.5 * b * h * sqrt(h * c0) * g1 / m ;
    *v_next = a * v_prev + 0.5 * h * ( a * force(x_prev, t_prev) + force(*x_next, *t_next)) / m + b * sqrt(h * c0) * g1 / m ;
}


void save_trajectory(double* t, double* x, double* v, double* E_kin, double* E_pot){
    FILE* f;
    int i;
    
    f = fopen("prueba.out", "w");
    
    int increase = n_steps / (n_steps / 1); // limita los puntos a plotear, para no saturar el plotter
    
    for (i = 0; i < n_steps; i += increase){
        fprintf(f, "%lf %lf %lf %lf %lf %lf\n", t[i], x[i], v[i], E_pot[i], E_kin[i], E_kin[i] + E_pot[i]);
    }
    
    fclose(f);
}

void control_parameters(double* x, double* v, double* E_kin, double* E_pot) {
    double mean_kin, var_kin, mean_pot, var_pot, mean_x, var_x, mean_v, var_v, delta_x, min_x, max_x, delta_v, min_v, max_v;
    int i;
    FILE* f_x;
    FILE* f_v;
    
    double Hist_x[n_intervalos], Hist_v[n_intervalos];
    
    f_x = fopen("histograma_x_runge_10.0.out", "w");
    f_v = fopen("histograma_v_runge_10.0.out", "w");
    
    estimadores_estadisticos(E_kin, n_steps, &mean_kin, &var_kin);
    estimadores_estadisticos(E_pot, n_steps, &mean_pot, &var_pot);
    estimadores_estadisticos(x, n_steps, &mean_x, &var_x);
    estimadores_estadisticos(v, n_steps, &mean_v, &var_v);
    
    histograma(x, Hist_x, n_steps, n_intervalos, &delta_x, &min_x, &max_x);
    histograma(v, Hist_v, n_steps, n_intervalos, &delta_v, &min_v, &max_v);
    
    for (i = 0; i < n_intervalos; i++) {
        fprintf(f_x, "%lf %lf\n", min_x + 0.5 * delta_x + i * delta_x, Hist_x[i]);
        fprintf(f_v, "%lf %lf\n", min_v + 0.5 * delta_v + i * delta_v, Hist_v[i]);
    }
    
    /// Los histogramas se han construido con un paso temporal h = 0.01 (al menos runge y euler, verlet en duda)
    
    fclose(f_x);
    fclose(f_v);
    
    
    printf("EnergÌa cinÈtica media: %lf\nEnergÌa potencial media: %lf\nEnergÌa tÈrmica: %lf\n", mean_kin, mean_pot, 0.5 * kb * T);
    printf("PosiciÛn media: %lf, Varianza: %lf\n", mean_x, var_x);
    printf("Velocidad media: %lf, Varianza: %lf\n", mean_v, var_v);
}

