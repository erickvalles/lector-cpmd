//
// Created by Erick on 15/06/2020.
//

#ifndef CPF_UTILS_H
#define CPF_UTILS_H
#include <vector>
#include "string"
#include "regex"
#include "../Atomo.h"
using std::vector;
using std::string;

/*
 * Este método calcula las posiciones periódicas para un átomo, recibe el tamaño de caja y el mitad del tamaño de la caja
 * Se pasa un arreglo de dobles por referencia ya que esas posiciones periódicas se utilizarán.*/
void calculaPosicionesPeriodicas(double *posPeriodicas,const double rx, const double ry, const double rz, double boxSize, double halfBox);
/*
 * Este mtodo */
void histograma(vector<Atomo> atomos,int n_atomos, double delta, vector<double> *histo, double boxSize, double mitadCaja);
//void histogramaMejorado(vector<Atomo> atomos,int n_atomos, double delta, vector<double> *histo, double boxSize, double mitadCaja);
double evaluaCajaRecursivo(double posicion, double boxSize, double halfBox);

double calculaDistancia(double *distancias);

double verificaComponenteR(double distancia, double boxSize, double halfBox);
double verificaComponenteParaImagen(double dif_componentes, double boxSize, double halfBox);


vector<double> componentes_restados(Atomo a1, Atomo a2);
void calcula_dist_componentes(double *distancias, Atomo a1, Atomo a2, double boxSize, double halfBox);

double calculaDelta(double boxSize, int tamHistograma);
void normalizarHistograma(vector<double> *histo, vector<double> *histo_norm, int n_atomos, int tam_histograma, double delta, double boxSize, int vecesHistograma);

double fr(double gr, double r);
double integral(vector <double> funcion_integrar, int n_divisiones);
void numerosCoordinacion(vector<double> *gdr, vector<double> *coordinaciones, double delta, double rho);

double integral_simple(double a, double b);

void factorEstructuraE(vector<double> *gdr, int tamHistograma, double delta_k, double rho, string dirSalida, double delta);
void sk(double *gdr, int tamHistograma, double delta_k, double rho, string dirSalida, double delta);
void distribucionAngulos(int numeroAtomos, std::string carpetaSalida, std::string file);
void listaVecinos(vector<Atomo> atomos, int n_atomos, float r_min, double mitadCaja, double boxSize,std::map<Atomo,vector<Atomo>> vecinos,vector<double> *histAngulos, double deltaAng, std::string trayectoria);
std::vector<Atomo> ordenarPorDistancia(std::vector<Atomo> vecinos);
#endif //CPF_UTILS_H