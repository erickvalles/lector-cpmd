//
// Created by Erick on 08/06/2020.
//


#include "Atomo.h"
#include "math.h"


std::string Atomo::enviarMensaje(){
    return "Listo";
}

void Atomo::setEspecie(int especie) {
    this->especie = especie;
}

void Atomo::setPos(double rx,double ry, double rz){
    this->rx = rx;
    this->ry = ry;
    this->rz = rz;
}

double Atomo::getrx(){
    return this->rx;
}
double Atomo::getry(){
    return this->ry;
}
double Atomo::getrz(){
    return this->rz;
}
void Atomo::setSpeeds(double vx, double vy, double vz){
    this->vx=vx;
    this->vy=vy;
    this->vz=vz;
}
void Atomo::setPeriodics(double prx, double pry, double prz){
    this->prx = prx;
    this->pry = pry;
    this->prz = prz;
}

double Atomo::getPrx(){
    return this->prx;
}
double Atomo::getPry(){
    return this->pry;
}
double Atomo::getPrz(){
    return this->prz;
}

double Atomo::distancia(double x, double y, double z){
    double dif_x = x-this->rx;
    double dif_y = y-this->ry;
    double dif_z = z-this->rz;
    double suma_cuadrados = (dif_x*dif_x)+(dif_y*dif_y)+(dif_z*dif_z);
    double dist = sqrt(suma_cuadrados);
    return dist;
}

