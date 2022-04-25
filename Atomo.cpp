//
// Created by Erick on 08/06/2020.
//


#include "Atomo.h"



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

