//
// Created by Erick on 08/06/2020.
//


#include "Atomo.h"
#include "math.h"
#include <map>
#include <algorithm>


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

int Atomo::getEspecie() const{
    return this->especie;
}

void Atomo::setId(int id){
    this->id = id;
}

int Atomo::getId() const{
    return this->id;
}
double Atomo::dotProduct(Atomo &b){
    return this->prx*b.prx + this->pry*b.pry + this->prz*b.prz;
}

double Atomo::absolute(){
    return sqrt(this->prx*this->prx + this->pry*this->pry + this->prz*this->prz);
}

vector <double> Atomo::vectorDifference(Atomo &b){
    vector <double> result;
    result.push_back(b.prx - this->prx );
    result.push_back(b.pry - this->pry);
    result.push_back(b.prz - this->prz);
    return result;
}


double Atomo::distancia(double x, double y, double z){
    double dif_x = x-this->rx;
    double dif_y = y-this->ry;
    double dif_z = z-this->rz;
    double suma_cuadrados = (dif_x*dif_x)+(dif_y*dif_y)+(dif_z*dif_z);
    double dist = sqrt(suma_cuadrados);
    return dist;
}
/*
std::vector<Atomo> ordenarPorDistancia(std::vector<Atomo> vecinos){
    vector <Atomo> atomos_copy = vecinos;
    std::sort(atomos_copy.begin(),atomos_copy.end());
    return atomos_copy;
}*/

