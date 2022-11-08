//
// Created by Erick on 08/06/2020.
//

#ifndef CPF_ATOMO_H
#define CPF_ATOMO_H
#include <string>
#include <vector>
#include <map>
/**
 * Agregar vector velocidades
 * Agregar propiedad "especie"->string <- qué tipo de átomo
 *
 * */
using std::vector;
using std::string;
class Atomo{
private:
    int id = 0;
    //double* pos;
    double rx = 0;//Angstroms
    double ry = 0;//Angstroms
    double rz = 0;//Angstroms
    //double* perPos;
    double prx=0;//Angstroms
    double pry=0;//Angstroms
    double prz=0;//Angstroms
    //double* vel;
    double vx=0;//Angstroms/picoSegs
    double vy=0;//Angstroms/picoSegs
    double vz=0;//Angstroms/picoSegs

    int especie = 0;//cambiar a int
public:
    std::string enviarMensaje();
    void setPos(double rx,double ry, double rz);
    void setPosArr(double pos[3]);
    double getrx();
    double getry();
    double getrz();
    double distanciaVecino=0;
    void setEspecie(int especie);
    int getEspecie() const;
    void setSpeeds(double vx, double vy, double vz);
    void setPeriodics(double prx, double pry, double prz);
    double getPrx();
    double getPry();
    double getPrz();
    double distancia(double x, double y, double z);
    void setId(int id);
    int getId() const;
    bool operator<(const Atomo& b) const{
        return distanciaVecino < b.distanciaVecino;
    }
    double dotProduct(Atomo &b);
    double absolute();
    //obtain a new vector from the difference between two atoms
    vector<double> vectorDifference(Atomo &b);
    //std::vector<Atomo> ordenarPorDistancia(std::vector<Atomo> vecinos);

};


#endif //CPF_ATOMO_H