#include <iostream>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <iterator>
#include "regex"
#include "Atomo.h"
#include "utils/utils.h"
#include <ctime>
#include <iomanip>
#include "cmath"
#include <map>

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::ifstream;
/*
* Encontrar la función de estructura estático final
* añadir el fscanf
* se ejecuta desde la terminal con esto:
* /usr/bin/g++ -O3 -g utils/utils.cpp utils/utils.h Atomo.cpp Atomo.h main.cpp -o main -lfftw3
*/

void obtenerArgumentos(int argc, char **argv, double &boxSize, int &numeroAtomos, int &tamHistograma, std::string &carpetaSalida, std::string &file, int &opc, float &r_min);

int dist_angulos(int numeroAtomos, std::string carpetaSalida, std::string file, double boxSize);

int main(int argc, char **argv) {//recibir como args
   
    /*
     * Los siguientes datos se pedirán desde la interfaz gráfica del programa
     * */
    double boxSize;
    int numeroAtomos;
    int tamHistograma;
    string carpetaSalida;
    string file;
    int opc;
    float r_min;
    
    obtenerArgumentos(argc, argv, boxSize, numeroAtomos, tamHistograma, carpetaSalida, file, opc, r_min);

    cout << "el directorio de salida es: "<<carpetaSalida << endl;
    cout << "el file a leer es: "<< file << endl;
    
   dist_angulos(numeroAtomos,carpetaSalida,file, boxSize);
    
    return 0;
}


int dist_angulos(int numeroAtomos, std::string carpetaSalida, std::string file, double boxSize){

    cout << "Se calcularán los ángulos" <<endl;
     // Se calcula la mitad de la caja aquí mismo por razones de rendimiento. 
    double mitadCaja = boxSize/2;
    
    
    // Este factor sirve para convertir de unidades atómicas a Angstroms
    const double factor_distancias =  0.52917720859;
    // Este factor sirve para convertir de unidades átomicas/us a Angstroms/us
    const double factor_velocidades = 21876.912541;
    //
    
     //El archivo del que leeremos las trayectorias
    
    
   
    
    vector<double> *gofr = {nullptr};

    //abriremos el archivo de la trayectoria en modo lectura
    std::fstream input_file;

    input_file.open(file, std::ios::in);

    
    //int contador;

    //Si no se logra abrir el archivo, mostramos un mensaje de error
    if(!input_file){
        std::cerr << "No se pudo abrir el archivo" << std::endl;
        return 0;
    }

    //Se utiliza un vector de objetos del tipo átomo
    vector<Atomo> *atomos = {nullptr};
    //Se reserva memoria para contener los átomos involucrados en la simulación
    atomos = new vector<Atomo>[numeroAtomos];

   
    
    // cuando esta variable llegue a n_atomos se reiniciará y el arreglo de átomos se vaciará
    int n_atoms = 0;
    //es una variable que nos sirve para saber si se están recorriendo todas las trayectorias
    int trayectorias = 0;
    // es una variable para saber cuántas veces se llama a la función histograma
    //int veces_histograma = 0;
    //crear un vector para almacenar el histograma de angulos
    
    double tamHistAngulos = 400.0;
    double deltaAng = 180.0/tamHistAngulos;
    std::cout << deltaAng << endl;
    vector<double> *histAngulos = {nullptr};
    double halfBox = boxSize/2;
    histAngulos = new vector<double>(tamHistAngulos,0.0);
    // Se crea una variable que leerá todo un renglón del archivo
    string trayectoria;
    std::map<Atomo,vector<Atomo>> vecinos;
    std::ofstream out_posiciones;
    std::ofstream out_tray_prob;

    out_posiciones.open(carpetaSalida+"periodic_pos_angulos.txt");
    out_tray_prob.open(carpetaSalida+"less_inputs.txt");
        // de lo contrario, se crea un archivo de prueba (que es más pequeño)

    int timestamp;
    double px, py, pz, vx, vy, vz;
    while (input_file >> timestamp >> px >> py >> pz >> vx >> vy >> vz){
        //se crea un vector llamado partes, ya que cada línea la dividiremos por los espacios en blanco que contiene
        
        //Creamos una instancia de la clase átomo
        Atomo atomo;//eliminarlos una vez que los saque del arreglo**
        
        /*
         * A la instancia de átomo le asignamos las posiciones que se obtienen directamente de la trayectoria
         * Utilizamos stod para convertir de string a double
         * Convertimos las posiciones a angstroms al asignarlas
         * */
        atomo.setPos(px*factor_distancias,py*factor_distancias,pz*factor_distancias);
        /*
         * A la instancia de átomo le asignamos las posiciones que se obtienen directamente de la trayectoria
         * Utilizamos stod para convertir de string a double
         * Convertimos las velocidades */
        //atomo.setSpeeds(stod(partes[5])*factor_velocidades,stod(partes[6])*factor_velocidades,stod(partes[7])*factor_velocidades);

        //La propiedad especie se asignará desde la interfaz gráfica también
        //Sirve para identificar qué tipo de átomo estamos utilizando
        atomo.setEspecie(1);
        atomo.setId(n_atoms);
        //Se crea un arreglo de doubles que almacenará las posiciones periódicas para el átomo
        double periodics[3] = {0.0,0.0,0.0};
        /* 
         * Se llama al método que calcula las posiciones periódicas, se le pasan las posiciones en x,y,z, el tamaño de la caja
         * y la mitad de la caja para calcularlas
         * */

        
        calculaPosicionesPeriodicas(periodics, atomo.getrx(), atomo.getry(), atomo.getrz(), boxSize,mitadCaja);
        
        /*
         * Como se pasó el arreglo periodics por referencia, podemos obtener las posiciones perióicas del mismo arreglo que mandamos
         * */
         
         
        atomo.setPeriodics(periodics[0],periodics[1],periodics[2]);
        out_posiciones << std::setprecision(20) << "         "<< periodics[0] << "        " << periodics[1] <<"     " << periodics[2]<<endl;
        /*
         * Una vez que tenemos el átomo con las pocisiones periódicas*/
        (*atomos).push_back(atomo);
        n_atoms++;
        if(n_atoms==numeroAtomos){
           //cout << "ya hay " << numeroAtomos << "hay que llamar a la funcion  \n" << endl;
           //aquí debo llamar a la función que calcule la lista de vecinos
           //listaVecinos(*atomos,n_atoms,3.2,mitadCaja, boxSize, vecinos, histAngulos, deltaAng, trayectoria);
           vecinosMejorada(atomos,r_min,boxSize, halfBox);
           
           (*atomos).clear();
            n_atoms = 0;

        }
        
        trayectorias++;


        
       // free(line);
    }

out_posiciones.close();
out_tray_prob.close();
    //crear un archivo para volcar el contenido del histograma
    std::ofstream salidaAngulos;
    double area_bajo_angulos = integral(*histAngulos,tamHistAngulos);
    string outAngulos = carpetaSalida+"/"+"hist_angulos.txt";
    salidaAngulos.open(outAngulos);

    for(int i=0; i<tamHistAngulos;i++){
        double angulo = deltaAng*i;
        double conteoAngulos = (*histAngulos)[i];
        //imprimir el contador con una mayor precision
        salidaAngulos << std::setprecision(10) << angulo << "                      " << conteoAngulos <<"            "<< conteoAngulos/area_bajo_angulos<< endl;
        //imprimir en consola la variable i
       // cout << (*histAngulos)[i] << endl;
        
    }

    cout << "angulos:"<<outAngulos<<endl;

    return 0;
}

