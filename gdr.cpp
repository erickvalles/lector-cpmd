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
    
    obtenerArgumentos(argc, argv, boxSize, numeroAtomos, tamHistograma, carpetaSalida, file, opc);

    cout << "el directorio de salida es: "<<carpetaSalida << endl;
    cout << "el file a leer es: "<< file << endl;

    // Se calcula la mitad de la caja aquí mismo por razones de rendimiento. 
    double mitadCaja = boxSize/2;
    
    
    // Este factor sirve para convertir de unidades atómicas a Angstroms
    const double factor_distancias =  0.52917720859;
    // Este factor sirve para convertir de unidades átomicas/us a Angstroms/us
    const double factor_velocidades = 21876.912541;
    //
    double rho = numeroAtomos/(boxSize*boxSize*boxSize);
    cout << "el valor de rho:" << rho << endl;
     //El archivo del que leeremos las trayectorias
    std::fstream input_file;

    input_file.open(file, std::ios::in);

    if(!input_file){
        std::cerr << "No se pudo abrir el archivo" << std::endl;
        return 0;
    }
    //Todos los array:
    int timestamp;
    double px, py, pz, vx, vy, vz;
    
    //Variables para almacenar los resultados
    vector<double> gofr(tamHistograma, 0.0);
    vector<double> coordinaciones(tamHistograma,0.0);
    vector<double> hist(tamHistograma,0.0);
    vector<Atomo> atomos;
    
    int n_atomos = 0;
    int trayectorias = 0;
    int veces_histograma = 0;

    /*
     * Delta es una cantidad que nos sirve para "clasificar" en qué parte del histograma estará un átomo (basado en las distancias)
     * */
    double delta = calculaDelta(boxSize,tamHistograma);
    cout << "el valor de delta:" << delta << endl;
     
    // Se crea una variable que leerá todo un renglón del archivo
    std::string trayectoria;
    //aquí, creo un archivo para ver qué posiciones periódicas son las que se obtuvieron (para fines de pruebas)
    //std::ofstream out_posiciones;

    //out_posiciones.open(carpetaSalida+"periodic_pos.txt");
    
    while (std::getline(input_file, trayectoria)){

        std::istringstream iss(trayectoria);
        //Creamos una instancia de la clase átomo
        Atomo atomo;//eliminarlos una vez que los saque del arreglo**
        
        /*
         * A la instancia de átomo le asignamos las posiciones que se obtienen directamente de la trayectoria
         * Utilizamos stod para convertir de string a double
         * Convertimos las posiciones a angstroms al asignarlas
         * */
        if(iss >> timestamp >> px >> py >> pz >> vx >> vy >> vz){
            atomo.setPos(px*factor_distancias,py*factor_distancias,pz*factor_distancias);
        }else{
            std::cerr << "Error al leer la línea: " << trayectoria << std::endl;
        }
        
        /*
         * A la instancia de átomo le asignamos las posiciones que se obtienen directamente de la trayectoria
         * Utilizamos stod para convertir de string a double
         * Convertimos las velocidades */
        //atomo.setSpeeds(stod(partes[5])*factor_velocidades,stod(partes[6])*factor_velocidades,stod(partes[7])*factor_velocidades);

        //La propiedad especie se asignará desde la interfaz gráfica también
        //Sirve para identificar qué tipo de átomo estamos utilizando
        atomo.setEspecie(1);
        //Se crea un arreglo de doubles que almacenará las posiciones periódicas para el átomo
        //double periodics[3] = {0.0,0.0,0.0};
        /*
         * Se llama al método que calcula las posiciones periódicas, se le pasan las posiciones en x,y,z, el tamaño de la caja
         * y la mitad de la caja para calcularlas
         * */
        // cout << "estoy en la trayectoria " << trayectorias << " positions_:  x"<< atomo.getrx()  << "y:"<< atomo.getry() << "z:" <<atomo.getrz() << endl;
        //calculaPosicionesPeriodicas(periodics, atomo.getrx(), atomo.getry(), atomo.getrz(), boxSize,mitadCaja);
        
        //out_posiciones << std::setprecision(20) << "         "<< periodics[0] << "        " << periodics[1] <<"     " << periodics[2]<<endl;
        //Posiciones periódicas ok!!
        
        /*
         * Como se pasó el arreglo periodics por referencia, podemos obtener las posiciones perióicas del mismo arreglo que mandamos
         * */
         
         
        //atomo.setPeriodics(periodics[0],periodics[1],periodics[2]);

        /*
         * Una vez que tenemos el átomo con las pocisiones periódicas*/
        atomos.push_back(atomo);
        n_atomos++;
        if(n_atomos==numeroAtomos){
          // cout << "ya hay " << numeroAtomos << "hay que llamar a la funcion  \n" << endl;
           histograma(atomos,numeroAtomos,delta,hist,boxSize,mitadCaja);
           //cout << "terminó con el histograma "<< veces_histograma << endl;
            n_atomos = 0;
            veces_histograma++;
            atomos.clear();
            //cout << "ya limpió el vector" << endl;
        }
        trayectorias++;
         //cout << "otra trayectoria:" << trayectorias << endl;
        


        //atomo.enviarMensaje();
       // free(line);
    }
    //out_posiciones.close();
    input_file.close();



//    std::ofstream output_file("./histograma.txt");
//    std::ostream_iterator<std::string> output_it(output_file,"\n");
//    std::copy(hist->begin(),hist->end(), output_it);
    //delete atomos;
    
    std::ofstream outfile1;
    
    outfile1.open(carpetaSalida+"/"+"histograma.txt");
    
    
    //cout << "Histograma: " << endl;
    for(int contador = 0;contador<tamHistograma;contador++){
        double r = (contador-0.5)*(delta);//es la distancia radial, máximo boxSize/2

       // cout << contador << r << val<< endl;

        //cout << val/veces_histograma << endl;//esto en otro archivo (normalizado)
        contador++;
        outfile1 << contador << "         "<< r << "        " << hist[contador] <<"     "<<endl;
    }
    outfile1.close();

    cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
    cout << "hay " <<trayectorias << " trayectorias" << endl;
    
    std::ofstream outfile2;
    
    

    
    outfile2.open(carpetaSalida+"/"+"histograma_normalizado.txt");    
    
    double hist_norm[tamHistograma];//almacenar el histograma normalizado
    //cout << "Histograma normalizado: " << endl;
    for(int contador = 0;contador<tamHistograma;contador++){
        double r = (contador-0.5)*(delta);//es la distancia radial, máximo boxSize/2



        double norm = hist[contador]/veces_histograma;//esto en otro archivo (normalizado)
        hist_norm[contador]=norm;
        //cout << contador2 << "         "<< r << "        " << norm <<"     "<<endl;
        //printf("%i       %2.9f             %f    \n",contador2,r,norm);
        
        outfile2 << contador<<  std::setprecision(20) << "         "<< r << "        " << norm <<"     "<<endl;

    }
    outfile2.close();

    
    /***********************************
        *     Gd(r)
        * 
        * 
        *************************************/


    normalizarHistograma(hist,gofr,numeroAtomos,tamHistograma,delta,boxSize,veces_histograma);
    //int contador3 = 0;
    std::ofstream outfile3;

    string outgdr = carpetaSalida+"/"+"gdr.txt";
    outfile3.open(outgdr);    
    
    
    //cout << "Histograma normalizado: " << endl;
    for(int contador=0;contador<tamHistograma;contador++){
        double r = (contador-0.5)*(delta);//es la distancia radial, máximo boxSize/2
         //cout << contador << "         "<< r << "        " << gofr[contador] <<"     "<<endl;
        
        outfile3 << contador <<  std::setprecision(20)<< "         "<< r << "        " << gofr[contador] <<"     "<<endl;

    }
    cout << "gdr:"<<outgdr<<endl;
    outfile3.close();


    /***********************
    -----------Coordinación promedio----------------
    
    **********************/
    
    
    
    
    
    
    numerosCoordinacion(gofr,coordinaciones,delta,rho);
    
    
    std::ofstream salidaCoordinacion;
    string outCoordinacion = carpetaSalida+"/"+"zdr.txt";
    salidaCoordinacion.open(outCoordinacion);

    
    
    
    for(int contador=0;contador<tamHistograma; contador++){
        double r = (contador-0.5)*(delta);//es la distancia radial, máximo boxSize/2
         //cout << contador3 << "         "<< r << "        " << norm <<"     "<<endl;
        salidaCoordinacion << std::setprecision(20)<< r << "        " << coordinaciones[contador] <<"     "<<endl;

    }

    cout << "coordinacion:" << outCoordinacion << endl;
    
    
    
    salidaCoordinacion.close();
    /*
     * Establecemos delta k 
     * vale 2*pi/(delta_r*tamaño del vector)
     * 
     * */
    double r_max = delta*16*tamHistograma;
    cout << "r max es: "<<r_max << endl;
     
    double delta_k = (2*M_PI)/r_max;
    
    std::cout << "el k_minimo es: "<< delta_k << endl;
    
    factorEstructuraE(gofr,tamHistograma,delta_k,rho,carpetaSalida,delta);
    //sk(gdr,tamHistograma,delta_k,rho,carpetaSalida,delta);

    std::cout << "Tengo el archivo" << std::endl;



    //input_file.seekg(0);





    input_file.close();
    return 1;

}