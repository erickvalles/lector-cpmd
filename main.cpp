#include <iostream>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <iterator>
#include "regex"
#include "Atomo.h"
#include "utils/utils.h"
#include <ctime>
#include <boost/algorithm/string.hpp>
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
*/
int gdr_main(double boxSize, int numeroAtomos, int tamHistograma, std::string carpetaSalida, std::string file, int prueba);
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
    
    if(argc < 2){
        // Es el tamaño de la caja de simulación que se utilizará
        boxSize = 22.3797;
        // Se pide el número de átomos que se utilizaron en la dimámica molecular
        numeroAtomos = 300;
        // Se pide el tamaño del histograma que se utilizará
        tamHistograma = 512;
        carpetaSalida = "/home/erick/data/salidas/";
        //file = "/home/erick/Ge00Sb00Te100T823K.cpmd";//"/home/erick/protocolo/copia.cpmd";////"/home/erick/protocolo/TRAJECTORY_00_Te_T823K.cpmd"; 
        file = "/home/erick/reduced.cpmd";
    }else{
        boxSize = std::stod(argv[1]);
        numeroAtomos = std::stoi(argv[2]);
        tamHistograma = std::stoi(argv[3]);
        file = argv[4];
        carpetaSalida = argv[5];
    }
    cout << "el directorio de salida es: "<<carpetaSalida << endl;
    
   
    
    
    //pediremos también el directorio para las salidas

    //Esta variable nos sirve para escribir dos tipos de archivos, unos para una trayectoria recortada
    bool prueba = false;
    
    
    //Si NO estamos haciendo pruebas, quitamos el comentario a la siguiente línea
    //const char *file = "/home/erick/protocolo/Ge00Sb00Te100T823K.cpmd";
    //string file = "/home/erick/protocolo/TRAJECTORY_00_Te_T823K.cpmd";    
    
    

    
   

    //menú
    int opc;
    cout << "¿Qué deseas calcular?"<<endl;
    cout << "1.- g(r) \n2.- ADF" << endl;

    std::cin >> opc;
    if(opc == 1){
        goto gdr;
    }else{
        goto angulos;
    }
    

    //Se lee el archivo línea por línea y se almacena en la variable trayectoria
    gdr:
        gdr_main(boxSize,numeroAtomos,tamHistograma,carpetaSalida,file,0);
    angulos:
        dist_angulos(numeroAtomos,carpetaSalida,file, boxSize);
    

    return 0;
}


int gdr_main(double boxSize, int numeroAtomos, int tamHistograma, std::string carpetaSalida, std::string file, int prueba){
     
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
    ifstream input_file;
    
    //Todos los array:
    
    double hist_norm[tamHistograma];//almacenar el histograma normalizado
    
    vector<double> *gofr = {nullptr};
    gofr = new vector<double>(tamHistograma,0.0);//almacenar gor
    
    vector<double> *coordinaciones = {nullptr};

    coordinaciones = new vector<double>(tamHistograma,0.0);/// almacenar los números de coordinacion

    //abriremos el archivo de la trayectoria en modo lectura
    input_file.open(file);
    
    int contador;

    //Si no se logra abrir el archivo, mostramos un mensaje de error
    if(!input_file){
        std::cerr << "No se pudo abrir el archivo" << std::endl;
        return 0;
    }

    //Se utiliza un vector de objetos del tipo átomo
    vector<Atomo> *atomos = {nullptr};
    //Se reserva memoria para contener los átomos involucrados en la simulación
    atomos = new vector<Atomo>[numeroAtomos];

    /*
     * Delta es una cantidad que nos sirve para "clasificar" en qué parte del histograma estará un átomo (basado en las distancias)
     * */
    double delta = calculaDelta(boxSize,tamHistograma);
    cout << "el valor de delta:" << delta << endl;

    //Tenemos un vector de dobles que representará el histograma
    vector<double> *hist = {nullptr};
    //Se aloja memoria suficiente para contener el histograma y todos los elementos se inicializan en 0
    hist = new vector<double>(tamHistograma,0.0);
    
    // cuando esta variable llegue a n_atomos se reiniciará y el arreglo de átomos se vaciará
    int n_atoms = 0;
    //es una variable que nos sirve para saber si se están recorriendo todas las trayectorias
    int trayectorias = 0;
    // es una variable para saber cuántas veces se llama a la función histograma
    int veces_histograma = 0;
    
    // Se crea una variable que leerá todo un renglón del archivo
    string trayectoria;
    //aquí, creo un archivo para ver qué posiciones periódicas son las que se obtuvieron (para fines de pruebas)
    std::ofstream out_posiciones;
    //Si no es una prueba, se crea un archivo "definitivo"
    if(!prueba){
        out_posiciones.open(carpetaSalida+"periodic_pos.txt");
    }else{
        // de lo contrario, se crea un archivo de prueba (que es más pequeño)
        out_posiciones.open(carpetaSalida+"periodic_pos_prueba.txt");
    }
    while (std::getline(input_file,trayectoria)){
        //se crea un vector llamado partes, ya que cada línea la dividiremos por los espacios en blanco que contiene
        vector<string> partes;
        //dividimos la variable trayectoria por espacios y la almacenamos en el vector partes utilizando la librería boost
        boost::split(partes,trayectoria, boost::is_any_of(" "), boost::token_compress_on);

        
        //Creamos una instancia de la clase átomo
        Atomo atomo;//eliminarlos una vez que los saque del arreglo**
        
        /*
         * A la instancia de átomo le asignamos las posiciones que se obtienen directamente de la trayectoria
         * Utilizamos stod para convertir de string a double
         * Convertimos las posiciones a angstroms al asignarlas
         * */
        atomo.setPos(stod(partes[2])*factor_distancias,stod(partes[3])*factor_distancias,stod(partes[4])*factor_distancias);
        /*
         * A la instancia de átomo le asignamos las posiciones que se obtienen directamente de la trayectoria
         * Utilizamos stod para convertir de string a double
         * Convertimos las velocidades */
        atomo.setSpeeds(stod(partes[5])*factor_velocidades,stod(partes[6])*factor_velocidades,stod(partes[7])*factor_velocidades);

        //La propiedad especie se asignará desde la interfaz gráfica también
        //Sirve para identificar qué tipo de átomo estamos utilizando
        atomo.setEspecie(1);
        //Se crea un arreglo de doubles que almacenará las posiciones periódicas para el átomo
        double periodics[3] = {0.0,0.0,0.0};
        /*
         * Se llama al método que calcula las posiciones periódicas, se le pasan las posiciones en x,y,z, el tamaño de la caja
         * y la mitad de la caja para calcularlas
         * */
        // cout << "estoy en la trayectoria " << trayectorias << " positions_:  x"<< atomo.getrx()  << "y:"<< atomo.getry() << "z:" <<atomo.getrz() << endl;
        calculaPosicionesPeriodicas(periodics, atomo.getrx(), atomo.getry(), atomo.getrz(), boxSize,mitadCaja);
        
        out_posiciones << std::setprecision(20) << "         "<< periodics[0] << "        " << periodics[1] <<"     " << periodics[2]<<endl;
        //Posiciones periódicas ok!!
        
        /*
         * Como se pasó el arreglo periodics por referencia, podemos obtener las posiciones perióicas del mismo arreglo que mandamos
         * */
         
         
        atomo.setPeriodics(periodics[0],periodics[1],periodics[2]);

        /*
         * Una vez que tenemos el átomo con las pocisiones periódicas*/
        (*atomos).push_back(atomo);
        n_atoms++;
        if(n_atoms==numeroAtomos){
           //cout << "ya hay " << numeroAtomos << "hay que llamar a la funcion  \n" << endl;
           histograma(*atomos,numeroAtomos,delta,hist,boxSize,mitadCaja);
           (*atomos).clear();
            n_atoms = 0;
            veces_histograma++;

        }
        trayectorias++;


        //atomo.enviarMensaje();
       // free(line);
    }
    out_posiciones.close();



//    std::ofstream output_file("./histograma.txt");
//    std::ostream_iterator<std::string> output_it(output_file,"\n");
//    std::copy(hist->begin(),hist->end(), output_it);
    //delete atomos;
    
    std::ofstream outfile1;
    if(!prueba){
        outfile1.open(carpetaSalida+"histograma.txt");
    }else{
        outfile1.open(carpetaSalida+"histograma_prueba.txt");
    }
    
    //cout << "Histograma: " << endl;
    for(int contador = 0;contador<tamHistograma;contador++){
        double r = (contador-0.5)*(delta);//es la distancia radial, máximo boxSize/2

       // cout << contador << r << val<< endl;

        //cout << val/veces_histograma << endl;//esto en otro archivo (normalizado)
        contador++;
        outfile1 << contador << "         "<< r << "        " << (*hist)[contador] <<"     "<<endl;
    }
    outfile1.close();

    cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
    cout << "hay " <<trayectorias << " trayectorias" << endl;
    
    std::ofstream outfile2;
    
    

    if(!prueba){
        outfile2.open(carpetaSalida+"histograma_normalizado.txt");    
    }else{
        outfile2.open(carpetaSalida+"histograma_normalizado_prueba.txt");
    }
    
    //cout << "Histograma normalizado: " << endl;
    for(contador = 0;contador<tamHistograma;contador++){
        double r = (contador-0.5)*(delta);//es la distancia radial, máximo boxSize/2



        double norm = (*hist)[contador]/veces_histograma;//esto en otro archivo (normalizado)
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
    if(!prueba){
        outfile3.open(carpetaSalida+"gdr.txt");    
    }else{
        outfile3.open(carpetaSalida+"gdr_prueba.txt");
    }
    
    //cout << "Histograma normalizado: " << endl;
    for(contador=0;contador<tamHistograma;contador++){
        double r = (contador-0.5)*(delta);//es la distancia radial, máximo boxSize/2
         //cout << contador3 << "         "<< r << "        " << norm <<"     "<<endl;
        
        outfile3 << contador <<  std::setprecision(20)<< "         "<< r << "        " << (*gofr)[contador] <<"     "<<endl;

    }
    cout << "En gdr hay: "<<contador<<endl;
    outfile3.close();


    /***********************
    -----------Coordinación promedio----------------
    
    **********************/
    
    
    
    
    
    
    numerosCoordinacion(gofr,coordinaciones,delta,rho);
    
    
    std::ofstream salidaCoordinacion;
    if(!prueba){
        salidaCoordinacion.open(carpetaSalida+"zdr.txt");    
    }else{
        salidaCoordinacion.open(carpetaSalida+"zdr_prueba.txt");
    }
    
    
    for(contador=0;contador<tamHistograma; contador++){
        double r = (contador-0.5)*(delta);//es la distancia radial, máximo boxSize/2
         //cout << contador3 << "         "<< r << "        " << norm <<"     "<<endl;
        salidaCoordinacion << std::setprecision(20)<< r << "        " << (*coordinaciones)[contador] <<"     "<<endl;

    }
    
    
    
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
}

int dist_angulos(int numeroAtomos, std::string carpetaSalida, std::string file, double boxSize){
    cout << "Se calcularán los ángulos" <<endl;
     // Se calcula la mitad de la caja aquí mismo por razones de rendimiento. 
    //double mitadCaja = boxSize/2;
    
    
    // Este factor sirve para convertir de unidades atómicas a Angstroms
    const double factor_distancias =  0.52917720859;
    // Este factor sirve para convertir de unidades átomicas/us a Angstroms/us
    const double factor_velocidades = 21876.912541;
    //
    
     //El archivo del que leeremos las trayectorias
    ifstream input_file;
    
   
    
    vector<double> *gofr = {nullptr};

    //abriremos el archivo de la trayectoria en modo lectura
    input_file.open(file);
    
    int contador;

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
    int veces_histograma = 0;
    
    // Se crea una variable que leerá todo un renglón del archivo
    string trayectoria;
    //aquí, creo un archivo para ver qué posiciones periódicas son las que se obtuvieron (para fines de pruebas)
    std::ofstream out_posiciones;
    //Si no es una prueba, se crea un archivo "definitivo"
    
        out_posiciones.open(carpetaSalida+"periodic_pos.txt");
    
    std::map<Atomo,vector<Atomo>> vecinos;
    
        // de lo contrario, se crea un archivo de prueba (que es más pequeño)
    while (std::getline(input_file,trayectoria)){
        //se crea un vector llamado partes, ya que cada línea la dividiremos por los espacios en blanco que contiene
        vector<string> partes;
        //dividimos la variable trayectoria por espacios y la almacenamos en el vector partes utilizando la librería boost
        boost::split(partes,trayectoria, boost::is_any_of(" "), boost::token_compress_on);
        double mitadCaja = boxSize/2;
        
        //Creamos una instancia de la clase átomo
        Atomo atomo;//eliminarlos una vez que los saque del arreglo**
        
        /*
         * A la instancia de átomo le asignamos las posiciones que se obtienen directamente de la trayectoria
         * Utilizamos stod para convertir de string a double
         * Convertimos las posiciones a angstroms al asignarlas
         * */
        atomo.setPos(stod(partes[2])*factor_distancias,stod(partes[3])*factor_distancias,stod(partes[4])*factor_distancias);
        /*
         * A la instancia de átomo le asignamos las posiciones que se obtienen directamente de la trayectoria
         * Utilizamos stod para convertir de string a double
         * Convertimos las velocidades */
        atomo.setSpeeds(stod(partes[5])*factor_velocidades,stod(partes[6])*factor_velocidades,stod(partes[7])*factor_velocidades);

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
        
        
        
        /*
         * Como se pasó el arreglo periodics por referencia, podemos obtener las posiciones perióicas del mismo arreglo que mandamos
         * */
         
         
        atomo.setPeriodics(periodics[0],periodics[1],periodics[2]);

        /*
         * Una vez que tenemos el átomo con las pocisiones periódicas*/
        (*atomos).push_back(atomo);
        n_atoms++;
        if(n_atoms==numeroAtomos){
           //cout << "ya hay " << numeroAtomos << "hay que llamar a la funcion  \n" << endl;
           //aquí debo llamar a la función que calcule la lista de vecinos
           listaVecinos(*atomos,n_atoms,3.2,mitadCaja, boxSize, vecinos);
           (*atomos).clear();
            n_atoms = 0;

        }
        trayectorias++;


        atomo.enviarMensaje();
       // free(line);
    }

    return 0;
}