//
// Created by Erick on 15/06/2020.
//
#include "utils.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include </usr/include/fftw3.h>
#include <map>

struct Celda
{
    std::vector<Atomo *> atomos;
};

void calculaPosicionesPeriodicas(double *posPeriodicas, const double rx, const double ry, const double rz, double boxSize, double halfBox)
{
    posPeriodicas[0] = evaluaCajaRecursivo(rx, boxSize, halfBox);
    posPeriodicas[1] = evaluaCajaRecursivo(ry, boxSize, halfBox);
    posPeriodicas[2] = evaluaCajaRecursivo(rz, boxSize, halfBox);
}

double evaluaCajaRecursivo(double posicion, double boxSize, double halfBox)
{
    double residuo = fmod(posicion, boxSize);
    if (residuo > halfBox)
    {
        residuo -= boxSize;
    }
    else if (residuo < -halfBox)
    {
        residuo += boxSize;
    }

    return residuo;
}

/*double evaluaCajaRecursivo(double posicion, double boxSize, double halfBox){
    bool reaplicar = false;
    if(posicion >= halfBox){
        posicion -= boxSize;
        posicion = evaluaCajaRecursivo(posicion,boxSize,halfBox);
    } else if(posicion < -halfBox){
        posicion += boxSize;
        posicion = evaluaCajaRecursivo(posicion,boxSize,halfBox);
    }else{
        return posicion;
    }
    return posicion;
}*/

void normalizarHistograma(vector<double> histo, vector<double> &histo_norm, int n_atomos, int tam_histograma, double delta, double boxSize, int vecesHistograma)
{
    double volumen = boxSize * boxSize * boxSize;
    // double pi = M_PI;
    double factor_normalizacion = volumen / (2 * M_PI * n_atomos * n_atomos * delta * vecesHistograma);
    std::cout << factor_normalizacion << std::endl;
    double r = 0;
    int i = 0;

    for (int el_hist : histo)
    {
        r = (i - 0.5) * delta;
        histo_norm[i] = (el_hist * factor_normalizacion) / (r * r);
        // std::cout << histo_norm[i] << std::endl;
        i++;
    }
}

void histograma(const vector<Atomo> &atomos, int n_atomos, double delta, vector<double> &histo, double boxSize, double mitadCaja)
{

    vector<Atomo> vecinos;

    int no_cumplen = 0;
    double distancias[3] = {0, 0, 0};
    for (int iteracion = 0; iteracion < n_atomos; iteracion++)
    {

        for (int it2 = iteracion + 1; it2 < n_atomos; it2++)
        {
            // std::cout << "Lectura   "<< it2 << " y " << iteracion << std::endl;

            calcula_dist_componentes(distancias, atomos.at(iteracion), atomos.at(it2), boxSize, mitadCaja);

            // std::cout << "Lectura   "<< j << std::endl;

            double distancia2 = calculaDistancia(distancias);
            //
            // double distancia = sqrtf(distancia2);
            double distancia = sqrt(distancia2);
            int position = 0;

            if (distancia < mitadCaja)
            { // falta obtener las imágenes también de aquí, restando a la distancia el tamaño de la caja
                position = distancia / delta;
                if (position < histo.size())
                {
                    histo[position] += 1;
                }
                else
                {
                    // std::cout << "La posición " << position <<"Se sale del histograma" << std::endl;
                }
            }
        }
    }
}

/*void histogramaMejorado(vector<Atomo> atomos,int n_atomos, double delta, vector<double> *histo, double boxSize, double mitadCaja){
    // std::cout << "Calculando histograma " << std::endl;
    //vector<int> hist(tamanio,0);
    //TimerBaseChrono tc("comienza histograma");
    //tc.Start();
    int no_cumple = 0;

    for (int iteracion = 0; iteracion<n_atomos;iteracion++) {
        //TimerBaseChrono tdist("recorrer atomo ");
        //tdist.Start();
        for ( int it2=iteracion+1;it2<n_atomos; it2++) {
            //std::cout << "Lectura   "<< i << std::endl;
            double distancias[3] = {0,0,0};

            Atomo atomo1 = atomos[iteracion];
            Atomo atomo2 = atomos[it2];
            calcula_dist_componentes(distancias,atomo1,atomo2,boxSize,mitadCaja);

            // std::cout << "Lectura   "<< j << std::endl;

            double distancia2 = calculaDistancia(distancias);
            //
            double distancia = sqrtf(distancia2);
            int position = 0;

            if(distancia<mitadCaja){
                position = distancia/delta;
                (*histo)[position] += 1;

            }else{
                atomo1.setPeriodics(atomo1.getPrx()+boxSize, atomo1.getPry()+boxSize, atomo1.getPrz()+boxSize);
                calcula_dist_componentes(distancias, atomo1, atomo2, boxSize, mitadCaja);
                distancia2 = calculaDistancia(distancias);
                distancia = sqrtf(distancia2);
                distancia = distancia-mitadCaja;
                if(distancia<mitadCaja){
                position = distancia/delta;
                (*histo)[position] += 1;

            }else{
                //std::cout << "la distancia " << distancia << "no cumple el criterio" << std::endl;
                no_cumple++;
            }

                // std::cout << "la distancia " << distancia << "no cumple el criterio" << std::endl;
            }

            //std::cout << "recorrer 300 tarda "<< tdist.GetMs() << std::endl;
        }
        //tc.getTime();
    }


}*/

void calcula_dist_componentes(double *distancias, Atomo a1, Atomo a2, double boxSize, double halfBox)
{

    double rx1 = a1.getrx();
    double rx2 = a2.getrx();
    double ry1 = a1.getry();
    double ry2 = a2.getry();
    double rz1 = a1.getrz();
    double rz2 = a2.getrz();

    double compx = rx1 - rx2;
    // std::cout << perx2 <<"-"<<perx1<<"=" <<compx << std::endl;

    double compy = ry1 - ry2;
    // std::cout << pery2 <<"-"<<pery1 <<"=" << compy << std::endl;

    double compz = rz1 - rz2;
    // std::cout <<" aquí sí llega"<< std::endl;

    distancias[0] = verificaComponenteParaImagen(compx, boxSize, halfBox);
    distancias[1] = verificaComponenteParaImagen(compy, boxSize, halfBox);
    distancias[2] = verificaComponenteParaImagen(compz, boxSize, halfBox);
}
/*void calcula_dist_componentesAnt(double *distancias, Atomo a1, Atomo a2, double boxSize, double halfBox){

    double perx1 = a1.getPrx();
    double perx2 = a2.getPrx();
    double pery1 = a1.getPry();
    double pery2 = a2.getPry();
    double perz1 = a1.getPrz();
    double perz2 = a2.getPrz();
    double compx = perx2-perx1;
    //std::cout << perx2 <<"-"<<perx1<<"=" <<compx << std::endl;


    double compy = pery2-pery1;
    //std::cout << pery2 <<"-"<<pery1 <<"=" << compy << std::endl;

    double compz = perz2-perz1;
    //std::cout <<" aquí sí llega"<< std::endl;

    distancias[0]=verificaComponenteParaImagen(compx,boxSize, halfBox);
    distancias[1]=verificaComponenteParaImagen(compy,boxSize, halfBox);
    distancias[2]=verificaComponenteParaImagen(compz,boxSize, halfBox);

}*/

double verificaComponenteParaImagen(double dif_componentes, double boxSize, double halfBox)
{ // usar esta angulos

    double residuo = fmod(dif_componentes, boxSize);
    // std::cout <<residuo << dif_componentes << std::endl;
    if (residuo > halfBox)
    {
        residuo -= boxSize;
    }
    else if (residuo < -halfBox)
    {
        residuo += boxSize;
    }

    return residuo;
}

/*double verificaComponenteParaImagen(double dif_componentes, double boxSize, double halfBox){//usar esta angulos
    int veces;
    veces = floor(dif_componentes/(boxSize));
    std::cout << "modulo: " << dif_componentes << "   " << boxSize <<"       "<< veces << std::endl;



    if(dif_componentes > 0){
        return dif_componentes>(halfBox)?dif_componentes-(boxSize*veces):dif_componentes;
    }
    if(dif_componentes < 0){
        return dif_componentes<(-halfBox)?dif_componentes+(boxSize*veces):dif_componentes;//esto es correcto, lo verifiqué y funciona
    }

    return dif_componentes;
}*/

double verificaComponenteR(double componente, double boxSize, double halfBox)
{
    int veces;
    veces = floor(componente / halfBox);
    if (componente > halfBox)
    {
        componente = componente - (boxSize * veces);
        componente = verificaComponenteR(componente, boxSize, halfBox);
    }
    else if (componente < -halfBox)
    {
        componente = componente + (boxSize * veces);
        componente = verificaComponenteR(componente, boxSize, halfBox);
    }

    return componente;
}

double calculaDelta(double boxSize, int tamHistograma)
{
    return (boxSize / 2) / tamHistograma;
    // return 10.0/tamHistograma;
}

double calculaDistancia(double *distancias)
{

    double sum2 = 0;
    for (int i = 0; i < 3; i++)
    {
        sum2 += distancias[i] * distancias[i];
    }

    return sum2;
}

double fr(double gr, double r)
{
    return gr * r * r;
}
double integral_simple(double gr, double r, double delta)
{
    double lim_inf = r;
    double lim_sup = lim_inf + delta;
    double h = lim_sup - lim_inf;
    return h * ((fr(gr, lim_inf) + fr(gr, lim_sup)) / 2);
}
double integral(vector<double> funcion_integrar, int n_divisiones)
{ // determinar quin será a y hasta dónde llegará b
    double tam_funcion = funcion_integrar.size();
    double a = 0.0;                    // el primer extremo en la función (valor en las x )
    double b = 180.0;                  // el valor en las x
    double h = (b - a) / n_divisiones; //(b-a)/n;

    double fa = funcion_integrar.front(); // son las amplitudes inicial
    double fb = funcion_integrar.back();  // amplitud final
    double sum = 0;

    for (double i = 1; i < tam_funcion; i++)
    {

        sum += funcion_integrar[i];
    }

    return h * (((fa + fb) / 2) + sum); //*/
    // return 0.0;
}

void numerosCoordinacion(vector<double> gdr, vector<double> &coordinaciones, double delta, double rho)
{
    double tam = gdr.size();

    double suma = 0;
    // vector de coordinaciones

    for (int i = 0; i < tam; i++)
    {
        double r = (i - 0.5) * delta;
        double n_c = integral_simple(gdr[i], r, delta); // densidad, n_partixulas/volumen
        double n_coor = n_c * 4 * M_PI * rho;

        // double tota = actual+anterior;
        suma += n_coor;
        coordinaciones[i] = suma;
        // anterior = n_coor;
        // imprimo a un archivo
    }
}

void factorEstructuraE(vector<double> gdr, int tamHistograma, double delta_k, double rho, string dirSalida, double delta)
{
    // std::cout << "Entramos a la transformada" << std::endl;

    // double factor = 16*rho*((2*M_PI)/delta_k*tamHistograma);
    int MAX = tamHistograma * 16; // el 16 será variable
    int N = MAX / 2;

    fftw_complex out[MAX], in[MAX];
    int i, j;

    for (i = 0; i < MAX; i++)
    {
        if (i < tamHistograma)
        {
            in[i][0] = gdr[i];
        }
        else
        {
            in[i][0] = 1;
        }
        // std::cout << in[i][0] << std::endl;
    }

    std::ofstream testFile;

    testFile.open(dirSalida + "/" + "test_file.txt");
    for (i = 0; i < MAX; i++)
    {
        double r = delta * i;
        double valAnterior = in[i][0];
        // std::cout << in[i][0] << std::endl;
        in[i][0] = (in[i][0] - 1) * r;
        // std::cout << in[i][0] << std::endl;

        testFile << std::setprecision(15) << "     " << r << "       " << valAnterior << "     " << in[i][0] << "     " << std::endl;
        // std::cout << in[i][0] << std::endl;
    }
    testFile.close();

    fftw_plan p;

    p = fftw_plan_dft_1d(MAX, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    // j=-(N);
    std::ofstream vOndaPos;
    string outVOndaPos = dirSalida + "/" + "vectores_onda_pos.txt";
    vOndaPos.open(outVOndaPos);

    vOndaPos << " #          k                            F(Re)                        F(Im)                      norm(Re)                                  norm(Im)            " << std::endl;
    for (i = 0; i < N; i++)
    {

        double k = i * delta_k;

        double real = ((out[i][0] * 32 * M_PI * rho) / (k * tamHistograma)) + 1.0;
        double im = -((out[i][1] * 32 * M_PI * rho) / (k * tamHistograma)) + 1.0;
        vOndaPos << std::setprecision(15) << i << "         " << k << "                " << out[i][0] << "               " << out[i][1] << "                    " << real << "      " << im << std::endl;
        //  j++;
        // en función de i, y en función de la frecuencia
    }
    std::cout << "vectores de onda pos:" << outVOndaPos << std::endl;

    j = N;
    std::ofstream vOndaNeg;
    string outVOndaNeg = dirSalida + "/" + "vectores_onda_neg.txt";
    vOndaNeg.open(outVOndaNeg);
    for (i = N; i < MAX; i++)
    {
        double k = -j * delta_k;
        // double w = 2*M_PI*f;
        j--;
        double val = out[i][0];
        vOndaNeg << std::setprecision(20) << i << "     " << k << "         " << val << "     " << out[i][1] << std::endl;
    }
    std::cout << "vectores onda neg:" << outVOndaNeg << std::endl;

    fftw_destroy_plan(p);
    fftw_cleanup();
}

void sk(double *gdr, int tamHistograma, double delta_k, double rho, string dirSalida, double delta)
{
    /*
    int MAX = 8192;
    int N=MAX/2;
    fftw_complex out[MAX], in[MAX];
    int i;

    for(i=0;i<MAX;i++){
        if(i<tamHistograma){
            in[i][0] = gdr[i];
            //std::cout << (*gdr)[i] << std::endl;
        }else{
        in[i][0]=1;

        }
        //std::cout << in[i][0] << std::endl;
    }

    for(i=0; i<MAX; i++){
        double r = i*delta;
        in[i][0] = (in[i][0]-1)*r;
    }


    fftw_plan p;

    p = fftw_plan_dft_1d(MAX,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(p);

    //j=-(N);
    std::ofstream vOndaPos;
    vOndaPos.open(dirSalida+"vectores_onda_pos.txt");

    for(i=0; i<N;i++){

        double k = i*delta_k;


        double real = ((out[i][0]*16*rho*2*M_PI)/(k*tamHistograma))+1.0;
        double im = -((out[i][1]*16*rho*2*M_PI)/(k*tamHistograma))+1.0;
        vOndaPos << std::setprecision(15) <<"     "<< k <<"       "<< real <<"     " << im << std::endl;

    }





    fftw_destroy_plan(p);
    fftw_cleanup();

    */
}

float distanciaAtomos(Atomo a1, Atomo a2)
{
    return a1.distancia(a2.getPrx(), a2.getPry(), a2.getPrz());
}

void vecinosMejorada(vector<Atomo> atomos, double r_min, double boxSize, double mitadCaja, vector<double> *histAngulos)
{
    int numCeldasPorDim = ceil(boxSize / r_min);
    double deltaAng = 180.0 / 500;
    std::vector<std::vector<std::vector<Celda>>> celdas(numCeldasPorDim, std::vector<std::vector<Celda>>(numCeldasPorDim, std::vector<Celda>(numCeldasPorDim)));

    for (Atomo &atomo : atomos)
    {
        int x = floor(atomo.getPrx() / r_min);
        int y = floor(atomo.getPry() / r_min);
        int z = floor(atomo.getPrz() / r_min);
        std::cout << "Antes: x=" << x << ", y=" << y << ", z=" << z << std::endl;
        x = (x + numCeldasPorDim) % numCeldasPorDim;
        y = (y + numCeldasPorDim) % numCeldasPorDim;
        z = (z + numCeldasPorDim) % numCeldasPorDim;
        std::cout << "Despues: x=" << x << ", y=" << y << ", z=" << z << std::endl;
        int otro = 0;
        celdas[x][y][z].atomos.push_back(&atomo);
    }

    // 183
    std::vector<std::vector<int>> offsets = {
        {-1, -1, -1},
        {0, -1, -1},
        {1, -1, -1},
        {-1, 0, -1},
        {0, 0, -1},
        {1, 0, -1},
        {-1, 1, -1},
        {0, 1, -1},
        {1, 1, -1},
        {-1, -1, 0},
        {0, -1, 0},
        {1, -1, 0},
        {-1, 0, 0},
        {0, 0, 0},
        {1, 0, 0},
        {-1, 1, 0},
        {0, 1, 0},
        {1, 1, 0},
        {-1, -1, 1},
        {0, -1, 1},
        {1, -1, 1},
        {-1, 0, 1},
        {0, 0, 1},
        {1, 0, 1},
        {-1, 1, 1},
        {0, 1, 1},
        {1, 1, 1},
    };
    std::map<Atomo *, std::vector<Atomo *>> listaVecinos;

    for (int i = 0; i < numCeldasPorDim; ++i)
    {
        for (int j = 0; j < numCeldasPorDim; ++j)
        {
            for (int k = 0; k < numCeldasPorDim; ++k)
            {
                Celda &celda = celdas[i][j][k];

                for (Atomo *atomo1 : celda.atomos)
                {
                    for (auto &offset : offsets)
                    {
                        int ii = (i + offset[0] + numCeldasPorDim) % numCeldasPorDim;
                        int jj = (j + offset[1] + numCeldasPorDim) % numCeldasPorDim;
                        int kk = (k + offset[2] + numCeldasPorDim) % numCeldasPorDim;

                        if (ii >= 0 && ii < numCeldasPorDim &&
                            jj >= 0 && jj < numCeldasPorDim &&
                            kk >= 0 && kk < numCeldasPorDim)
                        {
                            Celda &celda_vecina = celdas[ii][jj][kk];

                            for (Atomo *atomo2 : celda_vecina.atomos)
                            {
                                if (atomo1->getId() != atomo2->getId())
                                {
                                    if (verificaVecindad(*atomo1, *atomo2, r_min, mitadCaja, boxSize))
                                    {
                                        listaVecinos[atomo1].push_back(atomo2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    std::map<Atomo, vector<Atomo>>::iterator it;

    for (auto &pair : listaVecinos)
    {
        Atomo *atomoA = pair.first;
        auto &vecinos = pair.second;
        // Atomo atomoA = it->first;
        // obtener los vecinos de a1
        int id = atomoA->getId();
        // vector<Atomo> vecinosTriada = ordenarPorDistancia(it->second);
        // vector<Atomo> vecinos = it->second;
        // std::cout << "Atomo " << it->first.getId() << " tiene " << it->second.size() << " vecinos" << std::endl;
        for (int i = 0; i < vecinos.size() - 1; i++)
        {
            // iterar después del atomo actual
            Atomo *atomoB = vecinos[i];
            Atomo *atomoC = vecinos[i + 1];

            std::cout << "Base = " << id << "vecinos" << atomoB->getId() << "    " << atomoC->getId() << std::endl;

            vector<double> rAB = atomoA->vectorDifference(*atomoB);
            vector<double> rCB = atomoC->vectorDifference(*atomoB);

            // calcular el producto punto de los vectores v1 y v2
            double productoPunto = rAB[0] * rCB[0] + rAB[1] * rCB[1] + rAB[2] * rCB[2];
            // calcular el modulo del vector v1
            double mrAB = sqrt(rAB[0] * rAB[0] + rAB[1] * rAB[1] + rAB[2] * rAB[2]);
            // calcular el modulo del vector v2
            double mrCB = sqrt(rCB[0] * rCB[0] + rCB[1] * rCB[1] + rCB[2] * rCB[2]);
            // calcular el angulo entre los vectores v1 y v2
            double anguloRad = acos(productoPunto / (mrAB * mrCB));
            // convertir el angulo a grados
            double anguloGrad = anguloRad * (180 / M_PI);

            int pos = 0;

            pos = anguloGrad / deltaAng;
            // imprimir la posicion en la que quedará
            // std::cout << "Pos = "<< anguloGrad << "/" << deltaAng <<"="<<pos << std::endl;
            // incrementar esa posicion en el histograma
            (*histAngulos)[pos] += 1;
            // imprimir cuantos elementos en esa posicion hay
            std::cout << "Histograma[" << pos << "]=" << pos << std::endl;
        }
        vecinos.clear();
    }
}

void listaVecinos(vector<Atomo> atomos, int n_atomos, float r_min, double mitadCaja, double boxSize, std::map<Atomo, vector<Atomo>> vecinos, vector<double> *histAngulos, double deltaAng, std::string trayectoria)
{
    // file:///D:/tesis/libros/computer_simultation_of_liquids.pdf pp. 162
    // vector con vectores (triadas) que serán iterados después
    // vector<vector<Atomo>> vec;
    bool esVecino = false;
    int cont1 = 0;
    int cont2 = 0;
    for (auto &aa1 : atomos)
    {
        cont1++;
        for (auto aa2 : atomos)
        {
            cont2++;
            Atomo atomo1 = aa1;
            Atomo atomo2 = aa2;

            if (atomo1.getId() != atomo2.getId())
            {
                bool esVecino = verificaVecindad(atomo1, atomo2, r_min, mitadCaja, boxSize);
                if (esVecino)
                {
                    vecinos[atomo1].push_back(atomo2);
                }
            }
        }
    }
    std::map<Atomo, vector<Atomo>>::iterator it;
    std::ofstream problematicos;
    std::ofstream probs;
    // problematicos.open("/home/erick/data/salidas/input_problematicas.txt");
    // probs.open("/home/erick/data/salidas/reduced_problematics.txt");
    for (it = vecinos.begin(); it != vecinos.end(); it++)
    {
        Atomo atomoA = it->first;
        // obtener los vecinos de a1
        int id = atomoA.getId();
        // vector<Atomo> vecinosTriada = ordenarPorDistancia(it->second);
        vector<Atomo> vecinos = it->second;
        // std::cout << "Atomo " << it->first.getId() << " tiene " << it->second.size() << " vecinos" << std::endl;
        for (int i = 0; i < vecinos.size() - 1; i++)
        {
            // iterar después del atomo actual
            Atomo atomoB = vecinos[i];
            Atomo atomoC = vecinos[i + 1];

            // std::cout << "Base = "<< id << "vecinos" << atomoB.getId() <<"    "<< atomoC.getId() << std::endl;
            // comparar distancias y obtener el central aquí
            double anguloGrad = anguloConAtomoMasCercano(atomoA, atomoB, atomoC);
            problematicos << atomoA.getId() << "  " << atomoB.getId() << " " << atomoC.getId() << "   " << anguloGrad << std::endl;
            if (anguloGrad < 12 && anguloGrad > 10)
            {

                std::cout << "Problema aqui" << std::endl;
                std::cout << atomoA.enviarMensaje() << std::endl;
                std::cout << atomoB.enviarMensaje() << std::endl;
                std::cout << atomoC.enviarMensaje() << std::endl;
                // std::cout << "periodicas" << atomoA.getPrx() <<" " << atomoB.getPrx() << std::endl;
                // probs << trayectoria << std::endl;
            }
            int pos = 0;

            pos = anguloGrad / deltaAng;
            // imprimir la posicion en la que quedará
            // std::cout << "Pos = "<< anguloGrad << "/" << deltaAng <<"="<<pos << std::endl;
            // incrementar esa posicion en el histograma
            (*histAngulos)[pos] += 1;
            // imprimir cuantos elementos en esa posicion hay
            // std::cout << "Histograma["<<pos<<"]="<<(*histAngulos)[pos]<<std::endl;
        }
        vecinos.clear();
    }
    problematicos.close();
    probs.close();
}

bool verificaVecindad(Atomo a1, Atomo a2, double r_min, double mitadCaja, double boxSize)
{
    double distancias_comp[3] = {0, 0, 0};
    calcula_dist_componentes(distancias_comp, a1, a2, boxSize, mitadCaja);

    double distancia2 = calculaDistancia(distancias_comp);
    double distancia = sqrt(distancia2);

    bool esVecino = false;
    // std::cout << "Distancia="<<distancia<<std::endl;
    if (distancia <= mitadCaja)
    {
        if (distancia < r_min)
        {
            esVecino = true;
        }
        else
        {
            esVecino = false;
        }
    }

    return esVecino;
}

double calculaAngulo(Atomo &atomoCentral, Atomo &atomoB, Atomo &atomoC)
{

    vector<double> rAB = atomoCentral.vectorDifference(atomoB);
    vector<double> rAC = atomoCentral.vectorDifference(atomoC);

    // calcular el producto punto de los vectores v1 y v2
    double productoPunto = rAB[0] * rAC[0] + rAB[1] * rAC[1] + rAB[2] * rAC[2];
    // calcular el modulo del vector v1
    double mrAB = sqrt(rAB[0] * rAB[0] + rAB[1] * rAB[1] + rAB[2] * rAB[2]);
    // calcular el modulo del vector v2
    double mrAC = sqrt(rAC[0] * rAC[0] + rAC[1] * rAC[1] + rAC[2] * rAC[2]);
    // calcular el angulo entre los vectores v1 y v2
    double anguloRad = acos(productoPunto / (mrAB * mrAC));
    // convertir el angulo a grados
    double anguloGrad = anguloRad * 180 / M_PI;

    return anguloGrad;
}

void obtenerArgumentos(int argc, char **argv, double &boxSize, int &numeroAtomos, int &tamHistograma, std::string &carpetaSalida, std::string &file, int &opc)
{
    if (argc < 2)
    {
        // Es el tamaño de la caja de simulación que se utilizará
        boxSize = 22.3797;
        // Se pide el número de átomos que se utilizaron en la dimámica molecular
        numeroAtomos = 300;
        // Se pide el tamaño del histograma que se utilizará
        tamHistograma = 512;
        carpetaSalida = "/home/erick/data/salidas";
        file = "/home/erick/data/Ge00Sb00Te100T823K.cpmd"; //"/home/erick/protocolo/copia.cpmd";////"/home/erick/protocolo/TRAJECTORY_00_Te_T823K.cpmd"; Ge00Sb00Te100T823K.cpmd
        // file = "out.xx";
        opc = 1;
    }
    else
    {
        boxSize = std::stod(argv[1]);
        numeroAtomos = std::stoi(argv[2]);
        tamHistograma = std::stoi(argv[3]);
        file = argv[4];
        carpetaSalida = argv[5];
        opc = std::stoi(argv[6]);
    }
}

double distanciaEntreAtomos(Atomo &a, Atomo &b)
{

    return a.distancia(b.getPrx(), b.getPry(), b.getPrz());
}

double anguloConAtomoMasCercano(Atomo &a, Atomo &b, Atomo &c)
{
    double distanciaAB = distanciaEntreAtomos(a, b);
    double distanciaAC = distanciaEntreAtomos(a, c);
    double distanciaBC = distanciaEntreAtomos(b, c);

    double distanciaTotalA = distanciaAB + distanciaAC;
    double distanciaTotalB = distanciaAB + distanciaBC;
    double distanciaTotalC = distanciaAC + distanciaBC;

    /*std::cout << "Distancia A" << distanciaTotalA<<std::endl;
    std::cout << "Distancia B" << distanciaTotalB<<std::endl;
    std::cout << "Distancia C" << distanciaTotalC<<std::endl;*/
    // Devuelve el átomo con la menor distancia total.
    if (distanciaTotalA <= distanciaTotalB && distanciaTotalA <= distanciaTotalC)
    {
        // std::cout << "Central A"<<std::endl;
        return calculaAngulo(a, b, c);
    }
    else if (distanciaTotalB <= distanciaTotalA && distanciaTotalB <= distanciaTotalC)
    {
        // std::cout << "Central B"<<std::endl;
        return calculaAngulo(b, a, c);
    }
    else
    {
        //  std::cout << "Central C"<<std::endl;
        return calculaAngulo(c, a, b);
    }
}

void obtenerArgumentosAng(int argc, char **argv, double &boxSize, int &numeroAtomos, int &tamHistograma, std::string &carpetaSalida, std::string &file, int &opc, float r_min)
{
    if (argc < 2)
    {
        // Es el tamaño de la caja de simulación que se utilizará
        boxSize = 22.3797;
        // Se pide el número de átomos que se utilizaron en la dimámica molecular
        numeroAtomos = 300;
        // Se pide el tamaño del histograma que se utilizará
        tamHistograma = 512;
        carpetaSalida = "/home/erick/data/salidas";
        file = "/home/erick/data/Ge00Sb00Te100T823K.cpmd"; //"/home/erick/protocolo/copia.cpmd";////"/home/erick/protocolo/TRAJECTORY_00_Te_T823K.cpmd"; Ge00Sb00Te100T823K.cpmd
        // file = "out.xx";
        opc = 1;
        r_min = 3.0;
    }
    else
    {
        boxSize = std::stod(argv[1]);
        numeroAtomos = std::stoi(argv[2]);
        tamHistograma = std::stoi(argv[3]);
        file = argv[4];
        carpetaSalida = argv[5];
        opc = std::stoi(argv[6]);
        r_min = std::stod(argv[7]);
    }
}
