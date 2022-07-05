//
// Created by Erick on 15/06/2020.
//
#include "utils.h"
#include "cmath"
#include <iostream>
#include <fstream>
#include <iomanip>
#include </usr/include/fftw3.h>

void calculaPosicionesPeriodicas(double *posPeriodicas,const double rx, const double ry, const double rz, double boxSize, double halfBox){


        posPeriodicas[0] = evaluaCajaRecursivo(rx,boxSize, halfBox);
        posPeriodicas[1] = evaluaCajaRecursivo(ry,boxSize, halfBox);
        posPeriodicas[2] = evaluaCajaRecursivo(rz,boxSize, halfBox);


}


double evaluaCajaRecursivo(double posicion, double boxSize, double halfBox){
    
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
}

void normalizarHistograma(vector<double> *histo, vector<double> *histo_norm, int n_atomos, int tam_histograma, double delta, double boxSize, int vecesHistograma){
    double volumen = boxSize*boxSize*boxSize;
    //double pi = M_PI;
    double factor_normalizacion = volumen/(2*M_PI*n_atomos*n_atomos*delta*vecesHistograma);
    std::cout << factor_normalizacion << std::endl;
    double r = 0;
    int i = 0;


    for(int el_hist:*histo){
        r=(i-0.5)*delta;
        (*histo_norm)[i]=(el_hist*factor_normalizacion)/(r*r);
        i++;
    }
}

void histograma(vector<Atomo> atomos,int n_atomos, double delta, vector<double> *histo, double boxSize, double mitadCaja){
    // std::cout << "Calculando histograma " << std::endl;
    //vector<int> hist(tamanio,0);
    //TimerBaseChrono tc("comienza histograma");
    //tc.Start();
    vector<Atomo> vecinos;
    
    int no_cumplen = 0;
    double distancias[3] = {0,0,0};
    for (int iteracion = 0; iteracion<n_atomos;iteracion++) {
        //TimerBaseChrono tdist("recorrer atomo ");
        //tdist.Start();
        for ( int it2=iteracion+1;it2<n_atomos; it2++) {
            //std::cout << "Lectura   "<< i << std::endl;
            
                      
            calcula_dist_componentes(distancias,atomos[iteracion],atomos[it2],boxSize,mitadCaja);

            // std::cout << "Lectura   "<< j << std::endl;

            double distancia2 = calculaDistancia(distancias);
            //
            //double distancia = sqrtf(distancia2);
            double distancia = sqrt(distancia2);
            int position = 0;

            if(distancia<10.0){
                position = distancia/delta;
                (*histo)[position] += 1;

            }/*else{
                    no_cumplen++;
                    //std::cout << "la distancia " << distancia << "no cumple el criterio" << std::endl;
                
                // std::cout << "la distancia " << distancia << "no cumple el criterio" << std::endl;
            }*/

            //std::cout << "recorrer 300 tarda "<< tdist.GetMs() << std::endl;
        }
        //tc.getTime();
    }
    //std::cout << "por pasada no cumplen: "<< no_cumplen << std::endl;

}

void histogramaMejorado(vector<Atomo> atomos,int n_atomos, double delta, vector<double> *histo, double boxSize, double mitadCaja){
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
    

}

void calcula_dist_componentes(double *distancias, Atomo a1, Atomo a2, double boxSize, double halfBox){
    
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
    //std::cout << perz2 <<"-"<<perz1 <<"="<<compz << std::endl;
    
    distancias[0]=verificaComponenteR(compx,boxSize, halfBox);
    distancias[1]=verificaComponenteR(compy,boxSize, halfBox);
    distancias[2]=verificaComponenteR(compz,boxSize, halfBox);
    
}



double verificaDistancia(double distancia, double boxSize, double halfBox){
    return distancia>=(halfBox)?distancia-boxSize:distancia;
/*    if(distancia>halfBox){
        return distancia-boxSize;
    }else{
        return distancia;
    }*/
}
double verificaComponente(double distancia, double boxSize, double halfBox){
    
    if(distancia > 0){
        return distancia>(halfBox)?distancia-boxSize:distancia;
    }
    if(distancia < 0){
        return distancia<(-halfBox)?distancia+boxSize:distancia;
    }
    
    return distancia;
}

double verificaComponenteR(double componente, double boxSize, double halfBox){
    if(componente >= halfBox){
        componente = componente - boxSize;
        componente = verificaComponenteR(componente,boxSize,halfBox);
    } else if(componente <= -halfBox){
        componente = componente + boxSize;
        componente = evaluaCajaRecursivo(componente,boxSize,halfBox);
        
    }
    
    return componente;
}



double calculaDelta(double boxSize, int tamHistograma){
    //return (boxSize/2)/tamHistograma;
    return 10.0/tamHistograma;
}

double calculaDistancia(double *distancias){


    double sum2 = 0;
    for(int i=0;i<3;i++){
        sum2 += distancias[i]*distancias[i];
    }

    return sum2;

}

double fr(double gr, double r){
    return gr*r*r;
}
double integral_simple(double gr, double r, double delta){
    double lim_inf = r;
    double lim_sup = lim_inf+delta;
    double h = lim_sup-lim_inf;
    return h*((fr(gr,lim_inf)+fr(gr,lim_sup))/2);    
    
}
double integral(double gr, double r, int n, double delta){ //determinar quin será a y hasta dónde llegará b
    
    double h = delta;//(b-a)/n;
    double a = r;
    double b = a+delta;
    
    double fra = fr(gr,a);
    
    double sum = 0;
    
    for(double i = a+h; i<b;i+=h){
        
        sum+=fr(gr,i);
    }
    
    
    double frb = fr(gr,b);
    
    
    
    return (b-a)*((fra+2*sum+frb)/(2*n));//*/
    //return 0.0;

}


void numerosCoordinacion(vector<double> *gdr, vector<double> *coordinaciones, double delta, double rho) {
    double tam = gdr->size();


    double suma = 0;
    //vector de coordinaciones


    for(int i = 0; i<tam; i++){
        double r = (i-0.5)*delta;
        double n_c = integral_simple((*gdr)[i],r, delta);//densidad, n_partixulas/volumen
        double n_coor = n_c*4*M_PI*rho;

        //double tota = actual+anterior;
        suma += n_coor;
        (*coordinaciones)[i]=suma;
        //anterior = n_coor;
        //imprimo a un archivo

    }

}

void factorEstructuraE(vector<double> *gdr, int tamHistograma, double delta_k, double rho, string dirSalida, double delta){
    //std::cout << "Entramos a la transformada" << std::endl;
    
    //double factor = 16*rho*((2*M_PI)/delta_k*tamHistograma);
    int MAX = tamHistograma*16;//el 16 será variable
    int N=MAX/2;
    
    fftw_complex out[MAX], in[MAX];
    int i,j;
    
    for(i=0;i<MAX;i++){
        if(i<tamHistograma){
            in[i][0] = (*gdr)[i];
            
        }else{
        in[i][0]=1;
        
        }
        //std::cout << in[i][0] << std::endl;
    }
    
    std::ofstream testFile;
    
    testFile.open(dirSalida+"test_file.txt");
    for(i=0;i<MAX;i++){
        double r = delta*i;
        double valAnterior = in[i][0];
        //std::cout << in[i][0] << std::endl;
        in[i][0] = (in[i][0]-1)*r;
        //std::cout << in[i][0] << std::endl;
        
         testFile << std::setprecision(15) <<"     "<< r <<"       "<< valAnterior <<"     " << in[i][0] <<"     "  << std::endl; 
        //std::cout << in[i][0] << std::endl;
    }
    testFile.close();
    
    
    fftw_plan p;
    
    p = fftw_plan_dft_1d(MAX,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(p);
    
    //j=-(N);
    std::ofstream vOndaPos;
    vOndaPos.open(dirSalida+"vectores_onda_pos.txt");
    
    vOndaPos << " #          k                            F(Re)                        F(Im)                      norm(Re)                                  norm(Im)            " << std::endl; 
    for(i=0; i<N;i++){
        
        double k = i*delta_k;
        
        
        double real = ((out[i][0]*32*M_PI*rho)/(k*tamHistograma))+1.0;
        double im = -((out[i][1]*32*M_PI*rho)/(k*tamHistograma))+1.0;
        vOndaPos << std::setprecision(15) << i <<"         "<< k <<"                "<<out[i][0]<<"               "<<out[i][1]<<"                    "<< real <<"      " << im << std::endl; 
      //  j++;
        //en función de i, y en función de la frecuencia
    }
    
    
    j = N;
    std::ofstream vOndaNeg;
    vOndaNeg.open(dirSalida+"vectores_onda_neg.txt");
    for(i=N; i<MAX;i++){
        double k = -j*delta_k;
        //double w = 2*M_PI*f;
        j--;
        double val = out[i][0];
        vOndaNeg << std::setprecision(20) << i <<"     "<<k<<"         " << val <<"     " << out[i][1] << std::endl; 
    }
    
    
    fftw_destroy_plan(p);
    fftw_cleanup();
    
    
    
}

void sk(double *gdr, int tamHistograma, double delta_k, double rho, string dirSalida, double delta){
    
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
    
    
    
}

float distanciaAtomos(Atomo a1, Atomo a2){
    return a1.distancia(a2.getrx(),a2.getry(),a2.getrz());
}

void listaVecinos(vector<Atomo> atomos, int n_atomos, float r_min, double mitadCaja, double boxSize){
    //file:///D:/tesis/libros/computer_simultation_of_liquids.pdf pp. 162
    for(int i=0; i<n_atomos-1;i++){
        for(int j=0; j<n_atomos-1;j++){
            Atomo atomo1 = atomos[i];
            Atomo atomo2 = atomos[j];
            double distancia = distanciaAtomos(atomo1,atomo2);
            double distancias[3] = {0,0,0};
            if(i!=j){
                if(distancia<mitadCaja){
                    if(distancia<r_min){
                        std::cout << "atomo "<< i << "vecino " << j <<"Lo metemos a la lista" << std::endl;
                    }else{
                        atomo1.setPeriodics(atomo1.getPrx()+boxSize, atomo1.getPry()+boxSize, atomo1.getPrz()+boxSize);
                        distancia = distanciaAtomos(atomo1,atomo2);
                        distancia = distancia-mitadCaja;
                        if(distancia<r_min){
                            std::cout << "atomo imagen"<< i << "vecino " << j <<"Lo metemos a la lista" << std::endl;
                        }

                    }
                
            
                }
            }
        }
    }
}