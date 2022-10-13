//
// Created by Erick on 15/06/2020.
//
#include "utils.h"
#include "cmath"
#include <iostream>
#include <fstream>
#include <iomanip>
#include </usr/include/fftw3.h>
#include <map>

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

            if(distancia<mitadCaja){//falta obtener las imágenes también de aquí, restando a la distancia el tamaño de la caja
                position = distancia/delta;
                (*histo)[position] += 1;

            }/*else{
                double imgx = verificaComponenteParaImagen(distancias[0],boxSize,mitadCaja);
                double imgy = verificaComponenteParaImagen(distancias[1],boxSize,mitadCaja);
                double imgz = verificaComponenteParaImagen(distancias[2],boxSize,mitadCaja);

                double dist_comp_img[3] = {imgx,imgy,imgz};

                distancia2 = calculaDistancia(dist_comp_img);

                distancia = sqrt(distancia2);
                if(distancia<mitadCaja){
                    position = distancia/delta;
                    (*histo)[position] += 1;
                }
            }*/

            //std::cout << "recorrer 300 tarda "<< tdist.GetMs() << std::endl;
        }
        //tc.getTime();
    }
    //std::cout << "por pasada no cumplen: "<< no_cumplen << std::endl;

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
    
    distancias[0]=verificaComponenteParaImagen(compx,boxSize, halfBox);
    distancias[1]=verificaComponenteParaImagen(compy,boxSize, halfBox);
    distancias[2]=verificaComponenteParaImagen(compz,boxSize, halfBox);
    
}



double verificaDistancia(double distancia, double boxSize, double halfBox){
    return distancia>=(halfBox)?distancia-boxSize:distancia;
/*    if(distancia>halfBox){
        return distancia-boxSize;
    }else{
        return distancia;
    }*/
}
double verificaComponenteParaImagen(double dif_componentes, double boxSize, double halfBox){//usar esta angulos
    
    if(dif_componentes > 0){
        return dif_componentes>(halfBox)?dif_componentes-boxSize:dif_componentes;
    }
    if(dif_componentes < 0){
        return dif_componentes<(-halfBox)?dif_componentes+boxSize:dif_componentes;
    }
    
    return dif_componentes;
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
double integral(vector <double> funcion_integrar, int n_divisiones){ //determinar quin será a y hasta dónde llegará b
    double tam_funcion = funcion_integrar.size();
    double a = 0.0;//el primer extremo en la función (valor en las x )
    double b = 180.0;//el valor en las x
    double h = (b - a)/n_divisiones;//(b-a)/n;
    
    double fa = funcion_integrar.front();//son las amplitudes inicial
    double fb = funcion_integrar.back();//amplitud final
    double sum = 0;
    
    for(double i = 1; i<tam_funcion;i++){
        
        sum+=funcion_integrar[i];
    }
    
    return h*(((fa+fb)/2) + sum);//*/
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

float distanciaAtomos(Atomo a1, Atomo a2){
    return a1.distancia(a2.getPrx(),a2.getPry(),a2.getPrz());
}

void listaVecinos(vector<Atomo> atomos, int n_atomos, float r_min, double mitadCaja, double boxSize, std::map<Atomo,vector<Atomo>> vecinos, vector<double> *histAngulos, double deltaAng, std::string trayectoria){
    //file:///D:/tesis/libros/computer_simultation_of_liquids.pdf pp. 162
    for(auto &aa1: atomos){
        for(auto aa2:atomos){
            Atomo atomo1 = aa1;
            Atomo atomo2 = aa2;
            
            if(atomo1.getId()!=atomo2.getId()){
                double distancias_comp[3] = {0,0,0};
                calcula_dist_componentes(distancias_comp,atomo1,atomo2,boxSize,mitadCaja);

                // std::cout << "Lectura   "<< j << std::endl;

                double distancia2 = calculaDistancia(distancias_comp);
            //
            //double distancia = sqrtf(distancia2);
                double distancia = sqrt(distancia2);
            
                    if(distancia<=mitadCaja){
                        if(distancia<r_min){//imágenes
                            //if((atomo1.getId()==0 && atomo2.getId()==94) || (atomo1.getId()==0 && atomo2.getId()==45)){
                               // std::cout << distancia << std::endl;
                            //}
                        
                            vecinos[atomo1].push_back(atomo2);
                        }

                    }else{
                        double imgx = verificaComponenteParaImagen(distancias_comp[0],boxSize,mitadCaja);
                        double imgy = verificaComponenteParaImagen(distancias_comp[1],boxSize,mitadCaja);
                        double imgz = verificaComponenteParaImagen(distancias_comp[2],boxSize,mitadCaja);

                        double dist_comp_img[3] = {imgx,imgy,imgz};

                        double otraDistancia2 = calculaDistancia(dist_comp_img);

                        double OtraDistancia = sqrt(otraDistancia2);
                         
                        if(OtraDistancia<=mitadCaja){
                            if(OtraDistancia<r_min){
                                //if((atomo1.getId()==0 && atomo2.getId()==94) || (atomo1.getId()==0 && atomo2.getId()==45)){
                                    //std::cout << distancia << " imagen" << std::endl;
                                //}
                            
                            vecinos[atomo1].push_back(atomo2);
                            }
                        }
                        
                    }
                }
            }
        }
    std::map<Atomo,vector<Atomo>>::iterator it;
    std::ofstream problematicos;
    problematicos.open("/home/erick/data/salidas/input_problematicas.txt");
    for(it=vecinos.begin(); it!=vecinos.end(); it++){
        Atomo a1 = it->first;
        //obtener los vecinos de a1
        vector<Atomo> vecinos = it->second;
        
        //std::cout << "Atomo " << it->first.getId() << " tiene " << it->second.size() << " vecinos" << std::endl;
        for (int i=0; i<vecinos.size(); i++){
            //iterar después del atomo actual
            Atomo vFijo = vecinos[i];
            vector<double> v1 = a1.vectorDifference(vFijo);
            
            for(int j=i+1; j<vecinos.size(); j++){
                
                Atomo a2 = vecinos[j];
                //std::cout << "i: "<< a1.getId() << "j: " << vFijo.getId() <<"k: "<< a2.getId()<< std::endl;
                vector<double> v2 = a1.vectorDifference(a2);
                problematicos << a1.getId() <<"  "<< vFijo.getId() <<" "<<a2.getId() << std::endl; 
                //calcular el producto punto de los vectores v1 y v2
                double productoPunto = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
                //calcular el modulo del vector v1
                double moduloV1 = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
                //calcular el modulo del vector v2
                double moduloV2 = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
                //calcular el angulo entre los vectores v1 y v2
                double anguloRad = acos(productoPunto/(moduloV1*moduloV2));
                //convertir el angulo a grados
                double anguloGrad = anguloRad*180/M_PI;
                if(anguloGrad < 12 && anguloGrad > 10){
                    std::cout << trayectoria << std::endl;
                    std::cout << "periodicas" << a1.getPrx() <<" " << a2.getPrx() << std::endl;
                    //problematicos << trayectoria << std::endl;                    
                }
                int pos = 0;
                
                pos = anguloGrad/deltaAng;
                //imprimir la posicion en la que quedará
                //std::cout << "Pos = "<< anguloGrad << "/" << deltaAng <<"="<<pos << std::endl;
                //incrementar esa posicion en el histograma	
                (*histAngulos)[pos]+=1;
                //imprimir cuantos elementos en esa posicion hay
                //std::cout << "Histograma["<<pos<<"]="<<(*histAngulos)[pos]<<std::endl;
            }
            
            
        }
        
    }
    problematicos.close();

}