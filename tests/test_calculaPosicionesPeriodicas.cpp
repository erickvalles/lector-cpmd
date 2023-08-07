// tests/test_calculaPosicionesPeriodicas.cpp
#define CATCH_CONFIG_MAIN
#include "../catch.hpp"
//#include <catch2/catch_test_macros.hpp>
//#include <catch2/catch_approx.hpp>
#include "../Atomo.h" // Asegúrate de incluir correctamente Atomo.h según su ubicación
#include "../utils/utils.h"

TEST_CASE("Prueba de calculaPosicionesPeriodicas", "[calculaPosicionesPeriodicas]") {
    double posPeriodicas[3];
    // Llama a la función que deseas probar
    calculaPosicionesPeriodicas(posPeriodicas, 5.0, 6.0, 7.0, 12.0, 6.0);

    // Realiza las aserciones para verificar que la función se comporta correctamente
    REQUIRE(posPeriodicas[0] == Approx(5.0));
    /*REQUIRE(posPeriodicas[1] == Approx(6.0));
    REQUIRE(posPeriodicas[2] == Approx(7.0));*/
}

TEST_CASE("Prueba de evalúa caja", "[evaluaCajaRecursivo]"){


    double val = evaluaCajaRecursivo(6, 10, 5);
    /*double val2 = evaluaCajaRecursivo(7, 10, 5);
    double val3 = evaluaCajaRecursivo(-7, 10, 5);
    double val4 = evaluaCajaRecursivo(5, 10, 5);*/
    
    REQUIRE(val == Approx(-4.0));
    /*REQUIRE(val2 == Approx(-3.0));
    REQUIRE(val3 == Approx(3.0));
    REQUIRE(val3 == Approx(3.0));*/


}

TEST_CASE("Prueba verifica componente para imagen", "[verificaComponenteParaImagen]"){
    double val1 = verificaComponenteParaImagen(7,10,5);
    /*double val2 = verificaComponenteParaImagen(5,10,5);
    double val3 = verificaComponenteParaImagen(25,10,5);
    double val4 = verificaComponenteParaImagen(128.3,10,5);
    double val5 = verificaComponenteParaImagen(12878.38,10,5);*/
    
    REQUIRE(val1 == Approx(-3.0));
   /** REQUIRE(val2 == Approx(5.0));
    REQUIRE(val3 <= 5);
    REQUIRE(val3 >= -5);
    REQUIRE(val4 <= 5);
    REQUIRE(val4 >= -5);
    REQUIRE(val5 <= 5);
    REQUIRE(val5 >= -5);
    std::cout << "periodics: " << val1 << "   " << val2 <<"       "<< val3 << val4 << std::endl;*/
}

TEST_CASE("Prueba verifica si son vecinos", "[verificaVecindad]"){
    Atomo a0;
    Atomo a2;
    Atomo a3;
    Atomo a4;
    Atomo a5;
    Atomo a6;
    Atomo a45;
    Atomo a91;
    Atomo a94;
    Atomo a144;
    a0.setPos(1.89707997578496,43.35546232056073,-8.70900464433794);
    a2.setPos(78.73831638262645,90.22621865308219,17.94512167752482);
    a3.setPos(-6.46184284141017,47.41217064288488,-33.12068023480352);
    a4.setPos(26.48719104222151,-4.40497856534153,24.47295435059297);
    a5.setPos(40.42484067905382,33.71852966781326,11.37825900295537);
    a6.setPos(-2.32591888641725,54.42694255211948,53.72243699934348);
    a45.setPos(24.07170520664087,20.70426205957773,34.71702962482647);
    a94.setPos(-10.24779514520402,7.69050692061653,22.02463411676478);
    a91.setPos(7.42260755766238,26.80624524577102,43.39593342465682);
    a144.setPos(30.82944755047127,-3.25710163257071,30.59726560193713);
    
    double boxSize = 22.3797;
    
    double mitadCaja = boxSize/2;
    double r_min = 6;
    bool resultado1 = verificaVecindad(a0,a2,r_min,mitadCaja,boxSize);
    bool resultado2 = verificaVecindad(a0,a3,r_min,mitadCaja,boxSize);
    bool resultado3 = verificaVecindad(a0,a4,r_min,mitadCaja,boxSize);
    bool resultado4 = verificaVecindad(a0,a45,r_min,mitadCaja,boxSize);
    bool resultado5 = verificaVecindad(a0,a94,r_min,mitadCaja,boxSize);
    bool resultado6 = verificaVecindad(a0,a144,r_min,mitadCaja,boxSize);
    bool resultado7 = verificaVecindad(a94,a144,r_min,mitadCaja,boxSize);
    bool resultado8 = verificaVecindad(a2,a91,r_min,mitadCaja,boxSize);

    REQUIRE(resultado1 == false);
    REQUIRE(resultado2 == false);
    REQUIRE(resultado3 == false);
    REQUIRE(resultado4 == true);
    REQUIRE(resultado5 == false);
    REQUIRE(resultado6 == false);
    REQUIRE(resultado7 == true);
    REQUIRE(resultado8 == true);

}
