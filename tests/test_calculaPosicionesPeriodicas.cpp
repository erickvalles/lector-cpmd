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
    REQUIRE(posPeriodicas[1] == Approx(6.0));
    REQUIRE(posPeriodicas[2] == Approx(7.0));
}

TEST_CASE("Prueba de evalúa caja", "[evaluaCajaRecursivo]"){


    double val = evaluaCajaRecursivo(6, 10, 5);
    double val2 = evaluaCajaRecursivo(7, 10, 5);
    double val3 = evaluaCajaRecursivo(-7, 10, 5);
    double val4 = evaluaCajaRecursivo(5, 10, 5);
    
    REQUIRE(val == Approx(-4.0));
    REQUIRE(val2 == Approx(-3.0));
    REQUIRE(val3 == Approx(3.0));
    REQUIRE(val3 == Approx(3.0));


}
