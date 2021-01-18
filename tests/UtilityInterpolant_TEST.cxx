#include "gtest/gtest.h"

#include "PROPOSAL/propagation_utility/PropagationUtility.h"

using namespace PROPOSAL;

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

/* TEST(Calculate, Forward){ */
/*     auto integrand = [](double x)->double {return -1/x;}; */
/*     double lower_lim = 100; */

/*     auto integral = UtilityIntegral(integrand, lower_lim); */

/*     auto func = [&integral, lower_lim](double energy) { */
/*         return integral.Calculate(lower_lim, energy); */
/*     }; */
/*     Interpolant1DBuilder::Definition interpol_def; */
/*     interpol_def.function1d = func; */
/*     interpol_def.nodes = 200.; */
/*     interpol_def.xmin = 100.; */
/*     interpol_def.xmax = 1e14; */
/*     interpol_def.romberg = 5; */

/*     auto interpolant = UtilityInterpolant(integrand, lower_lim); */
/*     interpolant.BuildTables("unittest_interpolant", 4325465, interpol_def, false); */

/*     for (double logE_i = 2; logE_i < 14; logE_i+=1e-1) { */
/*         for (double logE_f = 2; logE_f < 14; logE_f+=1e-1){ */
/*             double E_i = std::pow(10, logE_i); */
/*             double E_f = std::pow(10, logE_f); */
/*             double analytical = std::log(E_i) - std::log(E_f); */
/*             if (E_i < E_f) { */
/*                 #ifndef NDEBUG */
/*                 EXPECT_DEATH(interpolant.Calculate(E_i, E_f), ""); */
/*                 #endif */
/*             } else { */
/*                 EXPECT_NEAR(interpolant.Calculate(E_i, E_f), analytical, analytical*1e-5); */
/*             } */
/*         } */
/*     } */
/* } */

/* TEST(Calculate, Reverse){ */
/*     auto integrand = [](double x)->double {return -1/x;}; */
/*     double lower_lim = 100; */

/*     auto integral = UtilityIntegral(integrand, lower_lim); */

/*     auto func = [&integral, lower_lim](double energy) { */
/*         return integral.Calculate(lower_lim, energy); */
/*     }; */
/*     Interpolant1DBuilder::Definition interpol_def; */
/*     interpol_def.function1d = func; */
/*     interpol_def.nodes = 200.; */
/*     interpol_def.xmin = 100.; */
/*     interpol_def.xmax = 1e14; */
/*     interpol_def.romberg = 5; */

/*     auto interpolant = UtilityInterpolant(integrand, lower_lim); */
/*     interpolant.BuildTables("unittest_interpolant", 4325465, interpol_def, true); */

/*     for (double logE_i = 2; logE_i < 14; logE_i+=1e-1) { */
/*         for (double logE_f = 2; logE_f < 14; logE_f+=1e-1){ */
/*             double E_i = std::pow(10, logE_i); */
/*             double E_f = std::pow(10, logE_f); */
/*             double analytical = std::log(E_i) - std::log(E_f); */
/*             if (E_i < E_f) { */
/*                 #ifndef NDEBUG */
/*                 EXPECT_DEATH(interpolant.Calculate(E_i, E_f), ""); */
/*                 #endif */
/*             } else { */
/*                 EXPECT_NEAR(interpolant.Calculate(E_i, E_f), analytical, analytical*1e-5); */
/*             } */
/*         } */
/*     } */
/* } */

/* TEST(Calculate, ApproximationCase){ */
/*     auto integrand = [](double x)->double {return -1/x;}; */
/*     double lower_lim = 100; */

/*     auto integral = UtilityIntegral(integrand, lower_lim); */

/*     auto func = [&integral, lower_lim](double energy) { */
/*         return integral.Calculate(lower_lim, energy); */
/*     }; */
/*     Interpolant1DBuilder::Definition interpol_def; */
/*     interpol_def.function1d = func; */
/*     interpol_def.nodes = 200.; */
/*     interpol_def.xmin = 0.; */
/*     interpol_def.xmax = 1e14; */
/*     interpol_def.romberg = 5; */

/*     auto interpolant = UtilityInterpolant(integrand, lower_lim); */
/*     interpolant.BuildTables("unittest_interpolant", 4325465, interpol_def, false); */

/*     double E_i = 1e7; */
/*     double E_f = 9999999; */
/*     EXPECT_TRUE(E_i - E_f < E_i * IPREC); //This is the condition we want to check in this UnitTest */
/*     double analytical = std::log(E_i) - std::log(E_f); */
/*     EXPECT_NEAR(interpolant.Calculate(E_i, E_f), analytical, analytical*1e-5); */
/* } */

/* TEST(Calculate, ApproximationCaseReverse){ */
/*     auto integrand = [](double x)->double {return -1/x;}; */
/*     double lower_lim = 100; */

/*     auto integral = UtilityIntegral(integrand, lower_lim); */

/*     auto func = [&integral, lower_lim](double energy) { */
/*         return integral.Calculate(lower_lim, energy); */
/*     }; */
/*     Interpolant1DBuilder::Definition interpol_def; */
/*     interpol_def.function1d = func; */
/*     interpol_def.nodes = 200.; */
/*     interpol_def.xmin = 0.; */
/*     interpol_def.xmax = 1e14; */
/*     interpol_def.romberg = 5; */

/*     auto interpolant = UtilityInterpolant(integrand, lower_lim); */
/*     interpolant.BuildTables("unittest_interpolant", 4325465, interpol_def, true); */

/*     double E_i = 1e7; */
/*     double E_f = 9999999; */
/*     EXPECT_TRUE(E_i - E_f < E_i * IPREC); //This is the condition we want to check in this UnitTest */
/*     double analytical = std::log(E_i) - std::log(E_f); */
/*     EXPECT_NEAR(interpolant.Calculate(E_i, E_f), analytical, analytical*1e-5); */
/* } */

/* TEST(GetUpperLimit, Forward){ */
/*     auto integrand = [](double x)->double {return -1/x;}; */
/*     double lower_lim = 100; */

/*     auto integral = UtilityIntegral(integrand, lower_lim); */

/*     auto func = [&integral, lower_lim](double energy) { */
/*         return integral.Calculate(lower_lim, energy); */
/*     }; */
/*     Interpolant1DBuilder::Definition interpol_def; */
/*     interpol_def.function1d = func; */
/*     interpol_def.nodes = 200.; */
/*     interpol_def.xmin = 100.; */
/*     interpol_def.xmax = 1e14; */
/*     interpol_def.romberg = 5; */

/*     auto interpolant = UtilityInterpolant(integrand, lower_lim); */
/*     interpolant.BuildTables("unittest_interpolant", 4325465, interpol_def, false); */

/*     for (double loglower = 2; loglower < 7; loglower+=1e-1) { */
/*         for (double logxi = -3; logxi < 3; logxi+=1e-1) { */
/*             double lower = std::pow(10, loglower); */
/*             double xi = std::pow(10, logxi); */
/*             double analytical_upper = lower * std::exp(-xi); */
/*             if (analytical_upper < 100) { */
/*                 #ifndef NDEBUG */
/*                 EXPECT_DEATH(interpolant.GetUpperLimit(lower, xi), ""); */
/*                 #endif */
/*             } else { */
/*                 EXPECT_NEAR(interpolant.GetUpperLimit(lower, xi), analytical_upper, analytical_upper*1e-5); */
/*             } */
/*         } */
/*     } */
/* } */

/* TEST(GetUpperLimit, Reverse){ */
/*     auto integrand = [](double x)->double {return -1/x;}; */
/*     double lower_lim = 100; */

/*     auto integral = UtilityIntegral(integrand, lower_lim); */

/*     auto func = [&integral, lower_lim](double energy) { */
/*         return integral.Calculate(lower_lim, energy); */
/*     }; */
/*     Interpolant1DBuilder::Definition interpol_def; */
/*     interpol_def.function1d = func; */
/*     interpol_def.nodes = 200.; */
/*     interpol_def.xmin = 100.; */
/*     interpol_def.xmax = 1e14; */
/*     interpol_def.romberg = 5; */

/*     auto interpolant = UtilityInterpolant(integrand, lower_lim); */
/*     interpolant.BuildTables("unittest_interpolant", 4325465, interpol_def, true); */

/*     for (double loglower = 2; loglower < 7; loglower+=1e-1) { */
/*         for (double logxi = -3; logxi < 3; logxi+=1e-1) { */
/*             double lower = std::pow(10, loglower); */
/*             double xi = std::pow(10, logxi); */
/*             double analytical_upper = lower * std::exp(-xi); */
/*             if (analytical_upper < 100) { */
/*                 continue; // no debug here */
/*             } else { */
/*                 EXPECT_NEAR(interpolant.GetUpperLimit(lower, xi), analytical_upper, analytical_upper*1e-5); */
/*             } */
/*         } */
/*     } */
/* } */

/* TEST(PlausabilityChecks, Forward){ */
/*     auto integrand = [](double x)->double {return -1/x;}; */
/*     double lower_lim = 100; */

/*     auto integral = UtilityIntegral(integrand, lower_lim); */

/*     auto func = [&integral, lower_lim](double energy) { */
/*         return integral.Calculate(lower_lim, energy); */
/*     }; */
/*     Interpolant1DBuilder::Definition interpol_def; */
/*     interpol_def.function1d = func; */
/*     interpol_def.nodes = 200.; */
/*     interpol_def.xmin = 100.; */
/*     interpol_def.xmax = 1e14; */
/*     interpol_def.romberg = 5; */

/*     auto interpolant = UtilityInterpolant(integrand, lower_lim); */
/*     interpolant.BuildTables("unittest_interpolant", 4325465, interpol_def, false); */

/*     auto energies = std::array<double, 6>{1e3, 1e5, 1e7, 1e9, 1e11, 1e13}; */
/*     for(double logE_i = 2; logE_i < 13; logE_i +=1e-1){ */
/*         for(double logE_f = 2; logE_f < 13; logE_f +=1e-1){ */
/*             double E_i = std::pow(10., logE_i); */
/*             double E_f = std::pow(10., logE_f); */
/*             if(E_i < E_f){ */
/*                 continue; */
/*             } */
/*             double xi = interpolant.Calculate(E_i, E_f); */
/*             double upper = interpolant.GetUpperLimit(E_i, xi); */
/*             EXPECT_NEAR(upper, E_f, E_f * 1e-5); */
/*         } */
/*     } */
/* } */

/* TEST(PlausabilityChecks, Reverse){ */
/*     auto integrand = [](double x)->double {return -1/x;}; */
/*     double lower_lim = 100; */

/*     auto integral = UtilityIntegral(integrand, lower_lim); */

/*     auto func = [&integral, lower_lim](double energy) { */
/*         return integral.Calculate(lower_lim, energy); */
/*     }; */
/*     Interpolant1DBuilder::Definition interpol_def; */
/*     interpol_def.function1d = func; */
/*     interpol_def.nodes = 200.; */
/*     interpol_def.xmin = 100.; */
/*     interpol_def.xmax = 1e14; */
/*     interpol_def.romberg = 5; */

/*     auto interpolant = UtilityInterpolant(integrand, lower_lim); */
/*     interpolant.BuildTables("unittest_interpolant", 4325465, interpol_def, true); */

/*     auto energies = std::array<double, 6>{1e3, 1e5, 1e7, 1e9, 1e11, 1e13}; */
/*     for(double logE_i = 2; logE_i < 13; logE_i +=1e-1){ */
/*         for(double logE_f = 2; logE_f < 13; logE_f +=1e-1){ */
/*             double E_i = std::pow(10., logE_i); */
/*             double E_f = std::pow(10., logE_f); */
/*             if(E_i < E_f){ */
/*                 continue; */
/*             } */
/*             double xi = interpolant.Calculate(E_i, E_f); */
/*             double upper = interpolant.GetUpperLimit(E_i, xi); */
/*             EXPECT_NEAR(upper, E_f, E_f * 1e-5); */
/*         } */
/*     } */
/* } */

/* TEST(IntegrationFromMax, Calculate){ */
/*     auto integrand = [](double x)->double {return -1/x;}; */
/*     double lower_lim = 100; */

/*     auto integral = UtilityIntegral(integrand, lower_lim); */
/*     double xmax = 1e14; */

/*     auto func = [&integral, lower_lim](double energy) { */
/*         return integral.Calculate(lower_lim, energy); */
/*     }; */
/*     Interpolant1DBuilder::Definition interpol_def; */
/*     interpol_def.function1d = func; */
/*     interpol_def.nodes = 200.; */
/*     interpol_def.xmin = 0.; */
/*     interpol_def.xmax = xmax; */
/*     interpol_def.romberg = 5; */

/*     auto interpolant = UtilityInterpolant(integrand, lower_lim); */
/*     interpolant.BuildTables("unittest_interpolant", 4325465, interpol_def); */

/*     double E_i = 1e5; */
/*     double E_f = 1e4; */
/*     double analytical = std::log(E_i) - std::log(E_f); */
/*     EXPECT_NEAR(interpolant.Calculate(E_i, E_f), analytical, analytical*1e-5); */
/* } */

/* TEST(IntegrationFromMax, ApproximationCase){ */
/*     auto integrand = [](double x)->double {return -1/x;}; */
/*     double lower_lim = 100; */

/*     auto integral = UtilityIntegral(integrand, lower_lim); */
/*     double xmax = 1e14; */

/*     auto func = [&integral, lower_lim](double energy) { */
/*         return integral.Calculate(lower_lim, energy); */
/*     }; */
/*     Interpolant1DBuilder::Definition interpol_def; */
/*     interpol_def.function1d = func; */
/*     interpol_def.nodes = 200.; */
/*     interpol_def.xmin = 0.; */
/*     interpol_def.xmax = xmax; */
/*     interpol_def.romberg = 5; */

/*     auto interpolant = UtilityInterpolant(integrand, lower_lim); */
/*     interpolant.BuildTables("unittest_interpolant", 4325465, interpol_def); */

/*     double E_i = 1e7; */
/*     double E_f = 9999999; */
/*     EXPECT_TRUE(E_i - E_f < E_i * IPREC); //This is the condition we want to check in this UnitTest */
/*     double analytical = std::log(E_i) - std::log(E_f); */
/*     EXPECT_NEAR(interpolant.Calculate(E_i, E_f), analytical, analytical*1e-5); */
/* } */

/* TEST(IntegrationFromMax, GetUpperLimit){ */
/*     auto integrand = [](double x)->double {return -1/x;}; */
/*     double lower_lim = 100; */

/*     auto integral = UtilityIntegral(integrand, lower_lim); */
/*     double xmax = 1e14; */

/*     auto func = [&integral, lower_lim](double energy) { */
/*         return integral.Calculate(lower_lim, energy); */
/*     }; */
/*     Interpolant1DBuilder::Definition interpol_def; */
/*     interpol_def.function1d = func; */
/*     interpol_def.nodes = 200.; */
/*     interpol_def.xmin = 0.; */
/*     interpol_def.xmax = xmax; */
/*     interpol_def.romberg = 5; */

/*     auto interpolant = UtilityInterpolant(integrand, lower_lim); */
/*     interpolant.BuildTables("unittest_interpolant", 4325465, interpol_def); */

/*     double E_i = 1e5; */
/*     double xi = 0.5; */
/*     double analytical_upper = E_i * std::exp(-xi); */
/*     EXPECT_NEAR(interpolant.GetUpperLimit(E_i, xi), analytical_upper, analytical_upper*1e-5); */
/* } */

/* TEST(IntegrationFromMax, PlausabilityChecks){ */
/*     auto integrand = [](double x)->double {return -1/x;}; */
/*     double lower_lim = 100; */

/*     auto integral = UtilityIntegral(integrand, lower_lim); */
/*     double xmax = 1e14; */

/*     auto func = [&integral, lower_lim](double energy) { */
/*         return integral.Calculate(lower_lim, energy); */
/*     }; */
/*     Interpolant1DBuilder::Definition interpol_def; */
/*     interpol_def.function1d = func; */
/*     interpol_def.nodes = 200.; */
/*     interpol_def.xmin = 0.; */
/*     interpol_def.xmax = xmax; */
/*     interpol_def.romberg = 5; */

/*     auto interpolant = UtilityInterpolant(integrand, lower_lim); */
/*     interpolant.BuildTables("unittest_interpolant", 4325465, interpol_def); */

/*     auto energies = std::array<double, 6>{1e3, 1e5, 1e7, 1e9, 1e11, 1e13}; */
/*     for(auto E_i : energies){ */
/*         for(auto E_f : energies){ */
/*             if(E_i < E_f){ */
/*                 continue; */
/*             } */
/*             double xi = interpolant.Calculate(E_i, E_f); */
/*             double upper = interpolant.GetUpperLimit(E_i, xi); */
/*             EXPECT_NEAR(upper, E_f, E_f * 1e-5); */
/*         } */
/*     } */
/* } */

/* TEST(PositiveFunction, Calculate){ */
/*     auto integrand = [](double x)->double {return 1/x;}; */
/*     double lower_lim = 100; */

/*     auto integral = UtilityIntegral(integrand, lower_lim); */

/*     auto func = [&integral, lower_lim](double energy) { */
/*         return integral.Calculate(lower_lim, energy); */
/*     }; */
/*     Interpolant1DBuilder::Definition interpol_def; */
/*     interpol_def.function1d = func; */
/*     interpol_def.nodes = 200.; */
/*     interpol_def.xmin = 0.; */
/*     interpol_def.xmax = 1e14; */
/*     interpol_def.romberg = 5; */

/*     auto interpolant = UtilityInterpolant(integrand, lower_lim); */
/*     interpolant.BuildTables("unittest_interpolant", 4325465, interpol_def); */

/*     double E_i = 1e5; */
/*     double E_f = 1e4; */
/*     double analytical = std::log(E_f) - std::log(E_i); */
/*     EXPECT_NEAR(interpolant.Calculate(E_i, E_f), analytical, std::abs(analytical*1e-5)); */
/* } */

/* TEST(PositiveFunction, ApproximationCase){ */
/*     auto integrand = [](double x)->double {return 1/x;}; */
/*     double lower_lim = 100; */

/*     auto integral = UtilityIntegral(integrand, lower_lim); */

/*     auto func = [&integral, lower_lim](double energy) { */
/*         return integral.Calculate(lower_lim, energy); */
/*     }; */
/*     Interpolant1DBuilder::Definition interpol_def; */
/*     interpol_def.function1d = func; */
/*     interpol_def.nodes = 200.; */
/*     interpol_def.xmin = 0.; */
/*     interpol_def.xmax = 1e14; */
/*     interpol_def.romberg = 5; */

/*     auto interpolant = UtilityInterpolant(integrand, lower_lim); */
/*     interpolant.BuildTables("unittest_interpolant", 4325465, interpol_def); */

/*     double E_i = 1e7; */
/*     double E_f = 9999999; */
/*     EXPECT_TRUE(E_i - E_f < E_i * IPREC); //This is the condition we want to check in this UnitTest */
/*     double analytical = std::log(E_f) - std::log(E_i); */
/*     EXPECT_NEAR(interpolant.Calculate(E_i, E_f), analytical, std::abs(analytical*1e-5)); */
/* } */

/* TEST(PositiveFunction, GetUpperLimit){ */
/*     auto integrand = [](double x)->double {return 1/x;}; */
/*     double lower_lim = 100; */

/*     auto integral = UtilityIntegral(integrand, lower_lim); */

/*     auto func = [&integral, lower_lim](double energy) { */
/*         return integral.Calculate(lower_lim, energy); */
/*     }; */
/*     Interpolant1DBuilder::Definition interpol_def; */
/*     interpol_def.function1d = func; */
/*     interpol_def.nodes = 200.; */
/*     interpol_def.xmin = 0.; */
/*     interpol_def.xmax = 1e14; */
/*     interpol_def.romberg = 5; */

/*     auto interpolant = UtilityInterpolant(integrand, lower_lim); */
/*     interpolant.BuildTables("unittest_interpolant", 4325465, interpol_def); */

/*     double E_i = 1e5; */
/*     double xi = 0.5; */
/*     double analytical_upper = E_i * std::exp(xi); */
/*     EXPECT_NEAR(interpolant.GetUpperLimit(E_i, xi), analytical_upper, analytical_upper*1e-5); */
/* } */
