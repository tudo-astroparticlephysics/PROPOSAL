
/* #pragma once */

/* #include "PROPOSAL/crosssection/CrossSection.h" */
/* #include "PROPOSAL/crosssection/parametrization/Parametrization.h" */

/* namespace PROPOSAL{ */

/*     class CrossSectionBuilder final : public CrossSection */
/*     { */
/*     public: */
/*         CrossSectionBuilder(std::string name, const ParticleDef &particle_def, std::shared_ptr<const Medium> medium = nullptr); */

/*         double CalculatedEdx(double energy) override; */
/*         double CalculatedE2dx(double energy) override; */
/*         double CalculatedNdx(double energy) override; */
/*         double CalculatedNdx(double energy, double rnd) override; */
/*         double CalculateStochasticLoss(double energy, double rnd1, double rnd2) override; */
/*         double CalculateCumulativeCrossSection(double energy, int component, double v) override; */

/*         void SetdEdx_function(Function func); */
/*         void SetdE2dx_function(Function func); */
/*         void SetdNdx_function(Function func); */
/*         void SetdNdx_rnd_function(std::function<double(double, double)> func); */
/*         void SetStochasticLoss_function(std::function<double(double, double, double)> func); */
/*         void SetCumulativeCrossSection_function(std::function<double(double, double, double)> func); */

/*         size_t GetHash() const override; */

/*     protected: */
/*         virtual bool compare(const CrossSection& cross) const override; */
/*         double CalculateStochasticLoss(double energy, double rnd1) override; */

/*     private: */
/*         std::string name; */
/*         Function dEdx_function; */
/*         Function dE2dx_function; */
/*         Function dNdx_function; */
/*         std::function<double(double, double)> dNdx_rnd_function; */
/*         std::function<double(double, double, double)> StochasticLoss_function; */
/*         std::function<double(double, double, double)> CumulativeCrossSection_function; */
/*     }; */

/*     std::ostream& operator<<(std::ostream&, PROPOSAL::CrossSection const&); */

/* } */
