
// unique_ptr<Interpolant> IonizInterpolant::build_dedx(
//     const ParticleDef& p_def, const Medium& medium)
// {
//     using ioniz = std::decay<typeid(medium)::value_type>::type;
//     std::unique_ptr<Ionization> param(new ioniz(parametrization_.get()));
//     IonizIntegral ioniz_integral(std::move(param), cuts_);
//     Interpolant1DBuilder::Definition interpol_def;
//     /* interpol_def.function1d = [this, &p_def, &medium](double energy) { */
//     /*     auto is_bhabha = parametrization_->name ==
//      * "IonizBergerSeltzerBhabha"; */
//     /*     auto is_moller = parametrization_->name ==
//      * "IonizBergerSeltzerMoller"; */
//     /*     if (is_bhabha || is_moller) */
//     /*         return parametrization_->FunctionToDEdxIntegral( */
//     /*             p_def, medium, energy, 0.); */
//     /*     auto physical_lim =
//      * reinterpret_cast<Ionization*>(parametrization_.get()) */
//     /*                             ->GetKinematicLimits(p_def, medium, energy);
//      */
//     /*     auto v_cut = GetEnergyCut(energy, physical_lim); */
//     /*     auto dEdx = [this, &p_def, &medium, energy](double v) { */
//     /*         return parametrization_->FunctionToDEdxIntegral( */
//     /*             p_def, medium, energy, v); */
//     /*     }; */
//     /*     return integral_.Integrate(physical_lim.vMin, v_cut, dEdx, 4); */
//     /* }; */
//     /* interpol_def.max = def.nodes_cross_section; */
//     interpol_def.xmin = parametrization_->GetLowerEnergyLim(p_def);
//     /* interpol_def.xmax = def.max_node_energy; */
//     /* interpol_def.romberg = def.order_of_interpolation; */
//     interpol_def.rational = true;
//     interpol_def.isLog = true;
//     /* interpol_def.rombergY = def.order_of_interpolation; */
//     interpol_def.logSubst = true;
//     return Helper::InitializeInterpolation("dEdx",
//         make_unique<Interpolant1DBuilder>(interpol_def), GetHash(p_def, medium),
//         def);
// }
