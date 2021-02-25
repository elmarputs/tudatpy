/*    Copyright (c) 2010-2020, Delft University of Technology
*    All rights reserved
		*
		*    This file is part of the Tudat. Redistribution and use in source and
*    binary forms, with or without modification, are permitted exclusively
*    under the terms of the Modified BSD license. You should have received
*    a copy of the license with this file. If not, please or visit:
*    http://tudat.tudelft.nl/LICENSE.
*/

#include "expose_low_thrust.h"

#include <tudat/astro/low_thrust.h>
#include <tudat/astro/low_thrust/shape_based/baseFunctionsHodographicShaping.h>
#include <tudat/astro/low_thrust/shape_based/hodographicShaping.h>
#include <tudat/astro/low_thrust/shape_based/hodographicShapingOptimisationSetup.h>

namespace py = pybind11;
namespace tltt = tudat::low_thrust_trajectories;
namespace tsbm = tudat::shape_based_methods;
namespace trf = tudat::root_finders;

namespace tudatpy{

void expose_low_thrust( py::module &m )
{

	py::class_<tltt::LowThrustLeg, std::shared_ptr<tltt::LowThrustLeg>>(
			m, "LowThrustLeg", "<no_doc>" )
			.def("get_trajectory", py::overload_cast<std::vector<double>&>(&tltt::LowThrustLeg::getTrajectory),
			        py::arg("epochs"))
			.def("get_state", &tltt::LowThrustLeg::computeCurrentStateVector, py::arg("epoch"));

// lowThrustLegSettings.h

	py::enum_<tltt::LowThrustLegTypes>( m, "LowThrustLegTypes", "<no_doc>" )
			.value( "hodographic_shaping_leg", tltt::hodographic_shaping_leg )
			.value( "spherical_shaping_leg", tltt::spherical_shaping_leg )
					//.value("sims_flanagan_leg", tltt::sims_flanagan_leg)
					//.value("hybrid_method_leg", tltt::hybrid_method_leg)
			.export_values( );

	py::class_<tltt::LowThrustLegSettings, std::shared_ptr<tltt::LowThrustLegSettings>>(
			m, "LowThrustLegSettings", "<no_doc>" )
			.def( py::init<const tltt::LowThrustLegTypes>( ),
				  py::arg( "low_thrust_leg_type" ));

	py::class_<tltt::HodographicShapingLegSettings, std::shared_ptr<tltt::HodographicShapingLegSettings>,
			tltt::LowThrustLegSettings>(
			m, "HodographicShapingLegSettings", "<no_doc>" )
			.def( py::init<const int, const double,
						  std::vector<std::shared_ptr<tsbm::BaseFunctionHodographicShaping> > &,
						  std::vector<std::shared_ptr<tsbm::BaseFunctionHodographicShaping> > &,
						  std::vector<std::shared_ptr<tsbm::BaseFunctionHodographicShaping> > &,
						  const Eigen::VectorXd,
						  const Eigen::VectorXd,
						  const Eigen::VectorXd>( ),
				  py::arg( "number_of_revolutions" ),
				  py::arg( "central_body_gravitational_parameter" ),
				  py::arg( "radial_velocity_function_components" ),
				  py::arg( "normal_velocity_function_components" ),
				  py::arg( "axial_velocity_function_components" ),
				  py::arg( "free_coefficients_radial_velocity_function" ),
				  py::arg( "free_coefficients_normal_velocity_function" ),
				  py::arg( "free_coefficients_axial_velocity_function" ));

	py::class_<tltt::SphericalShapingLegSettings, std::shared_ptr<tltt::SphericalShapingLegSettings>,
			tltt::LowThrustLegSettings>(
			m, "SphericalShapingLegSettings", "<no_doc>" )
			.def( py::init<const int, const double, const double,
						  const std::shared_ptr<trf::RootFinderSettings> &,
						  const std::pair<double, double>>( ),
				  py::arg( "number_of_revolutions" ),
				  py::arg( "central_body_gravitational_parameter" ),
				  py::arg( "initial_value_free_coefficient" ),
				  py::arg( "rootfinder_settings" ),
				  py::arg( "bounds_free_coefficient" ) = std::make_pair( TUDAT_NAN, TUDAT_NAN ));

	m.def( "create_low_thrust_leg",
		   &tltt::createLowThrustLeg,
		   py::arg( "low_thrust_leg_settings" ),
		   py::arg( "state_at_departure" ),
		   py::arg( "state_at_arrival" ),
		   py::arg( "time_of_flight" ));

	m.def( "hodographic",
		   &tltt::hodographicShapingLegSettings,
		   py::arg( "number_of_revolutions" ),
		   py::arg( "central_body_gravitational_parameter" ),
		   py::arg( "radial_velocity_function_components" ),
		   py::arg( "normal_velocity_function_components" ),
		   py::arg( "axial_velocity_function_components" ),
		   py::arg( "free_coefficients_radial_velocity_function" ),
		   py::arg( "free_coefficients_normal_velocity_function" ),
		   py::arg( "free_coefficients_axial_velocity_function" ));

	m.def( "spherical",
		   &tltt::sphericalShapingLegSettings,
		   py::arg( "number_of_revolutions" ),
		   py::arg( "central_body_gravitational_parameter" ),
		   py::arg( "initial_value_free_coefficient" ),
		   py::arg( "root_finder_settings" ),
		   py::arg( "bounds_free_coefficient" ) = std::make_pair( TUDAT_NAN, TUDAT_NAN ));

// shapeBasedMethod.h

	py::class_<tsbm::ShapeBasedMethod, std::shared_ptr<tsbm::ShapeBasedMethod>, tltt::LowThrustLeg>
			( m, "ShapeBasedMethod", "<no_doc>" )
			.def( "get_thrust_force_profile", &tsbm::ShapeBasedMethod::getThrustForceProfile,
				  py::arg( "epochs_vector" ),
				  py::arg( "thrust_profile" ),
				  py::arg( "specific_impulse_function" ),
				  py::arg( "integrator_Settings" ));

// sphericalShaping.h

	py::class_<tsbm::SphericalShaping, std::shared_ptr<tsbm::SphericalShaping>, tsbm::ShapeBasedMethod>(
			m, "SphericalShaping", "<no_doc>" )
			.def( "compute_delta_v", &tsbm::SphericalShaping::computeDeltaV )
			.def( "compute_time_of_flight", &tsbm::SphericalShaping::computeTimeOfFlight );

// hodographicShaping.h

	py::class_<tsbm::HodographicShaping, std::shared_ptr<tsbm::HodographicShaping>, tsbm::ShapeBasedMethod>(
			m, "HodographicShaping", "<no_doc>")
			.def(py::init<const Eigen::Vector6d&, const Eigen::Vector6d&, const double, const double,
					const int, const std::vector< std::shared_ptr< tsbm::BaseFunctionHodographicShaping > >&,
					const std::vector< std::shared_ptr< tsbm::BaseFunctionHodographicShaping > >&,
					const std::vector< std::shared_ptr< tsbm::BaseFunctionHodographicShaping > >&,
					const Eigen::VectorXd&,
					const Eigen::VectorXd&,
					const Eigen::VectorXd&,
					const double>())
			.def("compute_delta_v", &tsbm::HodographicShaping::computeDeltaV);

// hodographicShapingOptimisationSetup.h

	py::class_<tsbm::FixedTimeHodographicShapingOptimisationProblem, std::shared_ptr<tsbm::FixedTimeHodographicShapingOptimisationProblem>>(
			m, "FixedTimeHodographicShapingOptimisationProblem", "<no_doc>")
			.def(py::init<
						const Eigen::Vector6d&,
						const Eigen::Vector6d&,
						const double,
						const double,
						const int,
						const std::vector< std::shared_ptr< tsbm::BaseFunctionHodographicShaping > >&,
						const std::vector< std::shared_ptr< tsbm::BaseFunctionHodographicShaping > >&,
						const std::vector< std::shared_ptr< tsbm::BaseFunctionHodographicShaping > >&,
						const std::vector< std::vector< double > >&
						>())
			.def("fitness", &tsbm::FixedTimeHodographicShapingOptimisationProblem::fitness, py::arg("x"))
			.def("get_bounds", &tsbm::FixedTimeHodographicShapingOptimisationProblem::get_bounds)
			.def("get_nobj", &tsbm::FixedTimeHodographicShapingOptimisationProblem::get_nobj);

// baseFunctionHodographicShaping.h

	py::class_<tsbm::BaseFunctionHodographicShaping, std::shared_ptr<tsbm::BaseFunctionHodographicShaping>>(
			m, "BaseFunctionHodographicShaping", "<no_doc>");

	py::class_<tsbm::ConstantFunctionHodographicShaping, std::shared_ptr<tsbm::ConstantFunctionHodographicShaping>,
	        tsbm::BaseFunctionHodographicShaping>(m, "ConstantFunctionHodographicShaping", "<no_doc>")
	        .def(py::init<>())
	        .def("evaluate_function", &tsbm::ConstantFunctionHodographicShaping::evaluateFunction, py::arg("indepent_variable"))
			.def("evaluate_derivative", &tsbm::ConstantFunctionHodographicShaping::evaluateDerivative, py::arg("indepent_variable"))
			.def("evaluate_integral", &tsbm::ConstantFunctionHodographicShaping::evaluateIntegral, py::arg("indepent_variable"))
			.def(py::pickle(
					[](const tsbm::ConstantFunctionHodographicShaping &c)
					{
						return py::make_tuple();
					},
					[](py::tuple t)
					{
						return tsbm::ConstantFunctionHodographicShaping();
					}
					));

	py::class_<tsbm::ScaledPowerFunctionHodographicShaping, std::shared_ptr<tsbm::ScaledPowerFunctionHodographicShaping>,
			tsbm::BaseFunctionHodographicShaping>(m, "ScaledPowerFunctionHodographicShaping", "<no_doc>")
			.def(py::init<double, double>())
			.def("evaluate_function", &tsbm::ScaledPowerFunctionHodographicShaping::evaluateFunction, py::arg("indepent_variable"))
			.def("evaluate_derivative", &tsbm::ScaledPowerFunctionHodographicShaping::evaluateDerivative, py::arg("indepent_variable"))
			.def("evaluate_integral", &tsbm::ScaledPowerFunctionHodographicShaping::evaluateIntegral, py::arg("indepent_variable"))
			.def(py::pickle(
					[](const tsbm::ScaledPowerFunctionHodographicShaping &s)
					{
						return py::make_tuple(s.exponent_, s.scaleFactor_);
					},
					[](py::tuple t)
					{
						if(t.size() != 2)
						{
							throw std::runtime_error("Invalid state!");
						}

						return tsbm::ScaledPowerFunctionHodographicShaping(t[0].cast<double>(), t[1].cast<double>());
					}
			));

	py::class_<tsbm::SineFunctionHodographicShaping, std::shared_ptr<tsbm::SineFunctionHodographicShaping>,
			tsbm::BaseFunctionHodographicShaping>(m, "SineFunctionHodographicShaping", "<no_doc>")
			.def(py::init<double>())
			.def("evaluate_function", &tsbm::SineFunctionHodographicShaping::evaluateFunction, py::arg("indepent_variable"))
			.def("evaluate_derivative", &tsbm::SineFunctionHodographicShaping::evaluateDerivative, py::arg("indepent_variable"))
			.def("evaluate_integral", &tsbm::SineFunctionHodographicShaping::evaluateIntegral, py::arg("indepent_variable"))
			.def(py::pickle(
					[](const tsbm::SineFunctionHodographicShaping &s)
					{
						return py::make_tuple(s.frequency_);
					},
					[](py::tuple t)
					{
						return tsbm::SineFunctionHodographicShaping(t[0].cast<double>());
					}
			));

	py::class_<tsbm::CosineFunctionHodographicShaping, std::shared_ptr<tsbm::CosineFunctionHodographicShaping>,
			tsbm::BaseFunctionHodographicShaping>(m, "CosineFunctionHodographicShaping", "<no_doc>")
			.def(py::init<double>())
			.def("evaluate_function", &tsbm::CosineFunctionHodographicShaping::evaluateFunction, py::arg("indepent_variable"))
			.def("evaluate_derivative", &tsbm::CosineFunctionHodographicShaping::evaluateDerivative, py::arg("indepent_variable"))
			.def("evaluate_integral", &tsbm::CosineFunctionHodographicShaping::evaluateIntegral, py::arg("indepent_variable"))
			.def(py::pickle(
					[](const tsbm::CosineFunctionHodographicShaping &s)
					{
						return py::make_tuple(s.frequency_);
					},
					[](py::tuple t)
					{
						return tsbm::CosineFunctionHodographicShaping(t[0].cast<double>());
					}
			));

	py::class_<tsbm::ScaledPowerSineFunctionHodographicShaping, std::shared_ptr<tsbm::ScaledPowerSineFunctionHodographicShaping>,
			tsbm::BaseFunctionHodographicShaping>(m, "ScaledPowerSineFunctionHodographicShaping", "<no_doc>")
			.def(py::init<double, double, double>())
			.def("evaluate_function", &tsbm::ScaledPowerSineFunctionHodographicShaping::evaluateFunction, py::arg("indepent_variable"))
			.def("evaluate_derivative", &tsbm::ScaledPowerSineFunctionHodographicShaping::evaluateDerivative, py::arg("indepent_variable"))
			.def("evaluate_integral", &tsbm::ScaledPowerSineFunctionHodographicShaping::evaluateIntegral, py::arg("indepent_variable"))
			.def(py::pickle(
					[](const tsbm::ScaledPowerSineFunctionHodographicShaping &s)
					{
						return py::make_tuple(s.exponentPowerFunction_, s.frequencySineFunction_, s.scaleFactor_);
					},
					[](py::tuple t)
					{
						return tsbm::ScaledPowerSineFunctionHodographicShaping(t[0].cast<double>(), t[1].cast<double>(), t[2].cast<double>());
					}
			));

	py::class_<tsbm::ScaledPowerCosineFunctionHodographicShaping, std::shared_ptr<tsbm::ScaledPowerCosineFunctionHodographicShaping>,
			tsbm::BaseFunctionHodographicShaping>(m, "ScaledPowerCosineFunctionHodographicShaping", "<no_doc>")
			.def(py::init<double, double, double>())
			.def("evaluate_function", &tsbm::ScaledPowerCosineFunctionHodographicShaping::evaluateFunction, py::arg("indepent_variable"))
			.def("evaluate_derivative", &tsbm::ScaledPowerCosineFunctionHodographicShaping::evaluateDerivative, py::arg("indepent_variable"))
			.def("evaluate_integral", &tsbm::ScaledPowerCosineFunctionHodographicShaping::evaluateIntegral, py::arg("indepent_variable"))
			.def(py::pickle(
					[](const tsbm::ScaledPowerCosineFunctionHodographicShaping &s)
					{
						return py::make_tuple(s.exponentPowerFunction_, s.frequencyCosineFunction_, s.scaleFactor_);
					},
					[](py::tuple t)
					{
						return tsbm::ScaledPowerCosineFunctionHodographicShaping(t[0].cast<double>(), t[1].cast<double>(), t[2].cast<double>());
					}
			));

	m.def( "hodograph_constant",
		   &tsbm::hodographConstant );

	m.def( "hodograph_sine",
		   &tsbm::hodographSine,
		   py::arg( "frequency" ));

	m.def( "hodograph_cosine",
		   &tsbm::hodographCosine,
		   py::arg( "frequency" ));

	m.def( "hodograph_exponential",
		   &tsbm::hodographExponential,
		   py::arg( "exponent" ));

	m.def( "hodograph_scaled_exponential",
		   &tsbm::hodographScaledExponential,
		   py::arg( "exponent" ),
		   py::arg( "scale_factor" ));

	m.def( "hodograph_exponential_sine",
		   &tsbm::hodographExponentialSine,
		   py::arg( "exponent" ),
		   py::arg( "frequency" ));

	m.def( "hodograph_scaled_exponential_sine",
		   &tsbm::hodographScaledExponentialSine,
		   py::arg( "exponent" ),
		   py::arg( "frequency" ),
		   py::arg( "scale_factor" ));

	m.def( "hodograph_exponential_cosine",
		   &tsbm::hodographExponentialCosine,
		   py::arg( "exponent" ),
		   py::arg( "frequency" ));

	m.def( "hodograph_scaled_exponential_cosine",
		   &tsbm::hodographScaledExponentialCosine,
		   py::arg( "exponent" ),
		   py::arg( "frequency" ),
		   py::arg( "scale_factor" ));

	m.def( "hodograph_power",
		   &tsbm::hodographPower,
		   py::arg( "exponent" ));

	m.def( "hodograph_scaled_power",
		   &tsbm::hodographScaledPower,
		   py::arg( "exponent" ),
		   py::arg( "scale_factor" ));


	m.def( "hodograph_power_sine",
		   &tsbm::hodographPowerSine,
		   py::arg( "exponent" ),
		   py::arg( "frequency" ));

	m.def( "hodograph_scaled_power_sine",
		   &tsbm::hodographScaledPowerSine,
		   py::arg( "exponent" ),
		   py::arg( "frequency" ),
		   py::arg( "scale_factor" ));

	m.def( "hodograph_power_cosine",
		   &tsbm::hodographPowerCosine,
		   py::arg( "exponent" ),
		   py::arg( "frequency" ));

	m.def( "hodograph_scaled_power_cosine",
		   &tsbm::hodographScaledPowerCosine,
		   py::arg( "exponent" ),
		   py::arg( "frequency" ),
		   py::arg( "scale_factor" ));

}
}
