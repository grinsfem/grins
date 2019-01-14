//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
// Copyright (C) 2010-2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-


#ifndef GRINS_SPECTROSCOPIC_TEST_BASE_H
#define GRINS_SPECTROSCOPIC_TEST_BASE_H

#ifdef GRINS_HAVE_CPPUNIT

#include "test_comm.h"
#include "grins_test_paths.h"

#include "absorption_coeff_testing.h"

// GRINS
#include "grins/math_constants.h"
#include "grins/grins_enums.h"
#include "grins/simulation_builder.h"
#include "grins/simulation.h"
#include "grins/variable_warehouse.h"
#include "grins/single_variable.h"
#include "grins/multicomponent_variable.h"
#include "grins/chemistry_builder.h"
#include "grins/antioch_chemistry.h"

// libMesh
#include "libmesh/getpot.h"

#include <tuple>

namespace GRINSTesting
{
  class SpectroscopicTestBase
  {
  protected:
    std::shared_ptr<GRINS::Simulation> _sim;
    std::shared_ptr<GetPot> _input;

    //! Run the test on a given input file and calculated answer
    void run_test(std::stringstream & input_stream, libMesh::Real calc_answer)
    {
      this->init_sim(input_stream);

      _sim->run();

      CPPUNIT_ASSERT_DOUBLES_EQUAL( calc_answer, _sim->get_qoi_value(0),libMesh::TOLERANCE  );
    }

    //! Initialize the GetPot and Simulation class objects
    void init_sim(std::stringstream & input_stream)
    {
      _input.reset(new GetPot(std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/spectroscopic_base.in"));
      _input->parse_input_stream(input_stream);

      const char * const argv = "unit_driver";
      GetPot empty_command_line( (const int)1,&argv );
      GRINS::SimulationBuilder sim_builder;

      _sim.reset( new GRINS::Simulation(*_input,
                                        empty_command_line,
                                        sim_builder,
                                        *TestCommWorld ) );
    }

    void param_deriv_test(std::stringstream & input_stream)
    {
      this->init_sim(input_stream);

      std::string material = "TestMaterial";
      std::string hitran_data = "./test_data/CO2_data.dat";
      std::string hitran_partition = "./test_data/CO2_partition_function.dat";
      libMesh::Real T_min = 290.0,
        T_max = 310.0,
        T_step = 0.01;
      std::shared_ptr<GRINS::HITRAN> hitran( new GRINS::HITRAN(hitran_data,hitran_partition,T_min,T_max,T_step) );

      std::string species = "CO2";
      libMesh::Real thermo_pressure = 5066.25;
      libMesh::Real nu_desired = 3682.7649;
      libMesh::Real nu_min = 3682.69;
      libMesh::Real nu_max = 3682.8;
      GRINS::ChemistryBuilder chem_builder;
      std::unique_ptr<GRINS::AntiochChemistry> chem_ptr;
      chem_builder.build_chemistry(*(_input.get()),material,chem_ptr);
      std::shared_ptr<GRINS::AntiochChemistry> chem(chem_ptr.release());
      std::shared_ptr<AbsorptionCoeffTesting<GRINS::AntiochChemistry> >
        absorb( new AbsorptionCoeffTesting<GRINS::AntiochChemistry>(chem,hitran,nu_min,nu_max,nu_desired,species,thermo_pressure) );

      libMesh::Real T = 300.0;
      libMesh::Real P = 5066.25;

      std::vector<libMesh::Real> Y = {0.0763662233, 0.9236337767};

      for (unsigned int i=0; i<33; ++i)
        {
          this->T_param_derivatives(absorb,T,P,Y,i);
          this->P_param_derivatives(absorb,T,P,Y,i);
          this->Y_param_derivatives(absorb,T,P,Y,i);
        }
    }

    void elem_qoi_derivative_test(std::stringstream & input_stream)
    {
      // run the Simulation once to ensure everything is initialized and populated
      this->init_sim(input_stream);
      _sim->run();

      GRINS::MultiphysicsSystem * system = _sim->get_multiphysics_system();
      GRINS::AssemblyContext * context = libMesh::cast_ptr<GRINS::AssemblyContext*>(system->build_context().release());
      system->init_context(*context);
      context->pre_fe_reinit(*system,system->get_mesh().elem_ptr(0));

      libMesh::QoISet qs;
      qs.add_index(0);

      libMesh::DifferentiableQoI * qoi = system->get_qoi();
      qoi->element_qoi_derivative(*context,qs);

      this->T_elem_derivative(system,context);
      this->P_elem_derivative(system,context);
      this->Y_elem_derivative(system,context);
    }

    void assemble_qoi_derivative_test(std::stringstream & input_stream)
    {
      // run the Simulation once to ensure everything is initialized and populated
      this->init_sim(input_stream);
      _sim->run();

      GRINS::MultiphysicsSystem * system = _sim->get_multiphysics_system();

      libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> > & solution = (*system).solution;

      libMesh::QoISet qs;
      qs.add_index(0);

      system->assemble_qoi_derivative(qs,false,false);
      libMesh::NumericVector<libMesh::Number> & derivs = system->get_adjoint_rhs();

      libMesh::Real delta = 1.0e-6;

      for (unsigned int dof=0; dof<solution->size(); ++dof)
        {
          // analytical derivative
          libMesh::Real f_analytic = derivs(dof);

          libMesh::Number soln = (*solution)(dof);

          // ref minus delta
          solution->set(dof,soln-delta);
          solution->close();
          system->assemble_qoi(qs);
          libMesh::Real qoi_m1 = _sim->get_qoi_value(0);

          // ref plus delta
          solution->set(dof,soln+delta);
          solution->close();
          system->assemble_qoi(qs);
          libMesh::Real qoi_p1 = _sim->get_qoi_value(0);

          // return to ref
          solution->set(dof,soln);
          solution->close();

          // central difference
          libMesh::Real fd_approx = (qoi_p1 - qoi_m1)/(2.0*delta);

          CPPUNIT_ASSERT_DOUBLES_EQUAL(fd_approx,f_analytic,libMesh::TOLERANCE);
        }
    }

    std::string spectroscopic_string(const std::string & qoi_name, const std::string & qoi_string, unsigned int nx, unsigned int ny)
    {
      std::string text = "[Mesh]\n";
                  text +=  "[./Generation]\n";
                  text +=    "n_elems_x = '"+std::to_string(nx)+"'\n";
                  text +=    "n_elems_y = '"+std::to_string(ny)+"'\n";
                  text += "[]\n";
                  text += "[QoI]\n";
                  text +=   "enabled_qois = '"+qoi_name+"'\n";
                  text +=   "[./"+qoi_string+"]\n";
                  text +=     "material = 'TestMaterial'\n";
                  text +=     "species_of_interest = 'CO2'\n";
                  text +=     "hitran_data_file = './test_data/CO2_data.dat'\n";
                  text +=     "hitran_partition_function_file = './test_data/CO2_partition_function.dat'\n";
                  text +=     "partition_temperatures = '290 310 0.01'\n";
                  text +=     "desired_wavenumber = '3682.7649'\n";
                  text +=     "min_wavenumber = '3682.69'\n";
                  text +=     "max_wavenumber = '3682.8'\n";
                  text +=     "calc_thermo_pressure = 'false'\n";
                  text +=       "[./Rayfire]\n";
                  text +=         "origin = '0.0 0.025'\n";
                  text +=         "theta = '0.0'\n";
                  text += "[]";

      return text;
    }

    std::string laser_string_2D(const std::string & qoi_name, const std::string & qoi_string,
                             const std::string & top_origin, const std::string & centerline_origin, const std::string & bottom_origin,
                             libMesh::Real theta, unsigned int nx, unsigned int ny)
    {
      std::string text = "[Mesh]\n";
                  text +=  "[./Generation]\n";
                  text +=    "n_elems_x = '"+std::to_string(nx)+"'\n";
                  text +=    "n_elems_y = '"+std::to_string(ny)+"'\n";
                  text += "[]\n";
                  text += "[QoI]\n";
                  text +=   "enabled_qois = '"+qoi_name+"'\n";
                  text +=   "[./"+qoi_string+"]\n";
                  text +=     "material = 'TestMaterial'\n";
                  text +=     "species_of_interest = 'CO2'\n";
                  text +=     "hitran_data_file = './test_data/CO2_data.dat'\n";
                  text +=     "hitran_partition_function_file = './test_data/CO2_partition_function.dat'\n";
                  text +=     "partition_temperatures = '290 310 0.01'\n";
                  text +=     "desired_wavenumber = '3682.7649'\n";
                  text +=     "min_wavenumber = '3682.69'\n";
                  text +=     "max_wavenumber = '3682.8'\n";
                  text +=     "calc_thermo_pressure = 'false'\n";
                  text +=     "n_quadrature_points = '4'\n";
                  text +=     "intensity_profile = 'collimated gaussian'\n";
                  text +=     "w = '0.000848725'\n";
                  text +=     "top_origin = '"+top_origin+"'\n";
                  text +=     "centerline_origin = '"+centerline_origin+"'\n";
                  text +=     "bottom_origin = '"+bottom_origin+"'\n";
                  text +=       "[./Rayfire]\n";
                  text +=         "theta = '"+std::to_string(theta)+"'\n";
                  text += "[]";

      return text;
    }

  private:
    void T_param_derivatives(std::shared_ptr<AbsorptionCoeffTesting<GRINS::AntiochChemistry> >absorb, libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> & Y, unsigned int i)
    {
      libMesh::Real delta = 1.0e-6;

      std::vector<libMesh::Real> T_analytic;
      T_analytic.push_back(absorb->dS_dT(T,P,i));
      T_analytic.push_back(absorb->d_nuC_dT(T,P,Y,i));
      T_analytic.push_back(absorb->d_nuD_dT(T,P,i));
      T_analytic.push_back(absorb->d_voigt_a_dT(T,P,Y,i));
      T_analytic.push_back(absorb->d_voigt_w_dT(T,P,i));
      T_analytic.push_back(absorb->d_voigt_dT(T,P,Y,i));
      T_analytic.push_back(absorb->d_kv_dT(T,P,Y,i));

      std::vector<libMesh::Real> fd_plus;
      fd_plus.push_back(absorb->Sw(T+delta,P,i));
      fd_plus.push_back(absorb->nu_C(T+delta,P,Y,i));
      fd_plus.push_back(absorb->nu_D(T+delta,P,i));
      fd_plus.push_back(absorb->voigt_a(T+delta,P,Y,i));
      fd_plus.push_back(absorb->voigt_w(T+delta,P,i));
      fd_plus.push_back(absorb->voigt(T+delta,P,Y,i));
      fd_plus.push_back(absorb->kv(T+delta,P,Y,i));

      std::vector<libMesh::Real> fd_minus;
      fd_minus.push_back(absorb->Sw(T-delta,P,i));
      fd_minus.push_back(absorb->nu_C(T-delta,P,Y,i));
      fd_minus.push_back(absorb->nu_D(T-delta,P,i));
      fd_minus.push_back(absorb->voigt_a(T-delta,P,Y,i));
      fd_minus.push_back(absorb->voigt_w(T-delta,P,i));
      fd_minus.push_back(absorb->voigt(T-delta,P,Y,i));
      fd_minus.push_back(absorb->kv(T-delta,P,Y,i));

      check_param_derivatives(T_analytic,fd_plus,fd_minus,delta);
    }

    void P_param_derivatives(std::shared_ptr<AbsorptionCoeffTesting<GRINS::AntiochChemistry> >absorb, libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> & Y, unsigned int i)
    {
      libMesh::Real delta = 1.0e-3;

      std::vector<libMesh::Real> P_analytic;
      P_analytic.push_back(absorb->dS_dP(T,P,i));
      P_analytic.push_back(absorb->d_nuC_dP(T,Y,i));
      P_analytic.push_back(absorb->d_nuD_dP(T,i));
      P_analytic.push_back(absorb->d_voigt_a_dP(T,P,Y,i));
      P_analytic.push_back(absorb->d_voigt_w_dP(T,P,i));
      P_analytic.push_back(absorb->d_voigt_dP(T,P,Y,i));
      P_analytic.push_back(absorb->d_kv_dP(T,P,Y,i));

      std::vector<libMesh::Real> fd_plus;
      fd_plus.push_back(absorb->Sw(T,P+delta,i));
      fd_plus.push_back(absorb->nu_C(T,P+delta,Y,i));
      fd_plus.push_back(absorb->nu_D(T,P+delta,i));
      fd_plus.push_back(absorb->voigt_a(T,P+delta,Y,i));
      fd_plus.push_back(absorb->voigt_w(T,P+delta,i));
      fd_plus.push_back(absorb->voigt(T,P+delta,Y,i));
      fd_plus.push_back(absorb->kv(T,P+delta,Y,i));

      std::vector<libMesh::Real> fd_minus;
      fd_minus.push_back(absorb->Sw(T,P-delta,i));
      fd_minus.push_back(absorb->nu_C(T,P-delta,Y,i));
      fd_minus.push_back(absorb->nu_D(T,P-delta,i));
      fd_minus.push_back(absorb->voigt_a(T,P-delta,Y,i));
      fd_minus.push_back(absorb->voigt_w(T,P-delta,i));
      fd_minus.push_back(absorb->voigt(T,P-delta,Y,i));
      fd_minus.push_back(absorb->kv(T,P-delta,Y,i));

      check_param_derivatives(P_analytic,fd_plus,fd_minus,delta);
    }

    void Y_param_derivatives(std::shared_ptr<AbsorptionCoeffTesting<GRINS::AntiochChemistry> >absorb, libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> & Y, unsigned int i)
    {
      libMesh::Real delta = 1.0e-8;

      libMesh::Real y_base = Y[0];
      std::vector<libMesh::Real> Y_analytic;
      Y_analytic.push_back(absorb->dX_dY(Y,0));
      Y_analytic.push_back(absorb->d_nuC_dY(T,P,Y,0,i));
      Y_analytic.push_back(absorb->d_voigt_a_dY(T,P,Y,0,i));
      Y_analytic.push_back(absorb->d_voigt_dY(T,P,Y,0,i));
      Y_analytic.push_back(absorb->d_kv_dY(T,P,Y,0,i));

      std::vector<libMesh::Real> fd_plus;
      Y[0] = y_base + delta;
      fd_plus.push_back(absorb->_chemistry->X(absorb->_species_idx,absorb->_chemistry->M_mix(Y),Y[absorb->_species_idx]));
      fd_plus.push_back(absorb->nu_C(T,P,Y,i));
      fd_plus.push_back(absorb->voigt_a(T,P,Y,i));
      fd_plus.push_back(absorb->voigt(T,P,Y,i));
      fd_plus.push_back(absorb->kv(T,P,Y,i));

      std::vector<libMesh::Real> fd_minus;
      Y[0] = y_base - delta;
      fd_minus.push_back(absorb->_chemistry->X(absorb->_species_idx,absorb->_chemistry->M_mix(Y),Y[absorb->_species_idx]));
      fd_minus.push_back(absorb->nu_C(T,P,Y,i));
      fd_minus.push_back(absorb->voigt_a(T,P,Y,i));
      fd_minus.push_back(absorb->voigt(T,P,Y,i));
      fd_minus.push_back(absorb->kv(T,P,Y,i));

      Y[0] = y_base;

      check_param_derivatives(Y_analytic,fd_plus,fd_minus,delta);
    }

    void check_param_derivatives(std::vector<libMesh::Real> & analytic, std::vector<libMesh::Real> & fd_plus, std::vector<libMesh::Real> & fd_minus, libMesh::Real delta)
    {
      for (unsigned int d=0; d<analytic.size(); ++d)
        {
          libMesh::Real f_analytic = analytic[d];
          libMesh::Real fd_approx = (fd_plus[d] - fd_minus[d])/(2.0*delta);

          CPPUNIT_ASSERT_DOUBLES_EQUAL(fd_approx,f_analytic,libMesh::TOLERANCE);
        }
    }

    void T_elem_derivative(GRINS::MultiphysicsSystem * system, GRINS::AssemblyContext * context)
    {
      GRINS::PrimitiveTempFEVariables T_var( GRINS::GRINSPrivate::VariableWarehouse::get_variable_subclass<GRINS::PrimitiveTempFEVariables>("Temperature") );
      libMesh::Real delta = 1e-6;
      test_var_elem_derivs(system,context,T_var.T(),delta);
    }

    void P_elem_derivative(GRINS::MultiphysicsSystem * system, GRINS::AssemblyContext * context)
    {
      GRINS::PressureFEVariable P_var( GRINS::GRINSPrivate::VariableWarehouse::get_variable_subclass<GRINS::PressureFEVariable>("Pressure") );
      libMesh::Real delta = 1e-3;
      test_var_elem_derivs(system,context,P_var.p(),delta);
    }

    void Y_elem_derivative(GRINS::MultiphysicsSystem * system, GRINS::AssemblyContext * context)
    {
      GRINS::SpeciesMassFractionsVariable Y_var( GRINS::GRINSPrivate::VariableWarehouse::get_variable_subclass<GRINS::SpeciesMassFractionsVariable>("SpeciesMassFractions") );
      libMesh::Real delta = 1e-8;
      test_var_elem_derivs(system,context,Y_var.species(0),delta);
    }

    void test_var_elem_derivs(GRINS::MultiphysicsSystem * system, GRINS::AssemblyContext * context, unsigned int var_index, libMesh::Real delta)
    {
      libMesh::DifferentiableQoI * qoi = system->get_qoi();
      libMesh::Number & qoi_value = context->get_qois()[0];
      qoi_value = 0.0;

      libMesh::QoISet qs;
      qs.add_index(0);

      libMesh::DenseSubVector<libMesh::Number> deriv = context->get_qoi_derivatives(0,var_index);
      libMesh::DenseSubVector<libMesh::Number> & solution  = context->get_elem_solution(var_index);

      for (unsigned int d=0; d<solution.size(); ++d)
        {
          libMesh::Number soln = solution.el(d);

          solution(d) = soln+delta;
          qoi->element_qoi(*context,qs);
          libMesh::Number qoi_p1 = qoi_value;
          qoi_value = 0.0;

          solution(d) = soln-delta;
          qoi->element_qoi(*context,qs);
          libMesh::Number qoi_m1= qoi_value;
          qoi_value = 0.0;

          solution(d) = soln;

          libMesh::Real fd_approx = (qoi_p1 - qoi_m1)/(2.0*delta);

          CPPUNIT_ASSERT_DOUBLES_EQUAL(fd_approx,deriv(d),libMesh::TOLERANCE);
        }
    }

  };

} // namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT

#endif //GRINS_SPECTROSCOPIC_TEST_BASE_H
