//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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

#include "grins_config.h"

#ifdef GRINS_HAVE_CPPUNIT

#ifdef GRINS_HAVE_ANTIOCH

#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

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
#include "libmesh/parsed_function.h"

// Ignore warnings from auto_ptr in CPPUNIT_TEST_SUITE_END()
#include <libmesh/ignore_warnings.h>

namespace GRINSTesting
{
  class SpectroscopicAbsorptionTest : public CppUnit::TestCase
  {
  public:
    CPPUNIT_TEST_SUITE( SpectroscopicAbsorptionTest );
    CPPUNIT_TEST( single_elem_mesh );
    CPPUNIT_TEST( multi_elem_mesh );
    CPPUNIT_TEST( param_derivs );
    CPPUNIT_TEST( elem_qoi_derivatives );
    CPPUNIT_TEST_SUITE_END();

  public:

    void tearDown()
    {
      // Clear out the VariableWarehouse so it doesn't interfere with other tests.
      GRINS::GRINSPrivate::VariableWarehouse::clear();
    }

    //! Single QUAD4 elem, uniform T,P,Y
    void single_elem_mesh()
    {
      const std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/spectroscopic_absorption_qoi.in";
      libMesh::Real calc_answer = 0.520403290868787;

      this->run_test(filename,calc_answer);
    }

    //! 10x10 mesh, uniform T,P,Y
    void multi_elem_mesh()
    {
      const std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/spectroscopic_absorption_qoi_fine.in";
      libMesh::Real calc_answer = 0.520403290868787; // same physical conditions as single_elem_mesh, so answer should not change

      this->run_test(filename,calc_answer);
    }

    void param_derivs()
    {
      const std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/spectroscopic_absorption_qoi.in";
      this->init_sim(filename);

      std::string material = "TestMaterial";
      std::string hitran_data = "./test_data/CO2_data.dat";
      std::string hitran_partition = "./test_data/CO2_partition_function.dat";
      libMesh::Real T_min = 290.0,
                    T_max = 310.0,
                    T_step = 0.01;
      GRINS::SharedPtr<GRINS::HITRAN> hitran( new GRINS::HITRAN(hitran_data,hitran_partition,T_min,T_max,T_step) );

      std::string species = "CO2";
      libMesh::Real thermo_pressure = 5066.25;
      libMesh::Real nu_desired = 3682.7649;
      libMesh::Real nu_min = 3682.69;
      libMesh::Real nu_max = 3682.8;
      GRINS::ChemistryBuilder chem_builder;
      libMesh::UniquePtr<GRINS::AntiochChemistry> chem_ptr;
      chem_builder.build_chemistry(*(_input.get()),material,chem_ptr);
      GRINS::SharedPtr<GRINS::AntiochChemistry> chem(chem_ptr.release());
      GRINS::SharedPtr<AbsorptionCoeffTesting<GRINS::AntiochChemistry> > absorb = new AbsorptionCoeffTesting<GRINS::AntiochChemistry>(chem,hitran,nu_min,nu_max,nu_desired,species,thermo_pressure);

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

    void elem_qoi_derivatives()
    {
      const std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/spectroscopic_absorption_qoi.in";

      // run the Simulation once to ensure everything is initialized and populated
      this->init_sim(filename);
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

  private:
    GRINS::SharedPtr<GRINS::Simulation> _sim;
    GRINS::SharedPtr<GetPot> _input;

    //! Run the test on a given input file and calculated answer
    void run_test(const std::string filename, libMesh::Real calc_answer)
    {
      this->init_sim(filename);

      _sim->run();

      CPPUNIT_ASSERT_DOUBLES_EQUAL( calc_answer, _sim->get_qoi_value(0),libMesh::TOLERANCE  );
    }

    //! Initialize the GetPot and Simulation class objects
    void init_sim(const std::string & filename)
    {
      _input.reset(new GetPot(filename));

      const char * const argv = "unit_driver";
      GetPot empty_command_line( (const int)1,&argv );
      GRINS::SimulationBuilder sim_builder;

      _sim = new GRINS::Simulation(*_input,
                                   empty_command_line,
                                   sim_builder,
                                   *TestCommWorld );
    }

    void T_param_derivatives(GRINS::SharedPtr<AbsorptionCoeffTesting<GRINS::AntiochChemistry> >absorb, libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> & Y, unsigned int i)
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

    void P_param_derivatives(GRINS::SharedPtr<AbsorptionCoeffTesting<GRINS::AntiochChemistry> >absorb, libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> & Y, unsigned int i)
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

    void Y_param_derivatives(GRINS::SharedPtr<AbsorptionCoeffTesting<GRINS::AntiochChemistry> >absorb, libMesh::Real T, libMesh::Real P, std::vector<libMesh::Real> & Y, unsigned int i)
    {
      libMesh::Real delta = 1.0e-8;

      libMesh::Real y_base = Y[0];
      std::vector<libMesh::Real> Y_analytic;
      Y_analytic.push_back(absorb->dX_dY(Y));
      Y_analytic.push_back(absorb->d_nuC_dY(T,P,Y,i));
      Y_analytic.push_back(absorb->d_voigt_a_dY(T,P,Y,i));
      Y_analytic.push_back(absorb->d_voigt_dY(T,P,Y,i));
      Y_analytic.push_back(absorb->d_kv_dY(T,P,Y,i));

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

  CPPUNIT_TEST_SUITE_REGISTRATION( SpectroscopicAbsorptionTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_ANTIOCH
#endif // GRINS_HAVE_CPPUNIT
