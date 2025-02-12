/**
 * @file IDynamoBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Model/Boussinesq/Shell/Dynamo/IDynamoBackend.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/Insulating.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/Stencil.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperatorBoundary/FieldToRhs.hpp"
#include "QuICC/ModelOperatorBoundary/SolverHasBc.hpp"
#include "QuICC/ModelOperatorBoundary/SolverNoTau.hpp"
#include "QuICC/ModelOperatorBoundary/Stencil.hpp"
#include "QuICC/NonDimensional/CflAlfvenDamping.hpp"
#include "QuICC/NonDimensional/CflAlfvenScale.hpp"
#include "QuICC/NonDimensional/CflInertial.hpp"
#include "QuICC/NonDimensional/CflTorsional.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/Heating.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/MagneticPrandtl.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/RRatio.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/D1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/D2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/ICondition.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/InsulatingShell.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/Operator.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/R1D1DivR1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/Value.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y2SphLapl.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y3.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y3SphLapl.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y4.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y4SphLapl.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y4SphLapl2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Id.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/D1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/InsulatingShell.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/R1D1DivR1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/Value.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/ValueD1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/ValueD2.hpp"
#include "QuICC/Tools/IdToHuman.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace Dynamo {

std::vector<std::string> IDynamoBackend::fieldNames() const
{
   std::vector<std::string> names = {PhysicalNames::Velocity().tag(),
      PhysicalNames::Magnetic().tag(), PhysicalNames::Temperature().tag()};

   return names;
}

std::vector<std::string> IDynamoBackend::paramNames() const
{
   std::vector<std::string> names = {NonDimensional::Prandtl().tag(),
      NonDimensional::MagneticPrandtl().tag(), NonDimensional::Rayleigh().tag(),
      NonDimensional::Ekman().tag(), NonDimensional::Heating().tag(),
      NonDimensional::RRatio().tag()};

   return names;
}

std::vector<bool> IDynamoBackend::isPeriodicBox() const
{
   std::vector<bool> periodic = {false, false, false};

   return periodic;
}

std::map<std::string, MHDFloat> IDynamoBackend::automaticParameters(
   const std::map<std::string, MHDFloat>& cfg) const
{
   auto E = cfg.find(NonDimensional::Ekman().tag())->second;
   auto Pm = cfg.find(NonDimensional::MagneticPrandtl().tag())->second;
   auto rratio = cfg.find(NonDimensional::RRatio().tag())->second;

   std::map<std::string, MHDFloat> params = {
      {NonDimensional::CflInertial().tag(), 0.1 * E / Pm},
      {NonDimensional::CflTorsional().tag(), 0.1 * std::sqrt(E)},
      {NonDimensional::CflAlfvenScale().tag(), Pm / E},
      {NonDimensional::CflAlfvenDamping().tag(), (1.0 + Pm) / 2.0},
   };

   bool useGapWidth = true;
   if (useGapWidth)
   {
      params.emplace(NonDimensional::Lower1d().tag(), rratio / (1.0 - rratio));
      params.emplace(NonDimensional::Upper1d().tag(), 1.0 / (1.0 - rratio));
   }
   else
   {
      params.emplace(NonDimensional::Lower1d().tag(), rratio);
      params.emplace(NonDimensional::Upper1d().tag(), 1.0);
   }

   return params;
}

int IDynamoBackend::nBc(const SpectralFieldId& fId) const
{
   int nBc = 0;

   if (fId == std::make_pair(PhysicalNames::Velocity::id(),
                 FieldComponents::Spectral::TOR) ||
       fId == std::make_pair(PhysicalNames::Temperature::id(),
                 FieldComponents::Spectral::SCALAR))
   {
      nBc = 2;
   }
   else if (fId == std::make_pair(PhysicalNames::Velocity::id(),
                      FieldComponents::Spectral::POL))
   {
      nBc = 4;
   }
   else
   {
      nBc = 0;
   }

   return nBc;
}

void IDynamoBackend::applyTau(SparseMatrix& mat, const SpectralFieldId& rowId,
   const SpectralFieldId& colId, const int l,
   std::shared_ptr<details::BlockOptions> opts, const Resolution& res,
   const BcMap& bcs, const NonDimensional::NdMap& nds,
   const bool isSplitOperator) const
{
   auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);

   auto bcId = bcs.find(rowId.first)->second;

   auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
   auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();

   typedef SparseSM::Chebyshev::LinearMap::Boundary::ICondition::Position
      Position;

   SparseSM::Chebyshev::LinearMap::Boundary::Operator bcOp(nN, nN, ri, ro);

   if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::TOR) &&
       rowId == colId)
   {
      if (bcId == Bc::Name::NoSlip::id())
      {
         bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(
            Position::TOP);
         bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(
            Position::BOTTOM);
      }
      else if (bcId == Bc::Name::StressFree::id())
      {
         bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::R1D1DivR1>(
            Position::TOP);
         bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::R1D1DivR1>(
            Position::BOTTOM);
      }
      else
      {
         throw std::logic_error("Boundary conditions for Velocity Toroidal "
                                "component not implemented");
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
                        FieldComponents::Spectral::POL) &&
            rowId == colId)
   {
      if (this->useSplitEquation())
      {
         if (isSplitOperator)
         {
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(
               Position::TOP);
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(
               Position::BOTTOM);
         }
         else
         {
            if (bcId == Bc::Name::NoSlip::id())
            {
               bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::D1>(
                  Position::TOP);
               bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::D1>(
                  Position::BOTTOM);
            }
            else if (bcId == Bc::Name::StressFree::id())
            {
               bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::D2>(
                  Position::TOP);
               bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::D2>(
                  Position::BOTTOM);
            }
            else
            {
               throw std::logic_error("Boundary conditions for Velocity "
                                      "Poloidal component not implemented");
            }
         }
      }
      else
      {
         if (bcId == Bc::Name::NoSlip::id())
         {
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(
               Position::TOP);
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(
               Position::BOTTOM);
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::D1>(
               Position::TOP);
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::D1>(
               Position::BOTTOM);
         }
         else if (bcId == Bc::Name::StressFree::id())
         {
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(
               Position::TOP);
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(
               Position::BOTTOM);
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::D2>(
               Position::TOP);
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::D2>(
               Position::BOTTOM);
         }
         else
         {
            throw std::logic_error("Boundary conditions for Velocity Poloidal "
                                   "component not implemented");
         }
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Magnetic::id(),
                        FieldComponents::Spectral::TOR) &&
            rowId == colId)
   {
      if (bcId == Bc::Name::Insulating::id())
      {
         bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(
            Position::TOP);
         bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(
            Position::BOTTOM);
      }
      else
      {
         throw std::logic_error("Boundary conditions for Magnetic Toroidal "
                                "component not implemented");
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Magnetic::id(),
                        FieldComponents::Spectral::POL) &&
            rowId == colId)
   {
      if (bcId == Bc::Name::Insulating::id())
      {
         bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::InsulatingShell>(
            Position::TOP, l);
         bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::InsulatingShell>(
            Position::BOTTOM, l);
      }
      else
      {
         throw std::logic_error("Boundary conditions for Magnetic Poloidal "
                                "component not implemented");
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Temperature::id(),
                        FieldComponents::Spectral::SCALAR) &&
            rowId == colId)
   {
      if (bcId == Bc::Name::FixedTemperature::id())
      {
         bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(
            Position::TOP);
         bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(
            Position::BOTTOM);
      }
      else if (bcId == Bc::Name::FixedFlux::id())
      {
         bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::D1>(
            Position::TOP);
         bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::D1>(
            Position::BOTTOM);
      }
      else
      {
         throw std::logic_error(
            "Boundary conditions for Temperature not implemented (" +
            std::to_string(bcId) + ")");
      }
   }

   mat += bcOp.mat();
}

void IDynamoBackend::stencil(SparseMatrix& mat, const SpectralFieldId& fieldId,
   const int l, const Resolution& res, const bool makeSquare, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);

   auto bcId = bcs.find(fieldId.first)->second;

   auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
   auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();

   int s = this->nBc(fieldId);
   if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
                     FieldComponents::Spectral::TOR))
   {
      if (bcId == Bc::Name::NoSlip::id())
      {
         SparseSM::Chebyshev::LinearMap::Stencil::Value bc(nN, nN - s, ri, ro);
         mat = bc.mat();
      }
      else if (bcId == Bc::Name::StressFree::id())
      {
         SparseSM::Chebyshev::LinearMap::Stencil::R1D1DivR1 bc(nN, nN - s, ri,
            ro);
         mat = bc.mat();
      }
      else
      {
         throw std::logic_error("Galerkin boundary conditions for Velocity "
                                "Toroidal component not implemented");
      }
   }
   else if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
                          FieldComponents::Spectral::POL))
   {
      if (bcId == Bc::Name::NoSlip::id())
      {
         SparseSM::Chebyshev::LinearMap::Stencil::ValueD1 bc(nN, nN - s, ri,
            ro);
         mat = bc.mat();
      }
      else if (bcId == Bc::Name::StressFree::id())
      {
         SparseSM::Chebyshev::LinearMap::Stencil::ValueD2 bc(nN, nN - s, ri,
            ro);
         mat = bc.mat();
      }
      else
      {
         throw std::logic_error("Galerin boundary conditions for Velocity "
                                "Poloidal component not implemented");
      }
   }
   else if (fieldId == std::make_pair(PhysicalNames::Magnetic::id(),
                          FieldComponents::Spectral::TOR))
   {
      if (bcId == Bc::Name::Insulating::id())
      {
         SparseSM::Chebyshev::LinearMap::Stencil::Value bc(nN, nN - s, ri, ro);
         mat = bc.mat();
      }
      else
      {
         throw std::logic_error("Galerkin boundary conditions for Magnetic "
                                "Toroidal component not implemented");
      }
   }
   else if (fieldId == std::make_pair(PhysicalNames::Magnetic::id(),
                          FieldComponents::Spectral::POL))
   {
      if (bcId == Bc::Name::Insulating::id())
      {
         SparseSM::Chebyshev::LinearMap::Stencil::InsulatingShell bc(nN, nN - s,
            ri, ro, l);
         mat = bc.mat();
      }
      else
      {
         throw std::logic_error("Galerkin boundary conditions for Magnetic "
                                "Poloidal component not implemented");
      }
   }
   else if (fieldId == std::make_pair(PhysicalNames::Temperature::id(),
                          FieldComponents::Spectral::SCALAR))
   {
      if (bcId == Bc::Name::FixedTemperature::id())
      {
         SparseSM::Chebyshev::LinearMap::Stencil::Value bc(nN, nN - s, ri, ro);
         mat = bc.mat();
      }
      else if (bcId == Bc::Name::FixedFlux::id())
      {
         SparseSM::Chebyshev::LinearMap::Stencil::D1 bc(nN, nN - s, ri, ro);
         mat = bc.mat();
      }
      else
      {
         throw std::logic_error(
            "Galerkin boundary conditions for Temperature not implemented");
      }
   }

   if (makeSquare)
   {
      SparseSM::Chebyshev::LinearMap::Id qId(nN - s, nN, ri, ro);
      mat = qId.mat() * mat;
   }
}

void IDynamoBackend::applyGalerkinStencil(SparseMatrix& mat,
   const SpectralFieldId& rowId, const SpectralFieldId& colId, const int lr,
   const int lc, std::shared_ptr<details::BlockOptions> opts,
   const Resolution& res, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, lr)(0);

   auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
   auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();

   auto S = mat;
   this->stencil(S, colId, lc, res, false, bcs, nds);

   auto s = this->nBc(rowId);
   SparseSM::Chebyshev::LinearMap::Id qId(nN - s, nN, ri, ro, 0, s);
   mat = qId.mat() * (mat * S);
}

namespace implDetails {

MHDFloat effectiveRa(const NonDimensional::NdMap& nds)
{
   auto Pm = nds.find(NonDimensional::MagneticPrandtl::id())->second->value();
   auto effRa = nds.find(NonDimensional::Rayleigh::id())->second->value();
   auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();

   // Scaled on gap width
   if (ro != 1.0)
   {
      auto T = 1.0 / nds.find(NonDimensional::Ekman::id())->second->value();
      effRa *= Pm * Pm * T / ro;
   }

   return effRa;
}

MHDFloat effectiveBg(const NonDimensional::NdMap& nds)
{
   MHDFloat effBg = 1.0;
   auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();
   auto rratio = nds.find(NonDimensional::RRatio::id())->second->value();
   auto heatingMode = nds.find(NonDimensional::Heating::id())->second->value();

   if (ro == 1.0)
   {
      // Nothing
      //
   }
   // gap width and internal heating
   else if (heatingMode == 0)
   {
      effBg = 2.0 / (ro * (1.0 + rratio));
   }
   // gap width and differential heating
   else if (heatingMode == 1)
   {
      effBg = ro * ro * rratio;
   }

   return effBg;
}

} // namespace implDetails
} // namespace Dynamo
} // namespace Shell
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
