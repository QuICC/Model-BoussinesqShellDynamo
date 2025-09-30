/**
 * @file IDynamoModel.cpp
 * @brief Source of the Boussinesq thermal convection dynamo in a spherical
 * shell (Toroidal/Poloidal formulation)
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Shell/Dynamo/IDynamoModel.hpp"
#include "Model/Boussinesq/Shell/Dynamo/Induction.hpp"
#include "Model/Boussinesq/Shell/Dynamo/Momentum.hpp"
#include "Model/Boussinesq/Shell/Dynamo/Transport.hpp"
#include "Model/Boussinesq/Shell/Dynamo/gitHash.hpp"
#include "QuICC/Io/Variable/ShellNusseltWriter.hpp"
#include "QuICC/Io/Variable/ShellScalarEnergyWriter.hpp"
#include "QuICC/Io/Variable/ShellScalarLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/ShellScalarMSpectrumWriter.hpp"
#include "QuICC/Io/Variable/ShellTorPolEnergyWriter.hpp"
#include "QuICC/Io/Variable/ShellTorPolLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/ShellTorPolMSpectrumWriter.hpp"
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace Dynamo {

VectorFormulation::Id IDynamoModel::SchemeFormulation()
{
   return VectorFormulation::TORPOL;
}

std::string IDynamoModel::version() const
{
   return std::string(gitHash);
}

void IDynamoModel::addEquations(SharedSimulation spSim)
{
   // Add transport equation
   spSim->addEquation<Equations::Boussinesq::Shell::Dynamo::Transport>(
      this->spBackend());

   // Add Navier-Stokes equation
   spSim->addEquation<Equations::Boussinesq::Shell::Dynamo::Momentum>(
      this->spBackend());

   // Add induction equation
   spSim->addEquation<Equations::Boussinesq::Shell::Dynamo::Induction>(
      this->spBackend());
}

std::map<std::string, std::map<std::string, int>>
IDynamoModel::configTags() const
{
   std::map<std::string, int> onOff;
   onOff.emplace("enable", 1);

   std::map<std::string, int> options;
   options.emplace("enable", 0);
   options.emplace("numbered", 0);
   options.emplace("only_every", 1);

   std::map<std::string, std::map<std::string, int>> tags;
   // kinetic
   tags.emplace("kinetic_energy", onOff);
   tags.emplace("kinetic_l_spectrum", options);
   tags.emplace("kinetic_m_spectrum", options);
   // magnetic
   tags.emplace("magnetic_energy", onOff);
   tags.emplace("magnetic_l_spectrum", options);
   tags.emplace("magnetic_m_spectrum", options);
   // temperature
   tags.emplace("temperature_energy", onOff);
   tags.emplace("temperature_l_spectrum", options);
   tags.emplace("temperature_m_spectrum", options);
   tags.emplace("temperature_nusselt", onOff);

   return tags;
}

void IDynamoModel::addAsciiOutputFiles(SharedSimulation spSim)
{
   // Create temperature energy writer
   this->enableAsciiFile<Io::Variable::ShellScalarEnergyWriter>(
      "temperature_energy", "temperature", PhysicalNames::Temperature::id(),
      spSim);

   // Create temperature L energy spectrum writer
   this->enableAsciiFile<Io::Variable::ShellScalarLSpectrumWriter>(
      "temperature_l_spectrum", "temperature", PhysicalNames::Temperature::id(),
      spSim);

   // Create temperature M energy spectrum writer
   this->enableAsciiFile<Io::Variable::ShellScalarMSpectrumWriter>(
      "temperature_m_spectrum", "temperature", PhysicalNames::Temperature::id(),
      spSim);

   // Create kinetic energy writer
   this->enableAsciiFile<Io::Variable::ShellTorPolEnergyWriter>(
      "kinetic_energy", "kinetic", PhysicalNames::Velocity::id(), spSim);

   // Create kinetic L energy spectrum writer
   this->enableAsciiFile<Io::Variable::ShellTorPolLSpectrumWriter>(
      "kinetic_l_spectrum", "kinetic", PhysicalNames::Velocity::id(), spSim);

   // Create kinetic M energy spectrum writer
   this->enableAsciiFile<Io::Variable::ShellTorPolMSpectrumWriter>(
      "kinetic_m_spectrum", "kinetic", PhysicalNames::Velocity::id(), spSim);

   // Create magnetic energy writer
   this->enableAsciiFile<Io::Variable::ShellTorPolEnergyWriter>(
      "magnetic_energy", "magnetic", PhysicalNames::Magnetic::id(), spSim);

   // Create magnetic L energy spectrum writer
   this->enableAsciiFile<Io::Variable::ShellTorPolLSpectrumWriter>(
      "magnetic_l_spectrum", "magnetic", PhysicalNames::Magnetic::id(), spSim);

   // Create magnetic M energy spectrum writer
   this->enableAsciiFile<Io::Variable::ShellTorPolMSpectrumWriter>(
      "magnetic_m_spectrum", "magnetic", PhysicalNames::Magnetic::id(), spSim);

   // Create nusselt number writer
   this->enableAsciiFile<Io::Variable::ShellNusseltWriter>(
      "temperature_nusselt", "temperature_", PhysicalNames::Temperature::id(),
      spSim);
}

} // namespace Dynamo
} // namespace Shell
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
