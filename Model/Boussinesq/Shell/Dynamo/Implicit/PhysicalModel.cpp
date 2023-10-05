/**
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq thermal convection dynamo in a spherical shell (Toroidal/Poloidal formulation)
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq//Shell/Dynamo/Implicit/PhysicalModel.hpp"
#include "Model/Boussinesq//Shell/Dynamo/Implicit/ModelBackend.hpp"
#include "QuICC/Model/PyModelBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace Dynamo {

namespace Implicit {

   std::string PhysicalModel::PYMODULE()
   {
      return "boussinesq.shell.dynamo.implicit.physical_model";
   }

   void PhysicalModel::init()
   {
#ifdef QUICC_MODEL_BOUSSINESQSHELLDYNAMO_IMPLICIT_BACKEND_CPP
      IPhysicalModel<Simulation,StateGenerator,VisualizationGenerator>::init();

      this->mpBackend = std::make_shared<ModelBackend>();
#else
      IPhysicalPyModel<Simulation,StateGenerator,VisualizationGenerator>::init();

      this->mpBackend = std::make_shared<PyModelBackend>(this->PYMODULE(), this->PYCLASS());
#endif
   }

} // Implicit
} // Dynamo
} // Shell
} // Boussinesq
} // Model
} // QuICC
