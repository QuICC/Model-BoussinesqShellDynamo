/**
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq thermal convection dynamo in a spherical shell (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq//Shell/Dynamo/Explicit/PhysicalModel.hpp"
#include "Model/Boussinesq//Shell/Dynamo/Explicit/ModelBackend.hpp"
#include "QuICC/Model/PyModelBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace Dynamo {

namespace Explicit {

   std::string PhysicalModel::PYMODULE()
   {
      return "boussinesq.shell.dynamo.explicit.physical_model";
   }

   void PhysicalModel::init()
   {
#ifdef QUICC_MODEL_BOUSSINESQSHELLDYNAMO_EXPLICIT_BACKEND_CPP
      IPhysicalModel<Simulation,StateGenerator,VisualizationGenerator>::init();

      this->mpBackend = std::make_shared<ModelBackend>();
#else
      IPhysicalPyModel<Simulation,StateGenerator,VisualizationGenerator>::init();

      this->mpBackend = std::make_shared<PyModelBackend>(this->PYMODULE(), this->PYCLASS());
#endif
   }

} // Explicit
} // Dynamo
} // Shell
} // Boussinesq
} // Model
} // QuICC
