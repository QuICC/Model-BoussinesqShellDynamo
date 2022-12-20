/**
 * @file PhysicalModel.hpp
 * @brief Implementation of the Boussinesq thermal convection dynamo in a spherical shell model (Toroidal/Poloidal formulation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SHELL_DYNAMO_IMPLICIT_PHYSICALMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_SHELL_DYNAMO_IMPLICIT_PHYSICALMODEL_HPP

// System includes
//
#include <string>

// Project includes
//
#include "QuICC/Model/Boussinesq/Shell/Dynamo/IDynamoModel.hpp"
#include "QuICC/SpatialScheme/3D/SLFm.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace Dynamo {

namespace Implicit {

   /**
    * @brief Implementation of the Boussinesq thermal convection dynamo in a spherical shell model (Toroidal/Poloidal formulation)
    */
   class PhysicalModel: public IDynamoModel
   {
      public:
         /// Typedef for the spatial scheme used
         typedef SpatialScheme::SLFm SchemeType;

         /**
          * @brief Constructor
          */
         PhysicalModel() = default;

         /**
          * @brief Destructor
          */
         virtual ~PhysicalModel() = default;

         /// Python script/module name
         virtual std::string PYMODULE() override;

      protected:

      private:
   };

}
}
}
}
}
}

#endif // QUICC_MODEL_BOUSSINESQ_SHELL_DYNAMO_IMPLICIT_PHYSICALMODEL_HPP
