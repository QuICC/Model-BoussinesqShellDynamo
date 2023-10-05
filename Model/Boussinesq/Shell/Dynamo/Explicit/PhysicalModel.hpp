/**
 * @file PhysicalModel.hpp
 * @brief Implementation of the Boussinesq thermal convection dynamo in a spherical shell model (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SHELL_DYNAMO_EXPLICIT_PHYSICALMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_SHELL_DYNAMO_EXPLICIT_PHYSICALMODEL_HPP

// System includes
//
#include <string>

// Project includes
//
#include "Model/Boussinesq//Shell/Dynamo/IDynamoModel.hpp"
#include "QuICC/SpatialScheme/3D/SLFl.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace Dynamo {

namespace Explicit {

   /**
    * @brief Implementation of the Boussinesq thermal convection dynamo in a spherical shell model (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
    */
   class PhysicalModel: public IDynamoModel
   {
      public:
         /// Typedef for the spatial scheme used
         typedef SpatialScheme::SLFl SchemeType;

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

         /**
          * @brief Initialize specialized backend
          */
         void init() final;

      protected:

      private:
   };

}
}
}
}
}
}

#endif // QUICC_MODEL_BOUSSINESQ_SHELL_DYNAMO_EXPLICIT_PHYSICALMODEL_HPP
