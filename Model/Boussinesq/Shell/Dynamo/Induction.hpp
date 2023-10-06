/**
 * @file Induction.hpp
 * @brief Implementation of the vector induction equation for the Boussinesq thermal convection dynamo spherical shell
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SHELL_DYNAMO_INDUCTION_HPP
#define QUICC_MODEL_BOUSSINESQ_SHELL_DYNAMO_INDUCTION_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Equations/IVectorEquation.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Shell {

namespace Dynamo {

   /**
    * @brief Implementation of the vector induction equation for the Boussinesq thermal convection dynamo in a spherical shell
    */
   class Induction: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         Induction(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~Induction();

         /**
          * @brief Initialize nonlinear interaction kernel
          */
         virtual void initNLKernel(const bool force = false) override;

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements() override;

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling() override;

         /**
          * @brief Set the nonlinear integration components
          */
         virtual void setNLComponents() override;

      private:
   };

}
}
}
}
}

#endif // QUICC_MODEL_BOUSSINESQ_SHELL_DYNAMO_INDUCTION_HPP
