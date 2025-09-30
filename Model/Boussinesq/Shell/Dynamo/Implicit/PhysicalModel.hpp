/**
 * @file PhysicalModel.hpp
 * @brief Implementation of the Boussinesq thermal convection dynamo in a
 * spherical shell model (Toroidal/Poloidal formulation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SHELL_DYNAMO_IMPLICIT_PHYSICALMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_SHELL_DYNAMO_IMPLICIT_PHYSICALMODEL_HPP

// System includes
//
#include <string>

// Project includes
//
#include "QuICC/SpatialScheme/3D/SLFm.hpp"
#include "QuICC/Model/PyModelBackend.hpp"
#include "Model/Boussinesq/Shell/Dynamo/Implicit/ModelBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace Dynamo {

namespace Implicit {

/**
 * @brief Implementation of the Boussinesq thermal convection dynamo in a
 * spherical shell model (Toroidal/Poloidal formulation)
 */
template <typename TBuilder> class PhysicalModel : public TBuilder
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

   /**
    * @brief Initialize specialized backend
    */
   void init() final;

protected:
private:
};

template <typename TBuilder> void PhysicalModel<TBuilder>::init()
{
   TBuilder::init();
#ifdef QUICC_MODEL_BOUSSINESQSHELLDYNAMO_IMPLICIT_BACKEND_CPP

   this->mpBackend = std::make_shared<ModelBackend>();
#else
   std::string pyModule = "boussinesq.shell.dynamo.implicit.physical_model";
   std::string pyClass = "PhysicalModel";

   this->mpBackend =
      std::make_shared<PyModelBackend>(pyModule, pyClass);
#endif
}

} // namespace Implicit
} // namespace Dynamo
} // namespace Shell
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SHELL_DYNAMO_IMPLICIT_PHYSICALMODEL_HPP
