/** 
 * @file Induction.cpp
 * @brief Source of the implementation of the vector induction equation in the Boussinesq thermal convection dynamo in a spherical shell model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Model/Boussinesq/Shell/Dynamo/Induction.hpp"

// Project includes
//
#include "QuICC/Base/Typedefs.hpp"
#include "QuICC/Base/MathConstants.hpp"
#include "QuICC/Enums/NonDimensional.hpp"
#include "QuICC/PhysicalOperators/Cross.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Shell {

namespace Dynamo {

   Induction::Induction(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Induction::~Induction()
   {
   }

   void Induction::setCoupling()
   {
      #ifdef QUICC_SPATIALSCHEME_SLFL
         int start = 1;
      #else //if QUICC_SPATIALSCHEME_SLFM
         int start = 0;
      #endif //QUICC_SPATIALSCHEME_SLFL

      this->defineCoupling(FieldComponents::Spectral::TOR, CouplingInformation::PROGNOSTIC, start, true, false);

      this->defineCoupling(FieldComponents::Spectral::POL, CouplingInformation::PROGNOSTIC, start, true, false);
   }

   void Induction::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::POL,0);

      this->addNLComponent(FieldComponents::Spectral::TOR,1);
   }

   void Induction::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      ///
      /// Compute \f$\left(\vec u\wedge\vec B\right)\f$
      ///
      switch(compId)
      {
         case(FieldComponents::Physical::R):
            Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>::set(rNLComp, this->unknown().dom(0).phys(), this->vector(PhysicalNames::VELOCITY).dom(0).phys(), 1.0);
            break;
         case(FieldComponents::Physical::THETA):
            Physical::Cross<FieldComponents::Physical::PHI,FieldComponents::Physical::R>::set(rNLComp, this->unknown().dom(0).phys(), this->vector(PhysicalNames::VELOCITY).dom(0).phys(), 1.0);
            break;
         case(FieldComponents::Physical::PHI):
            Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::set(rNLComp, this->unknown().dom(0).phys(), this->vector(PhysicalNames::VELOCITY).dom(0).phys(), 1.0);
            break;
         default:
            assert(false);
            break;
      }
   }

   void Induction::setRequirements()
   {
      // Set velocity as equation unknown
      this->setName(PhysicalNames::MAGNETIC);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(PhysicalNames::MAGNETIC, FieldRequirement(false, true, true, false, false));

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, true, false, false));
   }

}
}
}
}
}
