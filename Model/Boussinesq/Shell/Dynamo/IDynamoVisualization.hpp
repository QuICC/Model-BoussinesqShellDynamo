/**
 * @file IDynamoVisualization.hpp
 * @brief Implementation of the Boussinesq thermal convection dynamo in a
 * spherical shell model (Toroidal/Poloidal formulation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SHELL_DYNAMO_IDYNAMOVISUALIZATION_HPP
#define QUICC_MODEL_BOUSSINESQ_SHELL_DYNAMO_IDYNAMOVISUALIZATION_HPP

// System includes
//
#include <string>

// Project includes
//
#include "QuICC/Generator/VisualizationGenerator.hpp"
#include "QuICC/Model/IVisualizationGeneratorBuilder.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace Dynamo {

/**
 * @brief Implementation of the Boussinesq thermal convection dynamo in a
 * spherical shell model (Toroidal/Poloidal formulation)
 */
class IDynamoVisualization : public IVisualizationGeneratorBuilder<VisualizationGenerator>
{
public:
   /**
    * @brief Constructor
    */
   IDynamoVisualization() = default;

   /**
    * @brief Destructor
    */
   virtual ~IDynamoVisualization() = default;

   /// Formulation used for vector fields
   virtual VectorFormulation::Id SchemeFormulation() override;

   /**
    * @brief Version string
    */
   std::string version() const final;

   /**
    * @brief Add the visualization generation equations
    *
    * @param spGen   Shared visualization generator
    */
   virtual void addVisualizers(SharedVisualizationGenerator spVis) override;

protected:
private:
};

} // namespace Dynamo
} // namespace Shell
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SHELL_DYNAMO_IDYNAMOVISUALIZATION_HPP
