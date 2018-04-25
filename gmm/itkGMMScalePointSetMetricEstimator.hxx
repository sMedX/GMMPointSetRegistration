#ifndef itkGMMScalePointSetMetricEstimator_hxx
#define itkGMMScalePointSetMetricEstimator_hxx

#include "itkGMMScalePointSetMetricEstimator.h"
#include "itkLBFGSBOptimizer.h"
#include "itkParticleSwarmOptimizer.h"

namespace itk
{
/**
 * Constructor
 */
template< typename TMetricType >
GMMScalePointSetMetricEstimator< TMetricType >::GMMScalePointSetMetricEstimator()
{
  m_Trace = false;
  m_UseInitialParameters = false;
}

/**
*
*/
template< typename TMetricType >
void
GMMScalePointSetMetricEstimator< TMetricType >
::SetInitialParameters(const ParametersType & parameters)
{
  m_InitialParameters = parameters;
  this->SetUseInitialParameters(true);
}

/**
*
*/
template< typename TMetricType >
void
GMMScalePointSetMetricEstimator< TMetricType >
::Estimate()
{
  if (!m_Metric) 
  {
    itkExceptionMacro(<< "Point set metric is not presented.");
  }

  if (!m_Metric->IsInitialized()) 
  {
    itkExceptionMacro(<< "Point set metric is presented but is not initialized.");
  }

  if (!m_UseInitialParameters) 
  {
    m_InitialParameters.set_size(1);

    this->ComputeMeanDistance();
    m_InitialParameters.Fill(m_RMSDistance/2);
  }

  itk::Array<long> boundSelection(1);
  boundSelection.Fill(1);

  itk::Array<double> lowerBound(1);
  lowerBound.Fill(0);

  itk::Array<double> upperBound(1);
  upperBound.Fill(0);

  typedef itk::GMMScalePointSetMetric<FixedPointSetType, MovingPointSetType> GMMScalePointSetMetric;
  GMMScalePointSetMetric::Pointer metric = GMMScalePointSetMetric::New();
  metric->SetPointSetMetric(m_Metric);
  metric->Initialize();

  typedef itk::LBFGSBOptimizer OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetCostFunction(metric);
  optimizer->SetMaximumNumberOfIterations(10);
  optimizer->SetMaximumNumberOfEvaluations(50);
  optimizer->SetBoundSelection(boundSelection);
  optimizer->SetLowerBound(lowerBound);
  optimizer->SetUpperBound(upperBound);
  optimizer->SetInitialPosition(m_InitialParameters);
  optimizer->SetTrace(m_Trace);
  optimizer->SetMinimize(true);
  optimizer->SetCostFunctionConvergenceFactor(1.0);
  optimizer->SetProjectedGradientTolerance(1.0e-05);
  optimizer->StartOptimization();

  m_Parameters = optimizer->GetCurrentPosition();

  if (m_Trace)
  {
    std::cout << optimizer->GetStopConditionDescription() << std::endl;
    std::cout << "Initial metric parameters   " << m_InitialParameters << std::endl;
    std::cout << "Estimated metric parameters " << m_Parameters << std::endl;
    std::cout << std::endl;
  }
}

/**
*
*/
template< typename TMetricType >
void
GMMScalePointSetMetricEstimator< TMetricType >
::ComputeMeanDistance()
{
  typename MetricType::FixedPointSetType::ConstPointer fixedPointSet = m_Metric->GetFixedPointSet();
  typename MetricType::MovingPointSetType::ConstPointer movingPointSet = m_Metric->GetMovingPointSet();

  typename MetricType::TransformType::ConstPointer transform = m_Metric->GetTransform();

  m_RMSDistance = 0;

  for (MovingPointIterator movingIt = movingPointSet->GetPoints()->Begin(); movingIt != movingPointSet->GetPoints()->End(); ++movingIt) {
    MovingPointSetType::PointType movingPoint = transform->TransformPoint(movingIt.Value());

    for (FixedPointIterator fixedIt = fixedPointSet->GetPoints()->Begin(); fixedIt != fixedPointSet->GetPoints()->End(); ++fixedIt) {
      m_RMSDistance += movingPoint.SquaredEuclideanDistanceTo(fixedIt.Value());
    }
  }

  m_RMSDistance = sqrt(m_RMSDistance / (fixedPointSet->GetNumberOfPoints() * movingPointSet->GetNumberOfPoints()));
}
} // end namespace itk
#endif
