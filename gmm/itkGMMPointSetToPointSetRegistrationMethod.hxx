#ifndef itkGMMPointSetToPointSetRegistrationMethod_hxx
#define itkGMMPointSetToPointSetRegistrationMethod_hxx

#include "itkGMMPointSetToPointSetRegistrationMethod.h"

namespace itk
{
/**
 * Constructor
 */
template< typename TFixedPointSet, typename TMovingPointSet >
GMMPointSetToPointSetRegistrationMethod< TFixedPointSet, TMovingPointSet >::GMMPointSetToPointSetRegistrationMethod()
{
  this->SetNumberOfRequiredOutputs(1);    // for the transform

  m_FixedPointSet = ITK_NULLPTR;
  m_MovingPointSet = ITK_NULLPTR;

  m_Transform = ITK_NULLPTR;
  m_Metric = ITK_NULLPTR;
  m_Optimizer = ITK_NULLPTR;

  m_NumberOfLevels = 1;
  m_GradientConvergenceTolerance = 1e-5;

  m_InitialTransformParameters = ParametersType(1);
  m_FinalTransformParameters = ParametersType(1);

  m_InitialTransformParameters.Fill(0);
  m_FinalTransformParameters.Fill(0);

  TransformOutputPointer transformDecorator = itkDynamicCastInDebugMode< TransformOutputType * >(this->MakeOutput(0).GetPointer() );
  this->ProcessObject::SetNthOutput(0, transformDecorator.GetPointer());
}

/*
 * Set the initial transform parameters
 */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMPointSetToPointSetRegistrationMethod< TFixedPointSet, TMovingPointSet >
::SetInitialTransformParameters(const ParametersType & param)
{
  m_InitialTransformParameters = param;
  this->Modified();
}

/**
 * Initialize by setting the interconnects between components.
 */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMPointSetToPointSetRegistrationMethod< TFixedPointSet, TMovingPointSet >
::Initialize() throw (ExceptionObject)
{
  if (!m_FixedPointSet) 
  {
    itkExceptionMacro(<< "FixedPointSet is not present");
  }

  if (!m_MovingPointSet) 
  {
    itkExceptionMacro(<< "MovingPointSet is not present");
  }

  if (!m_Metric) 
  {
    itkExceptionMacro(<< "Metric is not present");
  }

  if (!m_Optimizer) 
  {
    itkExceptionMacro(<< "Optimizer is not present");
  }

  if (!m_Transform) 
  {
    itkExceptionMacro(<< "Transform is not present");
  }

  // Validate initial transform parameters
  if (m_InitialTransformParameters.Size() != m_Transform->GetNumberOfParameters()) {
    m_InitialTransformParameters = m_Transform->GetParameters();
  }

  // Connect the transform to the Decorator.
  TransformOutputType *transformOutput = static_cast<TransformOutputType*>(this->ProcessObject::GetOutput(0));
  transformOutput->Set(m_Transform.GetPointer());

  // setup the transform
  m_Transform->SetParameters(m_InitialTransformParameters);

  m_InitialMetricValues.clear();
  m_InitialMetricValues.set_size(m_NumberOfLevels);
  m_InitialMetricValues.Fill(NAN);

  m_FinalMetricValues.clear();
  m_FinalMetricValues.set_size(m_NumberOfLevels);
  m_FinalMetricValues.Fill(NAN);

  m_Scales.clear();
  m_Scales.set_size(m_NumberOfLevels);
  m_Scales.Fill(NAN);
}

/**
 *  Get Output
 */
template< typename TFixedPointSet, typename TMovingPointSet >
const typename GMMPointSetToPointSetRegistrationMethod< TFixedPointSet, TMovingPointSet >::TransformOutputType *
GMMPointSetToPointSetRegistrationMethod< TFixedPointSet, TMovingPointSet >
::GetOutput() const
{
  return static_cast< const TransformOutputType * >( this->ProcessObject::GetOutput(0) );
}

template< typename TFixedPointSet, typename TMovingPointSet >
DataObject::Pointer
GMMPointSetToPointSetRegistrationMethod< TFixedPointSet, TMovingPointSet >
::MakeOutput(DataObjectPointerArraySizeType output)
{
  switch ( output )
    {
    case 0:
      return TransformOutputType::New().GetPointer();
      break;
    default:
      itkExceptionMacro("MakeOutput request for an output number larger than the expected number of outputs");
      return ITK_NULLPTR;
    }
}

/**
 *
 */
template< typename TFixedPointSet, typename TMovingPointSet >
ModifiedTimeType
GMMPointSetToPointSetRegistrationMethod< TFixedPointSet, TMovingPointSet >
::GetMTime() const
{
  ModifiedTimeType mtime = Superclass::GetMTime();
  ModifiedTimeType m;

  // Some of the following should be removed once ivars are put in the input and output lists
  if ( m_Transform )
    {
    m = m_Transform->GetMTime();
    mtime = ( m > mtime ? m : mtime );
    }

  if ( m_Metric )
    {
    m = m_Metric->GetMTime();
    mtime = ( m > mtime ? m : mtime );
    }

  if ( m_Optimizer )
    {
    m = m_Optimizer->GetMTime();
    mtime = ( m > mtime ? m : mtime );
    }

  if ( m_FixedPointSet )
    {
    m = m_FixedPointSet->GetMTime();
    mtime = ( m > mtime ? m : mtime );
    }

  if ( m_MovingPointSet )
    {
    m = m_MovingPointSet->GetMTime();
    mtime = ( m > mtime ? m : mtime );
    }

  return mtime;
}

/**
*
*/
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMPointSetToPointSetRegistrationMethod< TFixedPointSet, TMovingPointSet >
::GenerateData()
{
  try {
    // initialize the interconnects between components
    this->Initialize();
  }
  catch (ExceptionObject & excep) {
    m_FinalTransformParameters = ParametersType(1);
    m_FinalTransformParameters.Fill(0);

    // pass exception to caller
    throw excep;
  }

  // setup the metric
  m_Metric->SetFixedPointSet(m_FixedPointSet);
  m_Metric->SetMovingPointSet(m_MovingPointSet);
  m_Metric->SetTransform(m_Transform);
  m_Metric->Initialize();

  // setup the optimizer and metric estimator
  m_Optimizer->SetCostFunction(m_Metric);

  typedef typename FixedPointSetType::PointsContainer::ConstIterator   FixedPointIterator;
  typedef typename MovingPointSetType::PointsContainer::ConstIterator  MovingPointIterator;

  double RMSDistance = 0;

  for (MovingPointIterator movingIt = m_MovingPointSet->GetPoints()->Begin(); movingIt != m_MovingPointSet->GetPoints()->End(); ++movingIt) {
    MovingPointSetType::PointType movingPoint = m_Transform->TransformPoint(movingIt.Value());

    for (FixedPointIterator fixedIt = m_FixedPointSet->GetPoints()->Begin(); fixedIt != m_FixedPointSet->GetPoints()->End(); ++fixedIt) {
      RMSDistance += movingPoint.SquaredEuclideanDistanceTo(fixedIt.Value());
    }
  }

  RMSDistance = sqrt(RMSDistance / (m_FixedPointSet->GetNumberOfPoints() * m_MovingPointSet->GetNumberOfPoints()));
  m_Scales[0] = RMSDistance;

  for (size_t level = 0; level < m_NumberOfLevels; ++level) 
  {
    if (level == 0) {
      m_Metric->SetScale(m_Scales[level]);
    }
    else {
      typename MetricType::DerivativeType derivative;

      m_Scales[level] = m_Scales[level - 1] / 2;

      m_Metric->SetScale(m_Scales[level]);
      m_Metric->GetDerivative(m_FinalTransformParameters, derivative);

      // The optimization terminates when : || G || < gtol max(1, || X || ) where || . || denotes the Euclidean norm.
      if (derivative.two_norm() < m_GradientConvergenceTolerance * std::max(ParametersValueType(1), m_FinalTransformParameters.two_norm()))
        break;
    }

    m_Optimizer->SetInitialPosition(m_Transform->GetParameters());
    try {
      m_Optimizer->StartOptimization();
    }
    catch (ExceptionObject & excep) {
      std::cout << excep << std::endl;

      m_FinalTransformParameters = m_Optimizer->GetCurrentPosition();
      m_Transform->SetParameters(m_FinalTransformParameters);

      m_InitialMetricValues[level] = m_Metric->GetValue(m_Optimizer->GetInitialPosition());
      m_FinalMetricValues[level] = m_Metric->GetValue(m_Optimizer->GetCurrentPosition());
      break;
    }

    m_FinalTransformParameters = m_Optimizer->GetCurrentPosition();
    m_Transform->SetParameters(m_FinalTransformParameters);

    m_InitialMetricValues[level] = m_Metric->GetValue(m_Optimizer->GetInitialPosition());
    m_FinalMetricValues[level] = m_Metric->GetValue(m_Optimizer->GetCurrentPosition());
  }
}
} // end namespace itk
#endif
