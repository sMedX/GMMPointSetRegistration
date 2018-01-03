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
  m_FixedTransformedPointSet = ITK_NULLPTR;
  m_FixedInitialTransform = ITK_NULLPTR;
  m_MovingPointSet = ITK_NULLPTR;
  m_MovingTransformedPointSet = ITK_NULLPTR;
  m_MovingInitialTransform = ITK_NULLPTR;

  m_Transform = ITK_NULLPTR;
  m_Metric = ITK_NULLPTR;
  m_Optimizer = ITK_NULLPTR;

  m_NumberOfLevels = 0;

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
  if (!m_FixedPointSet) {
    itkExceptionMacro(<< "FixedPointSet is not present");
  }

  if (!m_MovingPointSet) {
    itkExceptionMacro(<< "MovingPointSet is not present");
  }

  if (!m_Metric) {
    itkExceptionMacro(<< "Metric is not present");
  }

  if (!m_Optimizer) {
    itkExceptionMacro(<< "Optimizer is not present");
  }

  if (!m_Transform) {
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

  if (m_FixedPointSetScale.size() != m_MovingPointSetScale.size()) {
    itkExceptionMacro(<< "Sizes of m_FixedPointSetScale and m_MovingPointSetScale are mismatched");
  }

  if (m_NumberOfLevels > m_FixedPointSetScale.size()) {
    itkExceptionMacro(<< "The number of levels is too large");
  }

  if (m_NumberOfLevels == 0) {
    m_NumberOfLevels = m_FixedPointSetScale.Size();
  }

  m_InitialMetricValue.clear();
  m_FinalMetricValue.clear();
}

template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMPointSetToPointSetRegistrationMethod< TFixedPointSet, TMovingPointSet >
::Preprocessing() throw (ExceptionObject)
{
  if ( m_FixedInitialTransform )
  {
    typename FixedPointsContainerType::Pointer points = FixedPointsContainerType::New();
    m_FixedTransformedPointSet = FixedPointSetType::New();

    for (FixedPointConstIterator it = m_FixedPointSet->GetPoints()->Begin(); it != m_FixedPointSet->GetPoints()->End(); ++it) {
      points->InsertElement(it.Index(), m_FixedInitialTransform->TransformPoint(it.Value()));
    }

    m_FixedTransformedPointSet->SetPoints(points);
  }

  if ( m_MovingInitialTransform ) 
  {
    typename MovingPointsContainerType::Pointer points = MovingPointsContainerType::New();
    m_MovingTransformedPointSet = MovingPointSetType::New();

    for (MovingPointConstIterator it = m_MovingPointSet->GetPoints()->Begin(); it != m_MovingPointSet->GetPoints()->End(); ++it) {
      points->InsertElement(it.Index(), m_MovingInitialTransform->TransformPoint(it.Value()));
    }

    m_MovingTransformedPointSet->SetPoints(points);
  }
}

/**
 * PrintSelf
 */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMPointSetToPointSetRegistrationMethod< TFixedPointSet, TMovingPointSet >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Metric: " << m_Metric.GetPointer() << std::endl;
  os << indent << "Optimizer: " << m_Optimizer.GetPointer() << std::endl;
  os << indent << "Transform: " << m_Transform.GetPointer() << std::endl;
  os << indent << " Fixed PointSet: " << m_FixedPointSet.GetPointer() << std::endl;
  os << indent << "Moving PointSet: " << m_MovingPointSet.GetPointer() << std::endl;
  os << indent << "Initial Transform Parameters: " << m_InitialTransformParameters << std::endl;
  os << indent << "  Final Transform Parameters: " << m_FinalTransformParameters << std::endl;

  if (m_InitialMetricValue.size() > 0) {
    os << indent << "Initial metric values ";
    for (size_t n = 0; n < m_InitialMetricValue.size(); ++n) {
      os << m_InitialMetricValue[n] << " ";
    }
    os << std::endl;
  }

  if (m_FinalMetricValue.size() > 0) {
    os << indent << "  Final metric values ";
    for (size_t n = 0; n < m_FinalMetricValue.size(); ++n) {
      os << m_FinalMetricValue[n] << " ";
    }
    os << std::endl;
  }
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
    m_FinalTransformParameters.Fill(0.0f);

    // pass exception to caller
    throw excep;
  }

  try {
    // perform preprocessing of fixed and point sets
    this->Preprocessing();
  }
  catch (ExceptionObject & excep) {
    // pass exception to caller
    throw excep;
  }

  // setup the metric
  if (m_FixedTransformedPointSet) {
    m_Metric->SetFixedPointSet(m_FixedTransformedPointSet);
  }
  else {
    m_Metric->SetFixedPointSet(m_FixedPointSet);
  }

  if (m_MovingTransformedPointSet) {
    m_Metric->SetMovingPointSet(m_MovingTransformedPointSet);
  }
  else {
    m_Metric->SetMovingPointSet(m_MovingPointSet);
  }

  m_Metric->SetTransform(m_Transform);
  m_Metric->Initialize();

  // setup the optimizer
  m_Optimizer->SetCostFunction(m_Metric);

  for (size_t level = 0; level < m_NumberOfLevels; ++level) {
    m_Metric->SetFixedPointSetScale(m_FixedPointSetScale[level]);
    m_Metric->SetMovingPointSetScale(m_MovingPointSetScale[level]);

    m_InitialMetricValue;
    m_Optimizer->SetInitialPosition(m_Transform->GetParameters());

    try {
      m_Optimizer->StartOptimization();
    }
    catch (ExceptionObject & excep) {
      std::cout << excep << std::endl;
    }

    // get the results
    m_FinalTransformParameters = m_Optimizer->GetCurrentPosition();
    m_Transform->SetParameters(m_FinalTransformParameters);

    m_InitialMetricValue.push_back(m_Metric->GetValue(m_Optimizer->GetInitialPosition()));
    m_FinalMetricValue.push_back(m_Metric->GetValue(m_Optimizer->GetCurrentPosition()));
  }
}
} // end namespace itk
#endif
