#ifndef itkGMMScalePointSetMetric_hxx
#define itkGMMScalePointSetMetric_hxx

#include "itkGMMScalePointSetMetric.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedPointSet, typename TMovingPointSet>
GMMScalePointSetMetric<TFixedPointSet, TMovingPointSet>::GMMScalePointSetMetric()
{
  m_NumberOfParameters = 1;
}

/** Initialize the metric */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMScalePointSetMetric< TFixedPointSet, TMovingPointSet >
::Initialize(void)
throw (ExceptionObject)
{
  if (!m_PointSetMetric) 
  {
    itkExceptionMacro(<< "Point set metric is not present.");
  }

  if (!m_PointSetMetric->IsInitialized()) 
  {
    itkExceptionMacro(<< "Point set metric is presented but is not initialized.");
  }

  m_PointSetMetric->InitializeForIteration(m_PointSetMetric->GetTransform()->GetParameters());
}

/**
* Get the match Measure
*/
template <typename TFixedPointSet, typename TMovingPointSet>
typename GMMScalePointSetMetric<TFixedPointSet, TMovingPointSet>::MeasureType
GMMScalePointSetMetric<TFixedPointSet, TMovingPointSet>::GetValue(const ParametersType & parameters) const
{
  typename PointSetMetricType::DerivativeType derivative;

  m_PointSetMetric->SetScale(parameters[0]);
  m_PointSetMetric->GetDerivative(m_PointSetMetric->GetTransform()->GetParameters(), derivative);

  MeasureType value = NumericTraits<MeasureType>::ZeroValue();

  for (size_t par = 0; par < m_PointSetMetric->GetNumberOfParameters(); ++par) 
  {
    value -= pow(derivative[par], 2);
  }

  return value;
}

/*
* Get both the match Measure and the Derivative Measure
*/
template <typename TFixedPointSet, typename TMovingPointSet>
void
GMMScalePointSetMetric<TFixedPointSet, TMovingPointSet>
::GetValueAndDerivative(const ParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  typename PointSetMetricType::DerivativeType derivative1;
  typename PointSetMetricType::DerivativeType derivative2;

  m_PointSetMetric->SetScale(parameters[0]);
  m_PointSetMetric->GetDerivatives(m_PointSetMetric->GetTransform()->GetParameters(), derivative1, derivative2);

  value = NumericTraits<DerivativeValueType>::ZeroValue();

  if (derivative.size() != this->m_NumberOfParameters) 
  {
    derivative.set_size(this->m_NumberOfParameters);
  }
  
  derivative.Fill(NumericTraits<DerivativeValueType>::ZeroValue());

  for (size_t n = 0; n < m_PointSetMetric->GetNumberOfParameters(); ++n) 
  {
    value -= pow(derivative1[n], 2);
    derivative[0] -= 2 * derivative1[n] * derivative2[n];
  }
}

/**
* Get the Derivative Measure
*/
template <typename TFixedPointSet, typename TMovingPointSet>
void 
GMMScalePointSetMetric<TFixedPointSet, TMovingPointSet>
::GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const
{
  itkExceptionMacro(<< "not implemented");
}

}

#endif
