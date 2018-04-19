#ifndef itkGMMSigmaMetric_hxx
#define itkGMMSigmaMetric_hxx

#include "itkGMMSigmaMetric.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedPointSet, typename TMovingPointSet>
GMMSigmaMetric<TFixedPointSet, TMovingPointSet>::GMMSigmaMetric()
{
  m_NumberOfParameters = 1;
}

/**
* Get the match Measure
*/
template <typename TFixedPointSet, typename TMovingPointSet>
typename GMMSigmaMetric<TFixedPointSet, TMovingPointSet>::MeasureType
GMMSigmaMetric<TFixedPointSet, TMovingPointSet>::GetValue(const ParametersType & parameters) const
{
  typename MetricType::DerivativeType derivative;
  MeasureType value = NumericTraits<MeasureType>::ZeroValue();

  m_Metric->SetScale(parameters[0]);
  m_Metric->GetDerivative(m_Metric->GetTransform()->GetParameters(), derivative);

  for (size_t par = 0; par < m_Metric->GetNumberOfParameters(); ++par) 
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
GMMSigmaMetric<TFixedPointSet, TMovingPointSet>
::GetValueAndDerivative(const ParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  if (derivative.size() != this->m_NumberOfParameters) 
  {
    derivative.set_size(this->m_NumberOfParameters);
  }

  value = NumericTraits<DerivativeValueType>::ZeroValue();
  derivative.Fill(NumericTraits<DerivativeValueType>::ZeroValue());

  typename MetricType::DerivativeType derivative1;
  typename MetricType::DerivativeType derivative2;

  m_Metric->SetScale(parameters[0]);
  m_Metric->GetDerivatives(m_Metric->GetTransform()->GetParameters(), derivative1, derivative2);

  for (size_t n = 0; n < m_Metric->GetNumberOfParameters(); ++n) 
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
GMMSigmaMetric<TFixedPointSet, TMovingPointSet>
::GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const
{
  itkExceptionMacro(<< "not implemented");
}

}

#endif
