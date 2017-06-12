#ifndef itkGMMKCPointSetToPointSetMetric_hxx
#define itkGMMKCPointSetToPointSetMetric_hxx

#include "itkGMMKCPointSetToPointSetMetric.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedPointSet, typename TMovingPointSet>
GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GMMKCPointSetToPointSetMetric()
{
}

/** Initialize the metric */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMKCPointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::Initialize(void)
throw (ExceptionObject)
{
  Superclass::Initialize();

  m_Gradient1.set_size(m_MovingPointSet->GetNumberOfPoints(), MovingPointSetDimension);
  m_Gradient2.set_size(m_MovingPointSet->GetNumberOfPoints(), MovingPointSetDimension);
}

/**
 * Get the match Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
typename GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::MeasureType
GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValue(const TransformParametersType & parameters) const
{
  itkExceptionMacro(<< "not implemented");
}
/**
 * Get the Derivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetDerivative(const TransformParametersType & parameters, DerivativeType & derivative) const
{
  itkExceptionMacro(<< "not implemented");
}

/*
 * Get both the match Measure and theDerivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValueAndDerivative(const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  m_Transform->SetParameters(parameters);

  if (derivative.size() != m_NumberOfParameters) {
    derivative.set_size(m_NumberOfParameters);
  }

  for (MovingPointIterator it = m_MovingPointSet->GetPoints()->Begin(); it != m_MovingPointSet->GetPoints()->End(); ++it) {
    const size_t row = it.Index();
    const typename MovingPointSetType::PointType transformedPoint = m_Transform->TransformPoint(it.Value());

    for (size_t n = 0; n < Self::MovingPointSetDimension; ++n) {
      m_TransformedPointMatrix(row, n) = transformedPoint[n];
    }
  }

  MeasureType value1 = GaussTransform(m_TransformedPointMatrix.data_block(),
    m_FixedPointMatrix.data_block(),
    m_TransformedPointMatrix.rows(),
    m_FixedPointMatrix.rows(),
    m_TransformedPointMatrix.cols(),
    m_MovingPointSetScale,
    m_FixedPointSetScale,
    m_Gradient1.data_block());

  MeasureType value2 = GaussTransform(m_TransformedPointMatrix.data_block(),
    m_TransformedPointMatrix.data_block(),
    m_TransformedPointMatrix.rows(),
    m_TransformedPointMatrix.rows(),
    m_TransformedPointMatrix.cols(),
    m_MovingPointSetScale,
    m_MovingPointSetScale,
    m_Gradient2.data_block());

  double ratio = value1 / value2;
  value = -value1 * ratio;

  for (MovingPointIterator it = m_MovingPointSet->GetPoints()->Begin(); it != m_MovingPointSet->GetPoints()->End(); ++it) {
    size_t row = it.Index();

    for (size_t dim = 0; dim < Self::MovingPointSetDimension; ++dim) {
      m_Gradient(row, dim) = 2.0 * ratio * (m_Gradient2(row, dim) * ratio - m_Gradient1(row, dim));
    }
  }

  // compute the derivatives
  derivative.Fill(NumericTraits<typename DerivativeType::ValueType>::ZeroValue());

  for (MovingPointIterator it = m_MovingPointSet->GetPoints()->Begin(); it != m_MovingPointSet->GetPoints()->End(); ++it) {
    size_t row = it.Index();
    m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(it.Value(), m_Jacobian, m_JacobianCache);

    for (size_t par = 0; par < m_NumberOfParameters; par++) {
      double sum = 0;

      for (size_t dim = 0; dim < Self::MovingPointSetDimension; dim++) {
        sum += m_Jacobian(dim, par) * m_Gradient(row, dim);
      }

      derivative[par] += sum;
    }
  }
}
}

#endif
