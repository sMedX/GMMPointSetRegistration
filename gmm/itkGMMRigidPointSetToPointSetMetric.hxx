#ifndef itkGMMRigidPointSetToPointSetMetric_hxx
#define itkGMMRigidPointSetToPointSetMetric_hxx

#include "itkGMMRigidPointSetToPointSetMetric.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedPointSet, typename TMovingPointSet>
GMMRigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GMMRigidPointSetToPointSetMetric()
{
}

/**
 * Get the match Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
typename GMMRigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::MeasureType
GMMRigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValue(const TransformParametersType & parameters) const
{

	m_Transform->SetParameters(parameters);

  for (MovingPointIterator it = m_MovingPointSet->GetPoints()->Begin(); it != m_MovingPointSet->GetPoints()->End(); ++it) {
    size_t row = it.Index();
    const typename MovingPointSetType::PointType transformedPoint = m_Transform->TransformPoint(it.Value());

    for (size_t n = 0; n < Self::MovingPointSetDimension; ++n) {
      m_TransformedPointMatrix(row, n) = transformedPoint[n];
    }
  }

  double value = GaussTransform(m_TransformedPointMatrix.data_block(),
    m_FixedPointMatrix.data_block(),
    m_TransformedPointMatrix.rows(),
    m_FixedPointMatrix.rows(),
    m_TransformedPointMatrix.cols(),
    m_MovingPointSetScale, m_FixedPointSetScale,
    m_Gradient.data_block());

  value *= -1;
  m_Gradient *= -1;

  return value;
}

/**
 * Get the Derivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void GMMRigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetDerivative(const TransformParametersType & parameters, DerivativeType & derivative) const
{

}

/*
 * Get both the match Measure and theDerivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void GMMRigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValueAndDerivative(const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  m_Transform->SetParameters(parameters);

  if (derivative.size() != m_NumberOfParameters) {
    derivative.set_size(m_NumberOfParameters);
  }

  for (MovingPointIterator it = m_MovingPointSet->GetPoints()->Begin(); it != m_MovingPointSet->GetPoints()->End(); ++it) {
    size_t row = it.Index();
    const typename MovingPointSetType::PointType transformedPoint = m_Transform->TransformPoint(it.Value());

    for (size_t n = 0; n < Self::MovingPointSetDimension; ++n) {
      m_TransformedPointMatrix(row, n) = transformedPoint[n];
    }
  }

  value = GaussTransform(m_TransformedPointMatrix.data_block(),
    m_FixedPointMatrix.data_block(),
    m_TransformedPointMatrix.rows(),
    m_FixedPointMatrix.rows(),
    m_TransformedPointMatrix.cols(),
    m_MovingPointSetScale,
    m_FixedPointSetScale,
    m_Gradient.data_block());

  value *= -1;
  m_Gradient *= -1;

  // compute the derivatives
  derivative.Fill(NumericTraits<typename DerivativeType::ValueType>::ZeroValue());

  for (MovingPointIterator it = m_MovingPointSet->GetPoints()->Begin(); it != m_MovingPointSet->GetPoints()->End(); ++it) {
    size_t row = it.Index();
    m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(it.Value(), m_Jacobian, m_JacobianCache);

    for (size_t par = 0; par < m_NumberOfParameters; par++) {
      double sum = 0;

      for (size_t dim = 0; dim < Self::FixedPointSetDimension; dim++) {
        sum += m_Jacobian(dim, par) * m_Gradient(row, dim);
      }

      derivative[par] += sum;
    }
  }
}
}

#endif
