#ifndef itkGMML2PointSetToPointSetMetric_hxx
#define itkGMML2PointSetToPointSetMetric_hxx

#include "itkGMML2PointSetToPointSetMetric.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedPointSet, typename TMovingPointSet>
GMML2PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GMML2PointSetToPointSetMetric()
{
}

/** Initialize the metric */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMML2PointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::Initialize() throw (ExceptionObject)
{
  Superclass::Initialize();

  m_Gradient1.set_size(this->m_MovingPointSet->GetNumberOfPoints(), this->MovingPointSetDimension);
  m_Gradient2.set_size(this->m_MovingPointSet->GetNumberOfPoints(), this->MovingPointSetDimension);
}

/**
 * Get the match Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
typename GMML2PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::MeasureType
GMML2PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValue(const TransformParametersType & parameters) const
{
  itkExceptionMacro(<< "not implemented");
}
/**
 * Get the Derivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void GMML2PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetDerivative(const TransformParametersType & parameters, DerivativeType & derivative) const
{
  itkExceptionMacro(<< "not implemented");
}

/*
 * Get both the match Measure and the Derivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void GMML2PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValueAndDerivative(const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  this->m_Transform->SetParameters(parameters);

  if (derivative.size() != this->m_NumberOfParameters) {
    derivative.set_size(this->m_NumberOfParameters);
  }

  for (MovingPointIterator it = this->m_MovingPointSet->GetPoints()->Begin(); it != this->m_MovingPointSet->GetPoints()->End(); ++it) {
    const size_t row = it.Index();
    const typename MovingPointSetType::PointType transformedPoint = this->m_Transform->TransformPoint(it.Value());

    for (size_t n = 0; n < Self::MovingPointSetDimension; ++n) {
      this->m_TransformedPointMatrix(row, n) = transformedPoint[n];
    }
  }

  MeasureType value1 = GaussTransform(this->m_TransformedPointMatrix.data_block(),
    this->m_FixedPointMatrix.data_block(),
    this->m_TransformedPointMatrix.rows(),
    this->m_FixedPointMatrix.rows(),
    this->m_TransformedPointMatrix.cols(),
    this->m_MovingPointSetScale,
    this->m_FixedPointSetScale,
    m_Gradient1.data_block());

  MeasureType value2 = GaussTransform(this->m_TransformedPointMatrix.data_block(),
    this->m_TransformedPointMatrix.data_block(),
    this->m_TransformedPointMatrix.rows(),
    this->m_TransformedPointMatrix.rows(),
    this->m_TransformedPointMatrix.cols(),
    this->m_MovingPointSetScale,
    this->m_MovingPointSetScale,
    m_Gradient2.data_block());

  value = -2 * value1 + value2;

  for (MovingPointIterator it = this->m_MovingPointSet->GetPoints()->Begin(); it != this->m_MovingPointSet->GetPoints()->End(); ++it) {
    size_t row = it.Index();

    for (size_t dim = 0; dim < Self::MovingPointSetDimension; ++dim) {
      this->m_Gradient(row, dim) = 2.0 * (-m_Gradient1(row, dim) + m_Gradient2(row, dim));
    }
  }

  // compute the derivatives
  derivative.Fill(NumericTraits<typename DerivativeType::ValueType>::ZeroValue());

  for (MovingPointIterator it = this->m_MovingPointSet->GetPoints()->Begin(); it != this->m_MovingPointSet->GetPoints()->End(); ++it) {
    size_t row = it.Index();
    this->m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(it.Value(), this->m_Jacobian, this->m_JacobianCache);

    for (size_t par = 0; par < this->m_NumberOfParameters; par++) {
      double sum = 0;

      for (size_t dim = 0; dim < Self::MovingPointSetDimension; dim++) {
        sum += this->m_Jacobian(dim, par) * this->m_Gradient(row, dim);
      }

      derivative[par] += sum;
    }
  }
}
}

#endif
