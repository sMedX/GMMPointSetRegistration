#ifndef itkGMMMLEPointSetToPointSetMetric_hxx
#define itkGMMMLEPointSetToPointSetMetric_hxx

#include "itkGMMMLEPointSetToPointSetMetric.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedPointSet, typename TMovingPointSet>
GMMMLEPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GMMMLEPointSetToPointSetMetric()
{
}

/** Initialize the metric */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMMLEPointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::Initialize() throw (ExceptionObject)
{
  Superclass::Initialize();
  m_Gradient.set_size(this->m_MovingPointSet->GetNumberOfPoints(), this->MovingPointSetDimension);
}

/**
 * Get the match Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
typename GMMMLEPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::MeasureType
GMMMLEPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValue(const TransformParametersType & parameters) const
{
  itkExceptionMacro(<< "not implemented");
}
/**
 * Get the Derivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void GMMMLEPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetDerivative(const TransformParametersType & parameters, DerivativeType & derivative) const
{
  itkExceptionMacro(<< "not implemented");
}

/*
 * Get both the match Measure and the Derivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void GMMMLEPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValueAndDerivative(const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  this->m_Transform->SetParameters(parameters);

  if (derivative.size() != this->m_NumberOfParameters) {
    derivative.set_size(this->m_NumberOfParameters);
  }

  // compute transformed point set
  this->m_TransformedPointSet = MovingPointSetType::New();
  for (MovingPointIterator iter = this->m_MovingPointSet->GetPoints()->Begin(); iter != this->m_MovingPointSet->GetPoints()->End(); ++iter) {
    this->m_TransformedPointSet->SetPoint(iter.Index(), this->m_Transform->TransformPoint(iter.Value()));
  }

  value = NumericTraits<MeasureType>::ZeroValue();
  double scale = 2*this->m_MovingPointSetScale*this->m_MovingPointSetScale;

  itk::Array<double> valuesOfProbability;
  valuesOfProbability.set_size(this->m_FixedPointSet->GetNumberOfPoints());

  for (FixedPointIterator fixedIter = this->m_FixedPointSet->GetPoints()->Begin(); fixedIter != this->m_FixedPointSet->GetPoints()->End(); ++fixedIter) {
    const typename FixedPointSetType::PointType fixedPoint = fixedIter.Value();
    double sum = 1.0e-05;

    for (MovingPointIterator movingIter = this->m_TransformedPointSet->GetPoints()->Begin(); movingIter != this->m_TransformedPointSet->GetPoints()->End(); ++movingIter) {
      const typename MovingPointSetType::PointType transformedPoint = movingIter.Value();
      const double distance = transformedPoint.SquaredEuclideanDistanceTo(fixedPoint);
      sum += exp(-distance / scale);
    }

    valuesOfProbability[fixedIter.Index()] = sum;
    value -= log(sum);
  }

  // compute the gradient
  this->m_Gradient.fill(0);

  for (MovingPointIterator movingIter = this->m_TransformedPointSet->GetPoints()->Begin(); movingIter != this->m_TransformedPointSet->GetPoints()->End(); ++movingIter) {
    const typename MovingPointSetType::PointType transformedPoint = movingIter.Value();
    
    for (FixedPointIterator fixedIter = this->m_FixedPointSet->GetPoints()->Begin(); fixedIter != this->m_FixedPointSet->GetPoints()->End(); ++fixedIter) {
      const typename FixedPointSetType::PointType fixedPoint = fixedIter.Value();
      const double distance = transformedPoint.SquaredEuclideanDistanceTo(fixedPoint);
      const double expval = exp(-distance / scale);
      const size_t row = fixedIter.Index();

      for (size_t dim = 0; dim < Self::MovingPointSetDimension; ++dim) {
        this->m_Gradient(row, dim) += 2.0 * expval * (transformedPoint[dim] - fixedPoint[dim]) / scale / valuesOfProbability[row];
      }
    }
  }

  // compute the derivatives
  derivative.Fill(NumericTraits<typename DerivativeType::ValueType>::ZeroValue());

  for (MovingPointIterator movingIter = this->m_MovingPointSet->GetPoints()->Begin(); movingIter != this->m_MovingPointSet->GetPoints()->End(); ++movingIter) {
    this->m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(movingIter.Value(), this->m_Jacobian, this->m_JacobianCache);
    const size_t row = movingIter.Index();

    for (size_t par = 0; par < this->m_NumberOfParameters; ++par) {
      double sum = 0;
      for (size_t dim = 0; dim < this->PointDimension; ++dim) {
        sum += this->m_Jacobian(dim, par) * this->m_Gradient(row, dim);
      }
      derivative[par] = sum;
    }
  }
}
}

#endif
