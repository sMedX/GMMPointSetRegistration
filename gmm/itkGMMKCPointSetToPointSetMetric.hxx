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
::Initialize() throw (ExceptionObject)
{
  Superclass::Initialize();
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
 * Get both the match Measure and the Derivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValueAndDerivative(const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  this->m_Transform->SetParameters(parameters);

  if (derivative.size() != this->m_NumberOfParameters) {
    derivative.set_size(this->m_NumberOfParameters);
  }

  this->m_TransformedPointSet = MovingPointSetType::New();
  for (MovingPointIterator iter = this->m_MovingPointSet->GetPoints()->Begin(); iter != this->m_MovingPointSet->GetPoints()->End(); ++iter) {
    this->m_TransformedPointSet->SetPoint(iter.Index(), this->m_Transform->TransformPoint(iter.Value()));
  }

  double value1 = 0;
  double value2 = 0;

  GradientType gradient1;
  GradientType gradient2;

  LocalDerivativeType derivative1(this->m_NumberOfParameters);
  derivative1.Fill(NumericTraits<typename DerivativeType::ValueType>::ZeroValue());

  LocalDerivativeType derivative2(this->m_NumberOfParameters);
  derivative2.Fill(NumericTraits<typename DerivativeType::ValueType>::ZeroValue());

  double scale1 = (this->m_FixedPointSetScale*this->m_FixedPointSetScale + this->m_MovingPointSetScale*this->m_MovingPointSetScale) / 2;
  double scale2 = this->m_MovingPointSetScale*this->m_MovingPointSetScale;

  for (MovingPointIterator movingIter1 = this->m_TransformedPointSet->GetPoints()->Begin(); movingIter1 != this->m_TransformedPointSet->GetPoints()->End(); ++movingIter1) {
    const typename MovingPointSetType::PointType transformedPoint1 = movingIter1.Value();

    //------------------------------------
    // compute gradient for the first part
    gradient1.Fill(0);

    for (FixedPointIterator fixedIter = this->m_FixedPointSet->GetPoints()->Begin(); fixedIter != this->m_FixedPointSet->GetPoints()->End(); ++fixedIter) {
      const typename FixedPointSetType::PointType fixedPoint = fixedIter.Value();
      const double distance = transformedPoint1.SquaredEuclideanDistanceTo(fixedPoint);
      const double expval = exp(-distance / scale1);
      value1 += expval;

      for (size_t dim = 0; dim < this->PointDimension; ++dim) {
        gradient1[dim] += (-2.0) * expval * (transformedPoint1[dim] - fixedPoint[dim]) / scale1;
      }
    }

    //-------------------------------------
    // compute gradient for the second part
    gradient2.Fill(0);

    for (MovingPointIterator movingIter2 = this->m_TransformedPointSet->GetPoints()->Begin(); movingIter2 != this->m_TransformedPointSet->GetPoints()->End(); ++movingIter2) {
      const typename MovingPointSetType::PointType transformedPoint2 = movingIter2.Value();
      const double distance = transformedPoint1.SquaredEuclideanDistanceTo(transformedPoint2);
      const double expval = exp(-distance / scale2);
      value2 += expval;

      for (size_t dim = 0; dim < this->PointDimension; ++dim) {
        gradient2[dim] += (-2.0) * expval * (transformedPoint1[dim] - transformedPoint2[dim]) / scale2;
      }
    }

    //-------------------------------------
    // compute derivatives
    this->m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(this->m_MovingPointSet->GetPoint(movingIter1.Index()), this->m_Jacobian, this->m_JacobianCache);

    for (size_t par = 0; par < this->m_NumberOfParameters; par++) {
      for (size_t dim = 0; dim < this->PointDimension; dim++) {
        derivative1[par] += this->m_Jacobian(dim, par) * gradient1[dim];
        derivative2[par] += this->m_Jacobian(dim, par) * gradient2[dim];
      }
    }
  }

  const double factor1 = this->m_TransformedPointSet->GetNumberOfPoints() * this->m_FixedPointSet->GetNumberOfPoints();
  const double factor2 = this->m_TransformedPointSet->GetNumberOfPoints() * this->m_TransformedPointSet->GetNumberOfPoints();

  value1 /= factor1;
  value2 /= factor2;

  const double ratio = value1 / value2;
  value = -value1 * ratio;

  for (size_t par = 0; par < this->m_NumberOfParameters; par++) {
    derivative[par] = (-2.0) * (derivative1[par] / factor1 - derivative2[par] * ratio / factor2) * ratio;
  }
}
}

#endif
