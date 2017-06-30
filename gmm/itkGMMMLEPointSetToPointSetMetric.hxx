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
  
  const double scale = 2.0 * this->m_MovingPointSetScale * this->m_MovingPointSetScale;

  itk::Array<double> valuesOfProbability;
  valuesOfProbability.set_size(this->m_NumberOfFixedPoints);

  for (FixedPointIterator fixedIter = this->m_FixedPointSet->GetPoints()->Begin(); fixedIter != this->m_FixedPointSet->GetPoints()->End(); ++fixedIter) {
    const typename FixedPointSetType::PointType fixedPoint = fixedIter.Value();
    double sum = 1.0e-05;

    for (MovingPointIterator transformedIter = this->m_TransformedPointSet->GetPoints()->Begin(); transformedIter != this->m_TransformedPointSet->GetPoints()->End(); ++transformedIter) {
      const typename MovingPointSetType::PointType transformedPoint = transformedIter.Value();
      const double distance = transformedPoint.SquaredEuclideanDistanceTo(fixedPoint);
      sum += exp(-distance / scale);
    }

    valuesOfProbability[fixedIter.Index()] = sum;
    value -= log(sum);
  }

  // compute the derivatives
  derivative.Fill(NumericTraits<typename DerivativeType::ValueType>::ZeroValue());
  GradientType gradient;

  for (MovingPointIterator transformedIter = this->m_TransformedPointSet->GetPoints()->Begin(); transformedIter != this->m_TransformedPointSet->GetPoints()->End(); ++transformedIter) {

    // compute gradient for the current transformed point
    const typename MovingPointSetType::PointType transformedPoint = transformedIter.Value();
    gradient.Fill(0);

    for (FixedPointIterator fixedIter = this->m_FixedPointSet->GetPoints()->Begin(); fixedIter != this->m_FixedPointSet->GetPoints()->End(); ++fixedIter) {
      const typename FixedPointSetType::PointType fixedPoint = fixedIter.Value();
      const double distance = transformedPoint.SquaredEuclideanDistanceTo(fixedPoint);
      const double expval = exp(-distance / scale);
      const double prbval = valuesOfProbability[fixedIter.Index()];

      for (size_t dim = 0; dim < this->PointDimension; ++dim) {
        gradient[dim] += 2.0 * expval * (transformedPoint[dim] - fixedPoint[dim]) / (scale * prbval);
      }
    }

    // compute derivatives for the current transformed point
    this->m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(this->m_MovingPointSet->GetPoint(transformedIter.Index()), this->m_Jacobian, this->m_JacobianCache);

    for (size_t par = 0; par < this->m_NumberOfParameters; ++par) {
      for (size_t dim = 0; dim < this->PointDimension; ++dim) {
        derivative[par] += this->m_Jacobian(dim, par) * gradient[dim];
      }
    }
  }
}
}

#endif
